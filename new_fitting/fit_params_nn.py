#!/usr/bin/env python3
"""
Neural network surrogate-based parameter fitting.

1. Generate training data: sample random params, run the numba pipeline
   to get EPSP ratios for all pair-protocol combos.
2. Train an MLP surrogate: params (10) -> mean EPSP ratios (5)
3. Optimize over the surrogate with differential_evolution (microsecond evals)
4. Verify the solution on the real forward model

Uses sklearn MLPRegressor (no PyTorch needed).

Usage:
    python fit_params_nn.py --n-samples 50000 --workers 40
    python fit_params_nn.py --data fitting_results/training_data.npz
    python fit_params_nn.py --model fitting_results/surrogate.pkl
"""

import os
import sys
import json
import time
import pickle
import argparse
import logging
import numpy as np
from pathlib import Path
from multiprocessing import Pool

from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error
from scipy.optimize import differential_evolution
from scipy.stats.qmc import Sobol

from fit_params import (
    preload_all_data, _init_objective, objective, evaluate_params,
    _compute_pair_ratio, _compute_rho_final,
    PARAM_NAMES, PARAM_BOUNDS, DEFAULT_PARAMS, DEFAULT_X0,
    EXPERIMENTAL_TARGETS, FIT_PARAMS,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

LOWER = np.array([b[0] for b in PARAM_BOUNDS])
UPPER = np.array([b[1] for b in PARAM_BOUNDS])

ALL_PROTO_NAMES = list(EXPERIMENTAL_TARGETS.keys())
ALL_TARGETS = np.array(list(EXPERIMENTAL_TARGETS.values()))

# These are set in main() based on --protocols selection
PROTO_NAMES = list(ALL_PROTO_NAMES)
TARGETS = ALL_TARGETS.copy()
PROTO_INDICES = list(range(len(ALL_PROTO_NAMES)))  # indices into the 5-column training data

# Module-level data for multiprocessing
_PROTOCOL_DATA = None
# Module-level weights for optimization (set in main)
OPT_WEIGHTS = np.ones(len(PROTO_NAMES))


# ══════════════════════════════════════════════════════════════════
# Data generation
# ══════════════════════════════════════════════════════════════════

def _compute_ratios(x):
    """Forward model: params (18,) -> mean EPSP ratios (N_ALL_PROTOS,).

    Always computes all 5 protocols so training data is reusable.
    Protocol selection happens at training/optimization time.
    """
    gamma_d, gamma_p = x[0], x[1]
    a00, a01, a02, d0 = x[2], x[3], x[4], x[5]
    a10, a11, a12, p0 = x[6], x[7], x[8], x[9]
    a20, a21, a22, d0_ap = x[10], x[11], x[12], x[13]
    a30, a31, a32, p0_ap = x[14], x[15], x[16], x[17]

    ratios = np.full(len(ALL_PROTO_NAMES), np.nan)
    for p_idx, proto in enumerate(ALL_PROTO_NAMES):
        if proto not in _PROTOCOL_DATA:
            continue
        pair_ratios = []
        for pd_item in _PROTOCOL_DATA[proto]:
            r = _compute_pair_ratio(
                pd_item["effcai"], pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                pd_item["is_apical"], pd_item["rho0"],
                pd_item["baseline_mean"], pd_item["singleton_means"],
                gamma_d, gamma_p,
                a00, a01, a02, d0,
                a10, a11, a12, p0,
                a20, a21, a22, d0_ap,
                a30, a31, a32, p0_ap,
            )
            if not np.isnan(r):
                pair_ratios.append(r)
        if pair_ratios:
            ratios[p_idx] = np.mean(pair_ratios)
    return ratios


def generate_training_data(n_samples, workers=1, seed=42):
    """Generate (params, ratios) pairs via Sobol quasi-random sampling."""
    logger.info(f"Sampling {n_samples} parameter vectors (Sobol)...")

    m = int(np.ceil(np.log2(n_samples)))
    n_sobol = 2 ** m

    sampler = Sobol(d=len(PARAM_NAMES), scramble=True, seed=seed)
    X_unit = sampler.random(n_sobol)[:n_samples]
    X = LOWER + X_unit * (UPPER - LOWER)
    X[0] = DEFAULT_X0

    logger.info(f"Running forward model for {n_samples} samples with {workers} workers...")
    t0 = time.time()

    if workers > 1:
        with Pool(workers) as pool:
            Y_list = list(pool.imap(_compute_ratios, X, chunksize=64))
    else:
        Y_list = []
        for i, xi in enumerate(X):
            Y_list.append(_compute_ratios(xi))
            if (i + 1) % 500 == 0:
                logger.info(f"  {i+1}/{n_samples} ({(i+1)/n_samples*100:.0f}%)")

    Y = np.array(Y_list)
    elapsed = time.time() - t0

    valid = ~np.any(np.isnan(Y), axis=1)
    X, Y = X[valid], Y[valid]

    logger.info(f"Generated {len(X)}/{n_samples} valid samples in {elapsed:.1f}s "
                f"({elapsed/n_samples*1000:.1f}ms/sample)")

    for i, proto in enumerate(ALL_PROTO_NAMES):
        logger.info(f"    {proto}: [{Y[:, i].min():.4f}, {Y[:, i].max():.4f}] "
                     f"mean={Y[:, i].mean():.4f}")

    return X, Y


# ══════════════════════════════════════════════════════════════════
# Surrogate training (sklearn)
# ══════════════════════════════════════════════════════════════════

def train_surrogate(X, Y, seed=42):
    """Train sklearn MLPRegressor surrogate."""

    # Normalize inputs to [0, 1]
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_val, Y_train, Y_val = train_test_split(
        X_scaled, Y, test_size=0.15, random_state=seed,
    )

    logger.info(f"Training MLP surrogate: {len(X_train)} train, {len(X_val)} val")

    model = MLPRegressor(
        hidden_layer_sizes=(256, 256, 128),
        activation="relu",
        solver="adam",
        alpha=1e-4,           # L2 regularization (prevents oscillation)
        batch_size=256,
        learning_rate="adaptive",
        learning_rate_init=1e-3,
        max_iter=1000,
        tol=1e-5,             # convergence tolerance
        early_stopping=True,
        validation_fraction=0.15,
        n_iter_no_change=20,  # stop if no val improvement for 20 iters
        random_state=seed,
        verbose=True,
    )

    model.fit(X_train, Y_train)

    # Validation metrics
    Y_pred_val = model.predict(X_val)
    mae_per_proto = np.mean(np.abs(Y_pred_val - Y_val), axis=0)
    mse_val = np.mean((Y_pred_val - Y_val) ** 2)

    logger.info(f"  Converged in {model.n_iter_} iterations")
    logger.info(f"  Validation MSE: {mse_val:.8f}")
    logger.info(f"  Validation MAE per protocol:")
    for i, proto in enumerate(PROTO_NAMES):
        logger.info(f"    {proto}: {mae_per_proto[i]:.6f}")

    return model, scaler


# ══════════════════════════════════════════════════════════════════
# Optimization over surrogate
# ══════════════════════════════════════════════════════════════════

def _run_single_lbfgsb(args):
    """Worker function for parallel L-BFGS-B starts."""
    x0, bounds, model_predict, scaler_transform, weights, targets = args
    from scipy.optimize import minimize

    def nn_objective(x):
        x_scaled = scaler_transform(x.reshape(1, -1))
        pred = np.atleast_1d(model_predict(x_scaled)[0])
        return float(np.sum(weights * (pred - targets) ** 2))

    result = minimize(
        nn_objective, x0, method="L-BFGS-B", bounds=bounds,
        options={"maxiter": 200, "ftol": 1e-15},
    )
    return result.fun, result.x.copy()


def optimize_surrogate(model, scaler, n_starts=500, top_k=20, seed=42, workers=1):
    """
    Optimize over the NN surrogate using multi-start L-BFGS-B.
    Keeps top_k diverse candidates for refinement on the real model.
    Parallelized with multiprocessing when workers > 1.
    """
    import heapq

    logger.info(f"Optimizing over surrogate ({n_starts} L-BFGS-B starts, "
                f"keeping top {top_k}, workers={workers})...")
    t0 = time.time()

    rng = np.random.default_rng(seed)
    x0_list = [DEFAULT_X0.copy()] + [rng.uniform(LOWER, UPPER) for _ in range(n_starts - 1)]

    if workers > 1:
        # Parallel: use joblib for pickling sklearn model
        from joblib import Parallel, delayed
        from scipy.optimize import minimize as _minimize

        def _run_one(x0):
            def nn_obj(x):
                x_scaled = scaler.transform(x.reshape(1, -1))
                pred = np.atleast_1d(model.predict(x_scaled)[0])
                return float(np.sum(OPT_WEIGHTS * (pred - TARGETS) ** 2))
            res = _minimize(nn_obj, x0, method="L-BFGS-B", bounds=PARAM_BOUNDS,
                            options={"maxiter": 200, "ftol": 1e-15})
            return res.fun, res.x.copy()

        results = Parallel(n_jobs=workers, verbose=5)(
            delayed(_run_one)(x0) for x0 in x0_list
        )

        top_candidates = []
        best_loss = np.inf
        for i, (fun, x) in enumerate(results):
            if len(top_candidates) < top_k:
                heapq.heappush(top_candidates, (-fun, x))
            elif fun < -top_candidates[0][0]:
                heapq.heapreplace(top_candidates, (-fun, x))
            if fun < best_loss:
                best_loss = fun
    else:
        # Sequential (original)
        from scipy.optimize import minimize
        top_candidates = []
        best_loss = np.inf

        def nn_objective(x):
            x_scaled = scaler.transform(x.reshape(1, -1))
            pred = np.atleast_1d(model.predict(x_scaled)[0])
            return float(np.sum(OPT_WEIGHTS * (pred - TARGETS) ** 2))

        for i, x0 in enumerate(x0_list):
            result = minimize(
                nn_objective, x0, method="L-BFGS-B", bounds=PARAM_BOUNDS,
                options={"maxiter": 200, "ftol": 1e-15},
            )

            if len(top_candidates) < top_k:
                heapq.heappush(top_candidates, (-result.fun, result.x.copy()))
            elif result.fun < -top_candidates[0][0]:
                heapq.heapreplace(top_candidates, (-result.fun, result.x.copy()))

            if result.fun < best_loss:
                best_loss = result.fun
                logger.info(f"  Start {i+1}/{n_starts}: SSE={result.fun:.12f} (new best)")

    elapsed = time.time() - t0
    logger.info(f"  {n_starts} starts completed in {elapsed:.1f}s")

    # Sort candidates by loss (best first)
    candidates = sorted([(-loss, x) for loss, x in top_candidates])
    best_params = candidates[0][1]
    best_surr_loss = candidates[0][0]

    # NN predictions at best optimum
    x_scaled = scaler.transform(best_params.reshape(1, -1))
    nn_preds = model.predict(x_scaled)[0]
    nn_preds = np.atleast_1d(nn_preds)

    logger.info(f"\n  Best surrogate SSE: {best_surr_loss:.12f}")
    logger.info(f"  Top {top_k} surrogate SSE range: [{candidates[0][0]:.8f}, {candidates[-1][0]:.8f}]")
    for i, proto in enumerate(PROTO_NAMES):
        logger.info(f"    {proto}: nn={nn_preds[i]:.6f} target={TARGETS[i]:.6f} "
                     f"err={nn_preds[i]-TARGETS[i]:+.6f}")

    return best_params, best_surr_loss, nn_preds, [x for _, x in candidates]


def refine_on_real_model(candidates, default_loss, max_time_per_candidate=120):
    """
    Take top-K candidates from NN and evaluate/refine on the real model.

    First picks the best candidate by direct evaluation (fast), then
    refines only the top 3 with Nelder-Mead (slow but gradient-free).
    """
    from scipy.optimize import minimize

    logger.info(f"\nEvaluating {len(candidates)} NN candidates on real model...")

    # Add default params as a candidate
    all_starts = [("default", DEFAULT_X0.copy())] + [
        (f"NN-{i+1}", c) for i, c in enumerate(candidates)
    ]

    # Phase 1: Quick eval of all candidates (just 1 objective call each)
    scored = []
    for label, x0 in all_starts:
        loss = objective(x0)
        scored.append((loss, label, x0))
        logger.info(f"  {label}: real SSE = {loss:.8f}")

    scored.sort(key=lambda t: t[0])
    logger.info(f"\n  Best 3 by direct eval:")
    for loss, label, _ in scored[:3]:
        logger.info(f"    {label}: {loss:.8f}")

    # Phase 2: Nelder-Mead refinement on top 3 only
    n_refine = min(3, len(scored))
    logger.info(f"\nRefining top {n_refine} candidates with Nelder-Mead "
                f"(max {max_time_per_candidate}s each)...")

    best_params = scored[0][2].copy()
    best_loss = scored[0][0]

    for rank in range(n_refine):
        init_loss, label, x0 = scored[rank]
        t0 = time.time()

        # Callback to enforce time limit and log progress
        _eval_count = [0]
        _best_so_far = [init_loss]

        def callback(xk):
            _eval_count[0] += 1
            elapsed = time.time() - t0
            if elapsed > max_time_per_candidate:
                raise StopIteration(f"Time limit {max_time_per_candidate}s reached")
            if _eval_count[0] % 200 == 0:
                curr = objective(xk)
                if curr < _best_so_far[0]:
                    _best_so_far[0] = curr
                logger.info(f"    {label}: {_eval_count[0]} iters, "
                            f"best={_best_so_far[0]:.8f}, {elapsed:.0f}s")

        try:
            result = minimize(
                objective, x0, method="Nelder-Mead", callback=callback,
                options={"maxiter": 2000, "xatol": 1e-6, "fatol": 1e-10, "adaptive": True},
            )
            final_loss = result.fun
            final_x = result.x.copy()
            nfev = result.nfev
        except StopIteration:
            # Time limit hit — evaluate current best
            final_loss = _best_so_far[0]
            final_x = x0  # Nelder-Mead doesn't expose intermediate x easily
            nfev = _eval_count[0]

        elapsed = time.time() - t0
        logger.info(f"  {label}: {init_loss:.8f} -> {final_loss:.8f} "
                     f"({nfev} evals, {elapsed:.0f}s)")

        if final_loss < best_loss:
            best_loss = final_loss
            best_params = final_x
            logger.info(f"    ^ new best real SSE: {best_loss:.10f}")

    logger.info(f"\n  Best real SSE after refinement: {best_loss:.10f} (default was {default_loss:.10f})")
    return best_params, best_loss


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="NN surrogate fitting for plasticity parameters (sklearn)",
    )
    parser.add_argument("--n-samples", type=int, default=50000)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--data", type=str, default=None,
                        help="Skip generation, load .npz")
    parser.add_argument("--model", type=str, default=None,
                        help="Skip training, load .pkl")
    parser.add_argument("--n-starts", type=int, default=500,
                        help="Number of L-BFGS-B random starts over surrogate")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save-dir", type=str, default="fitting_results")
    parser.add_argument("--weights", type=float, nargs="+", default=None,
                        help="Per-protocol weights for loss (must match number of selected protocols)")
    parser.add_argument("--protocols", type=str, nargs="+", default=None,
                        help=f"Which protocols to fit. Default: all. Available: {ALL_PROTO_NAMES}")

    args = parser.parse_args()
    save_dir = Path(args.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ── Set up protocol selection ──
    global PROTO_NAMES, TARGETS, PROTO_INDICES, OPT_WEIGHTS

    if args.protocols is not None:
        for p in args.protocols:
            if p not in ALL_PROTO_NAMES:
                parser.error(f"Unknown protocol '{p}'. Available: {ALL_PROTO_NAMES}")
        PROTO_INDICES = [ALL_PROTO_NAMES.index(p) for p in args.protocols]
        PROTO_NAMES = [ALL_PROTO_NAMES[i] for i in PROTO_INDICES]
        TARGETS = ALL_TARGETS[PROTO_INDICES]
    else:
        PROTO_INDICES = list(range(len(ALL_PROTO_NAMES)))
        PROTO_NAMES = list(ALL_PROTO_NAMES)
        TARGETS = ALL_TARGETS.copy()

    logger.info(f"Selected protocols: {PROTO_NAMES}")

    # ── Set up weights ──
    OPT_WEIGHTS = np.ones(len(PROTO_NAMES))
    if args.weights is not None:
        assert len(args.weights) == len(PROTO_NAMES), \
            f"Need {len(PROTO_NAMES)} weights (for {PROTO_NAMES}), got {len(args.weights)}"
        OPT_WEIGHTS = np.array(args.weights)
    logger.info(f"Protocol weights: {dict(zip(PROTO_NAMES, OPT_WEIGHTS))}")

    # Build the target dict for the selected protocols only
    selected_targets = {p: EXPERIMENTAL_TARGETS[p] for p in PROTO_NAMES}

    # ── 1. Preload traces + basis ──
    logger.info("Preloading trace and basis data...")
    global _PROTOCOL_DATA
    protocol_data = preload_all_data()
    _PROTOCOL_DATA = protocol_data
    _init_objective(protocol_data, selected_targets, weights=OPT_WEIGHTS)

    logger.info("Warming up numba...")
    _ = _compute_ratios(DEFAULT_X0)
    default_loss = objective(DEFAULT_X0)
    logger.info(f"Default params real loss: {default_loss:.10f}")

    # ── 2. Generate or load training data ──
    if args.data:
        logger.info(f"Loading data from {args.data}")
        npz = np.load(args.data)
        X, Y_all = npz["X"], npz["Y"]
        logger.info(f"Loaded {len(X)} samples (all {Y_all.shape[1]} protocols)")
        # Select only the columns for chosen protocols
        Y = Y_all[:, PROTO_INDICES]
        # Re-filter for NaN in selected columns only
        valid = ~np.any(np.isnan(Y), axis=1)
        X, Y = X[valid], Y[valid]
        logger.info(f"After selecting {PROTO_NAMES}: {len(X)} valid samples")
    else:
        X, Y_all = generate_training_data(args.n_samples, args.workers, args.seed)
        data_path = save_dir / "training_data.npz"
        np.savez(data_path, X=X, Y=Y_all)
        logger.info(f"Saved training data to {data_path}")
        # Select columns for chosen protocols
        Y = Y_all[:, PROTO_INDICES]
        valid = ~np.any(np.isnan(Y), axis=1)
        X, Y = X[valid], Y[valid]
        logger.info(f"After selecting {PROTO_NAMES}: {len(X)} valid samples")

    # ── 3. Train or load surrogate ──
    if args.model:
        logger.info(f"Loading model from {args.model}")
        with open(args.model, "rb") as f:
            saved = pickle.load(f)
        model, scaler = saved["model"], saved["scaler"]
        # Verify output dimension matches
        test_pred = model.predict(scaler.transform(DEFAULT_X0.reshape(1, -1)))
        n_out = test_pred.shape[1] if test_pred.ndim == 2 else 1
        if n_out != len(PROTO_NAMES):
            logger.error(f"Loaded model outputs {n_out} protocols "
                         f"but {len(PROTO_NAMES)} selected. Retrain without --model.")
            sys.exit(1)
    else:
        model, scaler = train_surrogate(X, Y, seed=args.seed)
        model_path = save_dir / "surrogate.pkl"
        with open(model_path, "wb") as f:
            pickle.dump({"model": model, "scaler": scaler,
                         "protocols": PROTO_NAMES}, f)
        logger.info(f"Saved model to {model_path}")

    # ── 4. Optimize over surrogate ──
    t0 = time.time()
    best_params, surr_loss, nn_preds, top_candidates = optimize_surrogate(
        model, scaler, n_starts=args.n_starts, seed=args.seed,
        workers=args.workers,
    )
    opt_time = time.time() - t0

    # ── 5. Refine top candidates on real forward model ──
    best_params, real_loss = refine_on_real_model(top_candidates, default_loss)
    real_eval = evaluate_params(best_params)

    print(f"\n{'='*70}")
    print(f"  Neural Network Surrogate Results (sklearn)")
    print(f"{'='*70}")
    print(f"  Protocols:         {PROTO_NAMES}")
    print(f"  Training samples:  {len(X)}")
    print(f"  Surrogate SSE:     {surr_loss:.10f}")
    print(f"  Real SSE:          {real_loss:.10f}")
    print(f"  Default SSE:       {default_loss:.10f}")
    print(f"  Surrogate opt:     {opt_time:.1f}s")
    print(f"{'='*70}")

    print("\n  Best parameters:")
    for i, name in enumerate(PARAM_NAMES):
        d = DEFAULT_PARAMS[name]
        b = best_params[i]
        pct = (b - d) / d * 100 if d != 0 else 0
        print(f"    {name:30s} = {b:10.4f}  (default: {d:.4f}, {pct:+.1f}%)")

    print(f"\n  {'Protocol':<15s} {'NN':>8s} {'Real':>8s} {'Exper':>8s} {'RealErr':>9s}")
    print(f"  {'-'*48}")
    for i, proto in enumerate(PROTO_NAMES):
        r = real_eval[proto]
        print(f"  {proto:<15s} {nn_preds[i]:8.4f} {r['predicted']:8.4f} "
              f"{r['experimental']:8.4f} {r['error']:+9.4f}")
    print(f"{'='*70}\n")

    # ── 6. Save ──
    results = {
        "method": "neural_network_surrogate",
        "x": best_params.tolist(),
        "fun_surrogate": float(surr_loss),
        "fun_real": float(real_loss),
        "fun_default": float(default_loss),
        "nn_predictions": nn_preds.tolist(),
        "per_protocol": real_eval,
        "param_names": PARAM_NAMES,
        "default_params": DEFAULT_PARAMS,
        "n_training_samples": len(X),
        "opt_time_seconds": opt_time,
        "selected_protocols": PROTO_NAMES,
        "weights": OPT_WEIGHTS.tolist(),
    }
    result_path = save_dir / "result_nn.json"
    with open(result_path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Results saved to {result_path}")


if __name__ == "__main__":
    main()
