"""
Pipeline: effcai traces -> rho prediction -> EPSP values.

Takes effective calcium traces (effcai_GB) as input and produces predicted EPSP
amplitudes as output by:
    1. Computing rho from effcai using the Graupner-Brunel ODE
    2. Binarizing rho (>= 0.5 -> potentiated, < 0.5 -> depressed)
    3. Extrapolating EPSP via linear superposition from basis results

Linear superposition (from verify_linearity.py):
    EPSP(config) = baseline + sum_i[ (singleton_i - baseline) ]  for each potentiated synapse i

where:
    baseline     = EPSP when all synapses are depressed  (config 0,0,...,0)
    singleton_i  = EPSP when only synapse i is potentiated (config 0,..,1_i,..,0)

Trial-to-trial variability:
    Basis results store mean and std across 10 trials (stochastic vesicle release).
    Variance is propagated through linear superposition:
        Var(EPSP) = (1-k)^2 * var_baseline + sum_i[ var_singleton_i ]
    where k = number of potentiated synapses.
    When n_trials > 0, samples are drawn from N(mean, std^2).
"""

import os
import argparse
import pickle
import logging
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ── Graupner-Brunel model parameters ────────────────────────────────────────
TAU_IND_GB = 70.0      # s
RHO_STAR_GB = 0.5

DEFAULT_PARAMS = {
    "gamma_d_GB_GluSynapse": 101.5,
    "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002, "a01": 1.954,
    "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456,
    "a30": 5.236, "a31": 1.782,
}

# Where per-pair basis CSVs live  (basis_{pre}_{post}.csv)
BASIS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/basis_results"


# ═══════════════════════════════════════════════════════════════════════════════
# Basis loading and linear extrapolation
# ═══════════════════════════════════════════════════════════════════════════════

def load_basis(pre_gid, post_gid, basis_dir=None):
    """
    Load basis results for a pre-post pair.

    Returns
    -------
    baseline_mean : float
        Mean EPSP amplitude when all synapses are depressed.
    baseline_std : float
        Std of EPSP amplitude when all synapses are depressed.
    singleton_means : list[float]
        Mean EPSP when only synapse i is potentiated (length = n_synapses).
    singleton_stds : list[float]
        Std of EPSP when only synapse i is potentiated (length = n_synapses).
    n_synapses : int
    """
    if basis_dir is None:
        basis_dir = BASIS_DIR

    csv_path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"No basis file for {pre_gid}->{post_gid} at {csv_path}"
        )

    df = pd.read_csv(csv_path)
    df["config_list"] = df["config"].apply(lambda x: [int(i) for i in x.split(",")])
    n_synapses = len(df["config_list"].iloc[0])

    # Baseline: all zeros
    baseline_row = df[df["config_list"].apply(lambda x: sum(x) == 0)]
    if baseline_row.empty:
        raise ValueError(f"Baseline config (all 0s) missing in {csv_path}")
    baseline_mean = baseline_row["mean"].values[0]
    baseline_std = baseline_row["std"].values[0]

    # Singletons: exactly one 1 at position i
    singleton_means = []
    singleton_stds = []
    for i in range(n_synapses):
        row = df[df["config_list"].apply(lambda x: sum(x) == 1 and x[i] == 1)]
        if row.empty:
            raise ValueError(f"Singleton config for synapse {i} missing in {csv_path}")
        singleton_means.append(row["mean"].values[0])
        singleton_stds.append(row["std"].values[0])

    return baseline_mean, baseline_std, singleton_means, singleton_stds, n_synapses


def extrapolate_epsp(rho_config, baseline_mean, baseline_std,
                     singleton_means, singleton_stds, n_trials=0, rng=None):
    """
    Predict EPSP for an arbitrary rho configuration using linear superposition.

    Mean:
        EPSP = baseline + sum_i[ (singleton_i - baseline) ]  for each i where rho_i == 1

    Variance propagation (independent basis simulations):
        EPSP = B*(1-k) + sum_{active} S_i
        Var(EPSP) = (1-k)^2 * Var(B) + sum_{active} Var(S_i)

    Parameters
    ----------
    rho_config : list[int] or str
    baseline_mean, baseline_std : float
    singleton_means, singleton_stds : list[float]
    n_trials : int
        If > 0, sample n_trials values from N(mean, std^2).
        If 0, return deterministic mean only.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    epsp_mean : float
    epsp_std : float
    trials : np.ndarray or None
        Sampled EPSP values if n_trials > 0, else None.
    """
    if isinstance(rho_config, str):
        rho_config = [int(x) for x in rho_config.split(",")]

    # Mean: baseline + sum of deltas
    epsp_mean = baseline_mean
    for i, rho_val in enumerate(rho_config):
        if rho_val >= 0.5:
            epsp_mean += (singleton_means[i] - baseline_mean)

    # Std propagation: EPSP = B*(1-k) + sum_{active} S_i
    k = sum(1 for r in rho_config if r >= 0.5)
    variance = (1 - k) ** 2 * baseline_std ** 2
    for i, rho_val in enumerate(rho_config):
        if rho_val >= 0.5:
            variance += singleton_stds[i] ** 2
    epsp_std = np.sqrt(variance)

    # Sample trials
    trials = None
    if n_trials > 0:
        if rng is None:
            rng = np.random.default_rng()
        trials = rng.normal(epsp_mean, epsp_std, size=n_trials)

    return epsp_mean, epsp_std, trials


# ═══════════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════════

def fetch_epsp_value(pre_gid, post_gid, freq, pre_post_displacement, rho_config,
                     basis_dir=None, n_trials=0, seed=None):
    """
    Predict EPSP amplitude for a given pair and rho configuration
    using linear extrapolation from basis results.

    Parameters
    ----------
    pre_gid : int
    post_gid : int
    freq : float
        Stimulation frequency in Hz (for context / logging).
    pre_post_displacement : float
        Pre-post spike timing displacement in ms (for context / logging).
    rho_config : list[int] or str
        Binary rho state per synapse, e.g. [0,1,1,0] or "0,1,1,0".
    basis_dir : str, optional
    n_trials : int
        If > 0, sample n_trials EPSP values from N(mean, std^2).
    seed : int, optional
        Random seed for trial sampling.

    Returns
    -------
    dict with keys:
        epsp_mean : float
        epsp_std : float
        trials : np.ndarray or None
    """
    baseline_mean, baseline_std, singleton_means, singleton_stds, n_syn = \
        load_basis(pre_gid, post_gid, basis_dir)

    rng = np.random.default_rng(seed) if seed is not None else np.random.default_rng()
    epsp_mean, epsp_std, trials = extrapolate_epsp(
        rho_config, baseline_mean, baseline_std,
        singleton_means, singleton_stds, n_trials, rng
    )

    logger.info(
        f"EPSP for {pre_gid}->{post_gid} ({freq}Hz, {pre_post_displacement}ms) "
        f"rho={rho_config}: mean={epsp_mean:.4f} +/- {epsp_std:.4f} mV"
    )
    if trials is not None:
        logger.info(f"  sampled {n_trials} trials: "
                     f"mean={np.mean(trials):.4f}, std={np.std(trials):.4f}")

    return {"epsp_mean": epsp_mean, "epsp_std": epsp_std, "trials": trials}


def fetch_epsp_ratio(pre_gid, post_gid, freq, pre_post_displacement,
                     initial_rho_config, final_rho_config, basis_dir=None,
                     n_trials=0, seed=None):
    """
    Predict EPSP ratio (after / before) for a given pair using linear
    extrapolation from basis results.

    Parameters
    ----------
    pre_gid, post_gid : int
    freq : float
    pre_post_displacement : float
    initial_rho_config : list[int] or str
    final_rho_config : list[int] or str
    basis_dir : str, optional
    n_trials : int
        If > 0, sample trial pairs and compute per-trial ratios.
    seed : int, optional

    Returns
    -------
    dict with keys:
        ratio_mean : float  (from mean EPSPs)
        epsp_before_mean, epsp_before_std : float
        epsp_after_mean, epsp_after_std : float
        trial_ratios : np.ndarray or None
    """
    baseline_mean, baseline_std, singleton_means, singleton_stds, _ = \
        load_basis(pre_gid, post_gid, basis_dir)

    rng = np.random.default_rng(seed) if seed is not None else np.random.default_rng()

    before_mean, before_std, before_trials = extrapolate_epsp(
        initial_rho_config, baseline_mean, baseline_std,
        singleton_means, singleton_stds, n_trials, rng
    )
    after_mean, after_std, after_trials = extrapolate_epsp(
        final_rho_config, baseline_mean, baseline_std,
        singleton_means, singleton_stds, n_trials, rng
    )

    ratio_mean = after_mean / before_mean if before_mean != 0 else float("nan")

    trial_ratios = None
    if before_trials is not None and after_trials is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            trial_ratios = np.where(before_trials != 0,
                                    after_trials / before_trials, np.nan)

    logger.info(
        f"EPSP ratio {pre_gid}->{post_gid} ({freq}Hz, {pre_post_displacement}ms): "
        f"{before_mean:.4f} +/- {before_std:.4f} -> "
        f"{after_mean:.4f} +/- {after_std:.4f} = {ratio_mean:.4f}"
    )
    if trial_ratios is not None:
        logger.info(f"  sampled {n_trials} trial ratios: "
                     f"mean={np.nanmean(trial_ratios):.4f}, "
                     f"std={np.nanstd(trial_ratios):.4f}")

    return {
        "ratio_mean": ratio_mean,
        "epsp_before_mean": before_mean,
        "epsp_before_std": before_std,
        "epsp_after_mean": after_mean,
        "epsp_after_std": after_std,
        "trial_ratios": trial_ratios,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Graupner-Brunel: effcai -> rho
# ═══════════════════════════════════════════════════════════════════════════════

def calculate_thresholds(c_pre, c_post, loc, params):
    """Calculate theta_d and theta_p from synapse location and calcium levels."""
    if loc == "basal":
        theta_d = params["a00"] * c_pre + params["a01"] * c_post
        theta_p = params["a10"] * c_pre + params["a11"] * c_post
    elif loc == "apical":
        theta_d = params["a20"] * c_pre + params["a21"] * c_post
        theta_p = params["a30"] * c_pre + params["a31"] * c_post
    else:
        raise ValueError(f"Unknown location: {loc}")
    return theta_d, theta_p


def compute_rho(effcai_trace, t, theta_d, theta_p, params, rho0=0.0):
    """
    Integrate the Graupner-Brunel rho ODE (Euler).

    rho' = (-rho*(1-rho)*(rho_star - rho)
            + pot*gamma_p*(1-rho)
            - dep*gamma_d*rho) / (tau_ind * 1e3)
    """
    n = len(t)
    rho = np.zeros(n)
    rho[0] = rho0

    gamma_d = params["gamma_d_GB_GluSynapse"]
    gamma_p = params["gamma_p_GB_GluSynapse"]
    inv_tau = 1.0 / (TAU_IND_GB * 1000.0)

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        r = rho[i]
        dep = 1.0 if effcai_trace[i] > theta_d else 0.0
        pot = 1.0 if effcai_trace[i] > theta_p else 0.0

        drho = (-r * (1.0 - r) * (RHO_STAR_GB - r)
                + pot * gamma_p * (1.0 - r)
                - dep * gamma_d * r) * inv_tau

        rho[i + 1] = np.clip(r + dt * drho, 0.0, 1.0)

    return rho


def binarize_rho(rho_values):
    """Binarize: >= 0.5 -> 1, < 0.5 -> 0."""
    return (np.asarray(rho_values) >= 0.5).astype(int)


# ═══════════════════════════════════════════════════════════════════════════════
# Full pipeline: effcai traces -> rho -> extrapolated EPSP
# ═══════════════════════════════════════════════════════════════════════════════

def effcai_to_epsp(effcai_traces, t, synprops, pre_gid, post_gid,
                   freq, pre_post_displacement, params=None, basis_dir=None,
                   n_trials=0, seed=None):
    """
    Full pipeline: effcai traces -> rho -> EPSP via linear extrapolation.

    Parameters
    ----------
    effcai_traces : np.ndarray, shape (n_synapses, n_timepoints)
    t : np.ndarray, shape (n_timepoints,)   — time in ms
    synprops : dict  — keys: "Cpre", "Cpost", "loc", "rho0_GB"
    pre_gid, post_gid : int
    freq : float
    pre_post_displacement : float
    params : dict, optional
    basis_dir : str, optional
    n_trials : int
        If > 0, sample trial EPSP values from propagated distribution.
    seed : int, optional

    Returns
    -------
    dict with rho_traces, rho_binary_initial, rho_binary_final,
         epsp_before_mean, epsp_before_std, epsp_after_mean, epsp_after_std,
         epsp_ratio_mean, before_trials, after_trials, ratio_trials
    """
    if params is None:
        params = DEFAULT_PARAMS

    n_synapses = effcai_traces.shape[0]

    # 1) effcai -> rho
    rho_traces = np.zeros_like(effcai_traces)
    for i in range(n_synapses):
        c_pre = synprops["Cpre"][i]
        c_post = synprops["Cpost"][i]
        loc = synprops["loc"][i]
        rho0 = synprops.get("rho0_GB", [0.0] * n_synapses)[i]

        theta_d, theta_p = calculate_thresholds(c_pre, c_post, loc, params)
        rho_traces[i] = compute_rho(effcai_traces[i], t, theta_d, theta_p, params, rho0)

        logger.info(
            f"  syn {i}: loc={loc}, theta_d={theta_d:.4f}, theta_p={theta_p:.4f}, "
            f"rho {rho_traces[i, 0]:.3f} -> {rho_traces[i, -1]:.3f}"
        )

    # 2) binarize
    rho_binary_initial = binarize_rho(rho_traces[:, 0])
    rho_binary_final = binarize_rho(rho_traces[:, -1])

    # 3) extrapolate EPSP
    baseline_mean, baseline_std, singleton_means, singleton_stds, _ = \
        load_basis(pre_gid, post_gid, basis_dir)

    rng = np.random.default_rng(seed) if seed is not None else np.random.default_rng()

    before_mean, before_std, before_trials = extrapolate_epsp(
        rho_binary_initial, baseline_mean, baseline_std,
        singleton_means, singleton_stds, n_trials, rng
    )
    after_mean, after_std, after_trials = extrapolate_epsp(
        rho_binary_final, baseline_mean, baseline_std,
        singleton_means, singleton_stds, n_trials, rng
    )

    ratio_mean = after_mean / before_mean if before_mean != 0 else float("nan")

    ratio_trials = None
    if before_trials is not None and after_trials is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_trials = np.where(before_trials != 0,
                                    after_trials / before_trials, np.nan)

    logger.info(f"  initial rho: {rho_binary_initial.tolist()}")
    logger.info(f"  final   rho: {rho_binary_final.tolist()}")
    logger.info(f"  EPSP before: {before_mean:.4f} +/- {before_std:.4f} mV")
    logger.info(f"  EPSP after:  {after_mean:.4f} +/- {after_std:.4f} mV")
    logger.info(f"  EPSP ratio:  {ratio_mean:.4f}")

    return {
        "rho_traces": rho_traces,
        "rho_binary_initial": rho_binary_initial,
        "rho_binary_final": rho_binary_final,
        "epsp_before_mean": before_mean,
        "epsp_before_std": before_std,
        "epsp_after_mean": after_mean,
        "epsp_after_std": after_std,
        "epsp_ratio_mean": ratio_mean,
        "before_trials": before_trials,
        "after_trials": after_trials,
        "ratio_trials": ratio_trials,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Data loading helper
# ═══════════════════════════════════════════════════════════════════════════════

def load_simulation_traces(pkl_path):
    """Load effcai traces and synapse properties from simulation_traces.pkl."""
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)

    t = data.get("t", np.arange(np.asarray(data["effcai_GB"]).shape[-1]))
    effcai = np.asarray(data["effcai_GB"])

    if effcai.shape[0] == len(t) and (effcai.ndim == 1 or effcai.shape[1] != len(t)):
        effcai = effcai.T if effcai.ndim > 1 else effcai.reshape(1, -1)

    synprops = data.get("synprop", {})
    n_syn = effcai.shape[0]

    defaults = {
        "Cpre": np.zeros(n_syn),
        "Cpost": np.zeros(n_syn),
        "loc": ["basal"] * n_syn,
        "rho0_GB": np.zeros(n_syn),
    }
    for key, default in defaults.items():
        if key not in synprops or synprops[key] is None or len(synprops[key]) == 0:
            logger.warning(f"synprop['{key}'] missing, using default")
            synprops[key] = default

    return effcai, t, synprops


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Predict EPSP values from effcai traces via linear extrapolation"
    )
    subparsers = parser.add_subparsers(dest="command")

    # ── fetch ────────────────────────────────────────────────────────────
    p = subparsers.add_parser("fetch", help="Extrapolate EPSP for a rho config")
    p.add_argument("pre_gid", type=int)
    p.add_argument("post_gid", type=int)
    p.add_argument("freq", type=float)
    p.add_argument("dt", type=float)
    p.add_argument("rho_config", help="e.g. 0,1,1,0")
    p.add_argument("--basis-dir", default=None)
    p.add_argument("--n-trials", type=int, default=0)
    p.add_argument("--seed", type=int, default=None)

    # ── ratio ────────────────────────────────────────────────────────────
    p = subparsers.add_parser("ratio", help="Extrapolate EPSP ratio")
    p.add_argument("pre_gid", type=int)
    p.add_argument("post_gid", type=int)
    p.add_argument("freq", type=float)
    p.add_argument("dt", type=float)
    p.add_argument("initial_rho", help="e.g. 0,0,1,0")
    p.add_argument("final_rho", help="e.g. 1,0,1,0")
    p.add_argument("--basis-dir", default=None)
    p.add_argument("--n-trials", type=int, default=0)
    p.add_argument("--seed", type=int, default=None)

    # ── pipeline ─────────────────────────────────────────────────────────
    p = subparsers.add_parser("pipeline", help="Full effcai->rho->EPSP pipeline")
    p.add_argument("pkl_path")
    p.add_argument("pre_gid", type=int)
    p.add_argument("post_gid", type=int)
    p.add_argument("freq", type=float)
    p.add_argument("dt", type=float)
    p.add_argument("--basis-dir", default=None)
    p.add_argument("--n-trials", type=int, default=0)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--output", default=None)
    for name in DEFAULT_PARAMS:
        p.add_argument(f"--{name}", type=float, default=None)

    args = parser.parse_args()

    if args.command == "fetch":
        result = fetch_epsp_value(args.pre_gid, args.post_gid,
                                  args.freq, args.dt, args.rho_config,
                                  args.basis_dir, args.n_trials, args.seed)
        print(f"EPSP = {result['epsp_mean']:.4f} +/- {result['epsp_std']:.4f} mV")
        if result["trials"] is not None:
            print(f"Trials ({args.n_trials}): {result['trials']}")

    elif args.command == "ratio":
        result = fetch_epsp_ratio(args.pre_gid, args.post_gid,
                                  args.freq, args.dt,
                                  args.initial_rho, args.final_rho,
                                  args.basis_dir, args.n_trials, args.seed)
        print(f"EPSP ratio = {result['ratio_mean']:.4f}")
        print(f"  before: {result['epsp_before_mean']:.4f} +/- {result['epsp_before_std']:.4f} mV")
        print(f"  after:  {result['epsp_after_mean']:.4f} +/- {result['epsp_after_std']:.4f} mV")
        if result["trial_ratios"] is not None:
            print(f"Trial ratios ({args.n_trials}): "
                  f"mean={np.nanmean(result['trial_ratios']):.4f}, "
                  f"std={np.nanstd(result['trial_ratios']):.4f}")

    elif args.command == "pipeline":
        effcai, t, synprops = load_simulation_traces(args.pkl_path)
        logger.info(f"Loaded {effcai.shape[0]} synapses, {effcai.shape[1]} timepoints")

        params = dict(DEFAULT_PARAMS)
        for name in DEFAULT_PARAMS:
            val = getattr(args, name, None)
            if val is not None:
                params[name] = val

        results = effcai_to_epsp(effcai, t, synprops,
                                 args.pre_gid, args.post_gid,
                                 args.freq, args.dt, params,
                                 args.basis_dir, args.n_trials, args.seed)

        print(f"\n{'='*60}")
        print(f"  Pair:         {args.pre_gid} -> {args.post_gid}")
        print(f"  Protocol:     {args.freq}Hz, {args.dt}ms")
        print(f"  Synapses:     {effcai.shape[0]}")
        print(f"  Initial rho:  {results['rho_binary_initial'].tolist()}")
        print(f"  Final rho:    {results['rho_binary_final'].tolist()}")
        print(f"  EPSP before:  {results['epsp_before_mean']:.4f} +/- {results['epsp_before_std']:.4f} mV")
        print(f"  EPSP after:   {results['epsp_after_mean']:.4f} +/- {results['epsp_after_std']:.4f} mV")
        print(f"  EPSP ratio:   {results['epsp_ratio_mean']:.4f}")
        if results["ratio_trials"] is not None:
            print(f"  Trial ratios: mean={np.nanmean(results['ratio_trials']):.4f}, "
                  f"std={np.nanstd(results['ratio_trials']):.4f}")
        print(f"{'='*60}")

        if args.output:
            with open(args.output, "wb") as f:
                pickle.dump(results, f)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
