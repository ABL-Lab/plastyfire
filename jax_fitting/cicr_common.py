#!/usr/bin/env python3
"""
Shared JAX-accelerated base class and utilities for CICR parameter fitting.

Provides:
  - Batched tensor collation for protocols
  - Pure functional JAX factories for ODE integration
  - Abstract CICRModel base class with vectorized optimization (DE, CMA-ES, NN)
"""

import os, sys, json, time, pickle, logging, argparse
from functools import partial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from pathlib import Path

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

EXPERIMENTAL_TARGETS = {
    "2Hz_5ms":    0.9886,
    "5Hz_5ms":    1.0161,
    "10Hz_10ms":  1.2013,
    "10Hz_-10ms": 0.7922,
    "50Hz_10ms":  1.06,
}

EXPERIMENTAL_ERRORS = {
    "2Hz_5ms":    0.04,
    "5Hz_5ms":    0.0727,
    "10Hz_10ms":  0.0626,
    "10Hz_-10ms": 0.0259,
    "50Hz_10ms":  0.09,
}

BASE_DIR = "/project/rrg-emuller/dhuruva/plastyfitting"
L5_TRACE_DIR = os.path.join(BASE_DIR, "trace_results/Chindemi_params")
L23_TRACE_DIR = os.path.join(BASE_DIR, "trace_results/L23PC_Chindemi_params")
L5_BASIS_DIR = os.path.join(BASE_DIR, "basis_results")
L23_BASIS_DIR = os.path.join(BASE_DIR, "basis_results_L23PC_L5TTPC")

PROTOCOL_PATHWAY = {
    "2Hz_5ms": "L5TTPC", "5Hz_5ms": "L5TTPC",
    "10Hz_10ms": "L5TTPC", "10Hz_-10ms": "L5TTPC",
    "50Hz_10ms": "L23PC", 
}

@partial(jax.jit, static_argnums=(3, 4))
def compute_effcai_piecewise_linear_jax(cai_trace, t, tau_effca=278.318, min_ca=70e-6, effcai0=0.0):
    dt_trace = jnp.diff(t, prepend=t[0])
    
    def scan_step(effcai, inputs):
        f0, f1, dt = inputs
        dt = jnp.where(dt <= 0, 1e-6, dt) 
        a = (f1 - f0) / dt
        decay = jnp.exp(-dt / tau_effca)
        
        term1 = f0 * tau_effca * (1.0 - decay)
        term2 = a * (tau_effca * dt - (tau_effca**2) * (1.0 - decay))
        
        effcai_new = jnp.where(
            inputs[2] <= 0,
            effcai,
            effcai * decay + term1 + term2
        )
        return effcai_new, effcai_new

    f0_arr = cai_trace[:-1] - min_ca
    f1_arr = cai_trace[1:] - min_ca
    dt_arr = dt_trace[1:]

    _, effcai_history = jax.lax.scan(scan_step, effcai0, (f0_arr, f1_arr, dt_arr))
    effcai_full = jnp.concatenate([jnp.array([effcai0]), effcai_history])
    return effcai_full

def _load_pkl(pkl_path):
    with open(pkl_path, "rb") as f: data = pickle.load(f)
    t = np.asarray(data["t"], dtype=np.float64)
    cai = np.asarray(data["cai_CR"], dtype=np.float64)
    if cai.ndim == 1: cai = cai.reshape(1, -1)
    elif cai.shape[0] == len(t) and cai.shape[1] != len(t): cai = cai.T
    
    n_syn = cai.shape[0]
    rho0 = np.zeros(n_syn, dtype=np.float64)
    if "rho_GB" in data:
        rho_gb = np.asarray(data["rho_GB"], dtype=np.float64)
        if rho_gb.ndim == 1: rho_gb = rho_gb.reshape(1, -1)
        elif rho_gb.shape[0] == len(t) and rho_gb.shape[1] != len(t): rho_gb = rho_gb.T
        rho0 = rho_gb[:, 0].copy()
        
    synprops = data.get("synprop", {})
    c_pre = np.asarray(synprops.get("Cpre", np.zeros(n_syn)), dtype=np.float64)
    c_post = np.asarray(synprops.get("Cpost", np.zeros(n_syn)), dtype=np.float64)
    loc_list = synprops.get("loc", ["basal"] * n_syn)
    is_apical = np.array([loc == "apical" for loc in loc_list], dtype=bool)
    
    # Load threshold traces
    pair_name = Path(pkl_path).parent.parent.name
    threshold_path = Path(BASE_DIR) / "trace_results/cpre_cpost_cai_raw_traces" / f"{pair_name}_threshold_traces.pkl"
    if threshold_path.exists():
        with open(threshold_path, "rb") as ft:
            thresh_data = pickle.load(ft)
        # Pad traces into equal matrices (cai_pre, cai_post, t_pre, t_post)
        pre_keys = list(thresh_data["pre"]["cai_CR"].keys())
        # We need to guarantee uniform sequence lengths here per synapse
        t_pre = np.asarray(thresh_data["pre"]["t"], dtype=np.float64)
        t_post = np.asarray(thresh_data["post"]["t"], dtype=np.float64)
        cai_pre = np.zeros((n_syn, len(t_pre)), dtype=np.float64)
        cai_post = np.zeros((n_syn, len(t_post)), dtype=np.float64)
        
        for idx in range(n_syn):
            if idx < len(pre_keys):
                s_key = pre_keys[idx]
                cai_pre[idx] = np.asarray(thresh_data["pre"]["cai_CR"][s_key])
                cai_post[idx] = np.asarray(thresh_data["post"]["cai_CR"][s_key])

        return {"cai": cai, "t": t, "c_pre": c_pre, "c_post": c_post, "is_apical": is_apical, "rho0": rho0, 
                "cai_pre": cai_pre, "t_pre": t_pre, "cai_post": cai_post, "t_post": t_post}
    
    return {"cai": cai, "t": t, "c_pre": c_pre, "c_post": c_post, "is_apical": is_apical, "rho0": rho0}

def _load_basis(pre_gid, post_gid, basis_dir):
    csv_path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    if not os.path.exists(csv_path): return None
    df = pd.read_csv(csv_path)
    configs = df["config"].apply(lambda x: [int(i) for i in x.split(",")])
    n_syn = len(configs.iloc[0])
    baseline_row = df[configs.apply(lambda x: sum(x) == 0)]
    if baseline_row.empty: return None
    
    baseline_mean = float(baseline_row["mean"].values[0])
    singleton_means = np.zeros(n_syn, dtype=np.float64)
    for i in range(n_syn):
        row = df[configs.apply(lambda x: sum(x) == 1 and x[i] == 1)]
        if row.empty: return None
        singleton_means[i] = row["mean"].values[0]
    return {"baseline_mean": baseline_mean, "singleton_means": singleton_means, "n_syn": n_syn}

def preload_all_data(max_pairs=None, protocols=None):
    if protocols is None: protocols = list(EXPERIMENTAL_TARGETS.keys())
    protocol_data = {proto: [] for proto in protocols}

    def _load_dir(trace_dir, basis_dir, protos):
        protos = [p for p in protos if p in protocol_data]
        if not protos or not Path(trace_dir).exists(): return
        dirs = sorted(d for d in Path(trace_dir).iterdir() if d.is_dir())
        for pair_dir in dirs:
            if max_pairs and all(len(protocol_data[p]) >= max_pairs for p in protos if p in protocol_data): break
            parts = pair_dir.name.split("-")
            if len(parts) != 2: continue
            basis = _load_basis(int(parts[0]), int(parts[1]), basis_dir)
            if not basis: continue
            for proto in protos:
                if max_pairs and len(protocol_data[proto]) >= max_pairs: continue
                pkl_path = pair_dir / proto / "simulation_traces.pkl"
                if pkl_path.exists():
                    pd_item = _load_pkl(str(pkl_path))
                    if pd_item["cai"].shape[0] == basis["n_syn"]:
                        pd_item.update(basis)
                        protocol_data[proto].append(pd_item)

    _load_dir(L5_TRACE_DIR, L5_BASIS_DIR, [p for p, pw in PROTOCOL_PATHWAY.items() if pw == "L5TTPC"])
    _load_dir(L23_TRACE_DIR, L23_BASIS_DIR, ["50Hz_10ms"])
    return protocol_data

def _interpolate_pair(pair, interp_dt):
    """Interpolate cai traces to uniform grid with interp_dt (ms) spacing."""
    t_old = pair['t']
    t_new = np.arange(t_old[0], t_old[-1], interp_dt)
    n_syn = pair['cai'].shape[0]
    cai_new = np.zeros((n_syn, len(t_new)))
    for s in range(n_syn):
        cai_new[s] = np.interp(t_new, t_old, pair['cai'][s])
    return {**pair, 'cai': cai_new, 't': t_new}

def collate_protocol_to_jax(pairs_list, dt_step=1, interp_dt=None):
    n_pairs = len(pairs_list)
    if n_pairs == 0: return None
    if interp_dt is not None:
        pairs_list = [_interpolate_pair(p, interp_dt) for p in pairs_list]
    elif dt_step > 1:
        pairs_list = [{**p, 'cai': p['cai'][:, ::dt_step], 't': p['t'][::dt_step]} for p in pairs_list]
    max_time = max(len(p['t']) for p in pairs_list)
    max_syn = max(p['cai'].shape[0] for p in pairs_list)
    
    cai = np.zeros((n_pairs, max_syn, max_time))
    t = np.zeros((n_pairs, max_time))
    c_pre, c_post, rho0, singletons = (np.zeros((n_pairs, max_syn)) for _ in range(4))
    is_apical, valid = (np.zeros((n_pairs, max_syn), dtype=bool) for _ in range(2))
    baseline = np.zeros(n_pairs)
    
    # Pre/Post Traces (assuming pre and post traces have max_time limits, or we use their own full lengths)
    max_t_pre = max(p.get('cai_pre', np.zeros((1,1))).shape[1] for p in pairs_list)
    max_t_post = max(p.get('cai_post', np.zeros((1,1))).shape[1] for p in pairs_list)
    
    # We pad the threshold traces
    cai_pre = np.zeros((n_pairs, max_syn, max_t_pre))
    cai_post = np.zeros((n_pairs, max_syn, max_t_post))
    t_pre = np.zeros((n_pairs, max_t_pre))
    t_post = np.zeros((n_pairs, max_t_post))
    
    for i, p in enumerate(pairs_list):
        ns, nt = p['cai'].shape[0], p['cai'].shape[1]
        cai[i, :ns, :nt] = p['cai']
        if nt < max_time: cai[i, :ns, nt:] = p['cai'][:, -1:]
        t[i, :nt] = p['t']
        if nt < max_time: t[i, nt:] = p['t'][-1]
        c_pre[i, :ns], c_post[i, :ns] = p['c_pre'], p['c_post']
        is_apical[i, :ns], rho0[i, :ns] = p['is_apical'], p['rho0']
        baseline[i], singletons[i, :ns], valid[i, :ns] = p['baseline_mean'], p['singleton_means'], True
        
        # Load threshold traces if they were extracted successfully
        if 'cai_pre' in p:
            nt_pre = p['cai_pre'].shape[1]
            cai_pre[i, :ns, :nt_pre] = p['cai_pre']
            if nt_pre < max_t_pre: cai_pre[i, :ns, nt_pre:] = p['cai_pre'][:, -1:]
            
            t_pre[i, :nt_pre] = p['t_pre']
            if nt_pre < max_t_pre: t_pre[i, nt_pre:] = p['t_pre'][-1]
            
            nt_post = p['cai_post'].shape[1]
            cai_post[i, :ns, :nt_post] = p['cai_post']
            if nt_post < max_t_post: cai_post[i, :ns, nt_post:] = p['cai_post'][:, -1:]
            
            t_post[i, :nt_post] = p['t_post']
            if nt_post < max_t_post: t_post[i, nt_post:] = p['t_post'][-1]
        
    return {
        'cai': jnp.array(cai), 't': jnp.array(t),
        'c_pre': jnp.array(c_pre), 'c_post': jnp.array(c_post),
        'is_apical': jnp.array(is_apical), 'rho0': jnp.array(rho0),
        'baseline': jnp.array(baseline), 'singletons': jnp.array(singletons),
        'valid': jnp.array(valid),
        'cai_pre': jnp.array(cai_pre), 't_pre': jnp.array(t_pre),
        'cai_post': jnp.array(cai_post), 't_post': jnp.array(t_post)
    }

class CICRModel(ABC):
    FIT_PARAMS = []        
    DEFAULT_PARAMS = {}    
    DESCRIPTION = ""       

    def __init__(self):
        self.PARAM_NAMES = [p[0] for p in self.FIT_PARAMS]
        self.PARAM_BOUNDS = [(p[1], p[2]) for p in self.FIT_PARAMS]
        self.LOWER = np.array([b[0] for b in self.PARAM_BOUNDS])
        self.UPPER = np.array([b[1] for b in self.PARAM_BOUNDS])
        self.DEFAULT_X0 = np.array([self.DEFAULT_PARAMS[n] for n in self.PARAM_NAMES])

    @abstractmethod
    def unpack_params(self, x): ...

    @abstractmethod
    def get_step_factory(self): ...

    def get_debug_sim_fn(self): return None

    def get_init_fn(self):
        def default_init(cai_first, rho0): return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rho0)
        return default_init

    def setup_jax(self, protocol_data, targets, lambda_reg=0.0, dt_step=1, interp_dt=None, use_loss_v2=False):
        self.targets_dict = targets
        self.proto_names = list(targets.keys())
        self.target_vals = jnp.array(list(targets.values()))
        self.target_errs = jnp.array([EXPERIMENTAL_ERRORS.get(p, 0.1) for p in self.proto_names])
        
        # Apply specific weights to protocols (e.g. 1.2 for 10Hz_-10ms LTD protocol)
        weights_list = [1.2 if p == "10Hz_-10ms" else 1.0 for p in self.proto_names]
        self.weights = jnp.array(weights_list)
        
        self.lambda_reg = lambda_reg
        
        raw_collated = {p: collate_protocol_to_jax(data, dt_step=dt_step, interp_dt=interp_dt) for p, data in protocol_data.items() if p in targets}
        self.collated_data = {p: d for p, d in raw_collated.items() if d is not None}
        self.proto_names = [p for p in self.proto_names if p in self.collated_data]
        self.target_vals = jnp.array([targets[p] for p in self.proto_names])
        self.weights = jnp.ones(len(self.proto_names))

        step_factory = self.get_step_factory()
        init_fn = self.get_init_fn()

        def sim_synapse(cai_trace, t_trace, c_pre, c_post, is_apical, rho0, sm, bmean, valid, cai_pre, t_pre, cai_post, t_post, params):
            dt_trace = jnp.diff(t_trace, prepend=t_trace[0])
            dt_trace = jnp.where(dt_trace <= 0, 1e-6, dt_trace)
            init = init_fn(cai_trace[0], rho0)
            
            # Pack extended syn_params
            syn_params = (c_pre, c_post, is_apical, cai_pre, t_pre, cai_post, t_post)
            scan_fn = step_factory(params, syn_params)
            
            final, _ = jax.lax.scan(scan_fn, init, (cai_trace, dt_trace))
            return jnp.where(valid & (final[-1] >= 0.5), sm - bmean, 0.0)

        vmap_syn = jax.vmap(sim_synapse, in_axes=(0, None, 0, 0, 0, 0, 0, None, 0, 0, None, 0, None, None))
        
        def sim_pair(cai_p, t_p, cpre_p, cpost_p, isapi_p, rho0_p, bmean, sm_p, valid_p, 
                     caipre_p, tpre_p, caipost_p, tpost_p, params):
                     
            contribs = vmap_syn(cai_p, t_p, cpre_p, cpost_p, isapi_p, rho0_p, sm_p, bmean, valid_p, 
                                caipre_p, tpre_p, caipost_p, tpost_p, params)
                                
            epsp_after = bmean + jnp.sum(contribs)
            contribs_before = jnp.where(valid_p & (rho0_p >= 0.5), sm_p - bmean, 0.0)
            epsp_before = bmean + jnp.sum(contribs_before)
            return jnp.where(epsp_before > 0, epsp_after / epsp_before, jnp.nan)

        vmap_pair = jax.vmap(sim_pair, in_axes=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, None))
        
        def sim_protocol(proto_data, params):
            ratios = vmap_pair(proto_data['cai'], proto_data['t'], proto_data['c_pre'], proto_data['c_post'], 
                               proto_data['is_apical'], proto_data['rho0'], proto_data['baseline'], 
                               proto_data['singletons'], proto_data['valid'], 
                               proto_data['cai_pre'], proto_data['t_pre'], proto_data['cai_post'], proto_data['t_post'],
                               params)
            return jnp.nanmean(ratios)

        def forward_single(x_array, collated):
            params = self.unpack_params(x_array)
            return jnp.stack([sim_protocol(collated[p], params) for p in self.proto_names])

        @jax.jit
        def forward_batch(x_matrix, collated):
            return jax.vmap(forward_single, in_axes=(0, None))(x_matrix, collated)

        idx_ltp = self.proto_names.index("10Hz_10ms") if "10Hz_10ms" in self.proto_names else -1
        idx_ltd = self.proto_names.index("10Hz_-10ms") if "10Hz_-10ms" in self.proto_names else -1

        @jax.jit
        def objective_single(x_array, collated):
            preds = forward_single(x_array, collated)
            
            if use_loss_v2:
                # error = |R_insilico - R_invitro| / Standard_Error(R_invitro)
                base_loss = jnp.sum(self.weights * (jnp.abs(preds - self.target_vals) / self.target_errs))
            else:
                base_loss = jnp.sum(self.weights * (preds - self.target_vals)**2)
                
            sep_penalty = 0.0
            if idx_ltp >= 0 and idx_ltd >= 0:
                actual_sep = preds[idx_ltp] - preds[idx_ltd]
                sep_penalty = 5.0 * jnp.maximum(0.0, 0.41 - actual_sep)**2
            loss = base_loss + sep_penalty
            if self.lambda_reg > 0: loss += self.lambda_reg * jnp.sum((x_array - self.DEFAULT_X0)**2)
            return loss

        @jax.jit
        def objective_batch(x_matrix, collated):
            return jax.vmap(objective_single, in_axes=(0, None))(x_matrix, collated)

        self.forward_batch = partial(forward_batch, collated=self.collated_data)
        self.objective_single = partial(objective_single, collated=self.collated_data)
        self.objective_batch = partial(objective_batch, collated=self.collated_data)

    def run_de(self, max_iter=1000, seed=42, popsize=15, patience=100, **kw):
        from scipy.optimize import differential_evolution
        logger.info(f"Running DE Vectorized (maxiter={max_iter}, popsize={popsize}, patience={patience})")
        def vectorized_obj(x_matrix): return np.array(self.objective_batch(jnp.array(x_matrix.T)))
        stagnation = {"best": np.inf, "count": 0}
        def _callback(intermediate_result):
            f = intermediate_result.fun
            
            # Absolute convergence threshold
            if f < 0.01:
                logger.info(f"Early stopping: Reached target loss ({f:.6f} < 0.01)")
                return True
                
            if stagnation["best"] - f > 1e-8:
                stagnation["best"], stagnation["count"] = f, 0
            else:
                stagnation["count"] += 1
                if stagnation["count"] >= patience:
                    logger.info(f"Early stopping: no improvement for {patience} steps.")
                    return True
            return False

        n_pop = popsize * len(self.PARAM_NAMES)
        rng = np.random.default_rng(seed)
        init_pop = self.LOWER + rng.random((n_pop, len(self.PARAM_NAMES))) * (self.UPPER - self.LOWER)
        init_pop[0] = np.clip(self.DEFAULT_X0, self.LOWER, self.UPPER)

        result = differential_evolution(vectorized_obj, self.PARAM_BOUNDS, maxiter=max_iter, popsize=popsize,
                                        tol=1e-6, seed=seed, disp=True, polish=True, vectorized=True,
                                        callback=_callback, init=init_pop)
        return {"method": "differential_evolution", "x": result.x.tolist(), "fun": float(result.fun)}

    def run_cmaes(self, max_iter=1000, seed=42, popsize=15, patience=100, **kw):
        import cma
        logger.info(f"Running CMA-ES (maxiter={max_iter}, popsize={popsize})")
        
        # Scale to bounds strictly
        sigma_0 = np.mean((self.UPPER - self.LOWER) / 4.0)
        options = {
            'bounds': [self.LOWER.tolist(), self.UPPER.tolist()],
            'maxiter': max_iter,
            'popsize': popsize * len(self.PARAM_NAMES),
            'seed': seed,
            'verbose': -1
        }
        
        x0 = np.clip(self.DEFAULT_X0, self.LOWER, self.UPPER)
        es = cma.CMAEvolutionStrategy(x0, sigma_0, options)
        
        def vectorized_obj(x_list):
            X_mat = jnp.array(x_list)
            # objective_batch expects shape (batch_size, n_params)
            f_jax = self.objective_batch(X_mat)
            return np.array(f_jax).tolist()

        best_f = np.inf
        stag_count = 0
        
        while not es.stop():
            X = es.ask()
            f_vals = vectorized_obj(X)
            es.tell(X, f_vals)
            
            current_best = np.min(f_vals)
            if best_f - current_best > 1e-8:
                best_f = current_best
                stag_count = 0
            else:
                stag_count += 1
                
            es.disp()
            if stag_count >= patience:
                logger.info(f"Early stopping CMA-ES: no improvement for {patience} steps.")
                break
                
        res = es.result
        return {"method": "cmaes", "x": res.xbest.tolist(), "fun": float(res.fbest)}
        
    def run_optax(self, max_iter=1000, seed=42, **kw):
        import optax
        logger.info(f"Running Optax AdamW (maxiter={max_iter})")
        
        # Soften objective wrapper
        # We need a pure scalar function of a 1D array
        idx_ltp = self.proto_names.index("10Hz_10ms") if "10Hz_10ms" in self.proto_names else -1
        idx_ltd = self.proto_names.index("10Hz_-10ms") if "10Hz_-10ms" in self.proto_names else -1

        @jax.jit
        def loss_fn(x_unconstrained):
            # Sigmoid map unconstrained x onto bounds
            lower, upper = jnp.array(self.LOWER), jnp.array(self.UPPER)
            x_val = lower + (upper - lower) * jax.nn.sigmoid(x_unconstrained)
            
            # objective_single natively expects a 1D array of length n_params
            return self.objective_single(x_val)

        optimizer = optax.adamw(learning_rate=0.01)
        
        # Inverse sigmoid for initial params
        x0_val = np.clip(self.DEFAULT_X0, self.LOWER + 1e-6, self.UPPER - 1e-6)
        x_init = np.log((x0_val - self.LOWER) / (self.UPPER - x0_val))
        
        params = jnp.array(x_init)
        opt_state = optimizer.init(params)
        
        value_and_grad_fn = jax.value_and_grad(loss_fn)
        
        best_x = params
        best_f = np.inf
        stag_count = 0
        
        for i in range(max_iter):
            loss_val, grads = value_and_grad_fn(params)
            updates, opt_state = optimizer.update(grads, opt_state, params)
            params = optax.apply_updates(params, updates)
            
            if best_f - loss_val > 1e-8:
                best_f = float(loss_val)
                best_x = params
                stag_count = 0
            else:
                stag_count += 1
                
            if i % 100 == 0:
                logger.info(f"Iter {i}: Loss = {loss_val:.6f}")
                
            if stag_count >= 100:
                logger.info(f"Optax Early Stopping after {i} iterations.")
                break
                
        lower, upper = self.LOWER, self.UPPER
        final_x = lower + (upper - lower) * np.array(jax.nn.sigmoid(best_x))
        return {"method": "optax", "x": final_x.tolist(), "fun": best_f}

    def run_optuna(self, max_iter=1000, seed=42, **kw):
        import optuna
        optuna.logging.set_verbosity(optuna.logging.INFO)
        logger.info(f"Running Optuna TPE (maxiter={max_iter})")
        
        def objective(trial):
            x = []
            for i, name in enumerate(self.PARAM_NAMES):
                x.append(trial.suggest_float(name, self.LOWER[i], self.UPPER[i]))
            loss = float(self.objective_single(jnp.array(x)))
            return loss

        sampler = optuna.samplers.TPESampler(seed=seed)
        study = optuna.create_study(direction='minimize', sampler=sampler)
        
        # Enqueue the default starting guess
        study.enqueue_trial({name: self.DEFAULT_X0[i] for i, name in enumerate(self.PARAM_NAMES)})
        
        study.optimize(objective, n_trials=max_iter)
        
        best = study.best_trial
        best_x = [best.params[name] for name in self.PARAM_NAMES]
        return {"method": "optuna", "x": best_x, "fun": best.value}

    def run_pso(self, max_iter=1000, seed=42, popsize=15, **kw):
        import pyswarms as ps
        logger.info(f"Running PySwarms PSO (maxiter={max_iter}, popsize={popsize})")
        
        options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
        n_particles = popsize * len(self.PARAM_NAMES)
        bounds = (self.LOWER, self.UPPER)
        
        def f_wrapper(x_matrix):
            return np.array(self.objective_batch(jnp.array(x_matrix)))
            
        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles, dimensions=len(self.PARAM_NAMES),
                                            options=options, bounds=bounds)
                                            
        cost, pos = optimizer.optimize(f_wrapper, iters=max_iter)
        return {"method": "pso", "x": pos.tolist(), "fun": float(cost)}

    def print_results(self, opt_result):
        x = np.array(opt_result["x"])
        preds = np.array(self.forward_batch(jnp.array([x])))[0]
        print(f"\n{'='*70}\n  Method: {opt_result['method']}\n  Loss:   {opt_result['fun']:.10f}\n{'='*70}")
        print("\n  Best parameters:")
        for i, name in enumerate(self.PARAM_NAMES):
            d = self.DEFAULT_PARAMS[name]; b = x[i]
            pct = (b - d) / d * 100 if d != 0 else 0
            print(f"    {name:30s} = {b:10.4f}  (default: {d:.4f}, {pct:+.1f}%)")
        print(f"\n  {'Protocol':<15s} {'Predicted':>10s} {'Experiment':>11s} {'Error':>10s}")
        print(f"  {'-'*49}")
        for p, pred, targ in zip(self.proto_names, preds, self.target_vals):
            print(f"  {p:<15s} {pred:10.4f} {targ:11.4f} {pred-targ:+10.4f}")
        print(f"{'='*70}\n")

    def run(self):
        parser = argparse.ArgumentParser(description=f"JAX GB Params ({self.DESCRIPTION})")
        parser.add_argument("--method", choices=["de", "cmaes", "optax", "optuna", "pso"], default="de")
        parser.add_argument("--max-iter", type=int, default=1000)
        parser.add_argument("--protocols", nargs="+", choices=list(EXPERIMENTAL_TARGETS.keys()))
        parser.add_argument("--max-pairs", type=int, default=None)
        parser.add_argument("--dt-step", type=int, default=1)
        parser.add_argument("--popsize", type=int, default=15)
        parser.add_argument("--early-stopping", type=int, default=100, help="Patience for early stopping")
        parser.add_argument("--interp-dt", type=float, default=None, help="Interpolate traces to this dt (ms) for CICR accuracy")
        parser.add_argument("--eval", type=str, default=None, help="JSON file or string of params to evaluate instead of fitting")
        parser.add_argument("--loss-v2", action="store_true", help="Use error = |R_insilico - R_invitro| / StdErr")
        parser.add_argument("--n-plots", type=int, default=1, help="Number of random diagnostic plots to generate")
        args = parser.parse_args()

        targets = {p: v for p, v in EXPERIMENTAL_TARGETS.items() if p in args.protocols} if args.protocols else EXPERIMENTAL_TARGETS
        protocol_data = preload_all_data(max_pairs=args.max_pairs, protocols=list(targets.keys()))
        self.setup_jax(protocol_data, targets, dt_step=args.dt_step, interp_dt=args.interp_dt, use_loss_v2=args.loss_v2)
        
        t0 = time.time()
        ts = time.strftime("%Y%m%d_%H%M%S")
        model_tag = self.DESCRIPTION.lower().split("(")[0].strip().replace(" ", "_")

        if args.eval:
            if os.path.isfile(args.eval):
                with open(args.eval) as f: ep = json.load(f)
            else:
                ep = json.loads(args.eval)
            if "best_parameters" in ep:
                ep = {k: (v["value"] if isinstance(v, dict) else v) for k, v in ep["best_parameters"].items()}
            dp = dict(self.DEFAULT_PARAMS)
            dp.update(ep)
            x_eval = np.array([dp[n] for n in self.PARAM_NAMES])
            loss = float(self.objective_single(jnp.array(x_eval)))
            res = {"method": "eval", "x": x_eval.tolist(), "fun": loss, "time": time.time() - t0}
        else:
            if args.method == "cmaes":
                res = self.run_cmaes(max_iter=args.max_iter, popsize=args.popsize, patience=args.early_stopping)
            elif args.method == "optax":
                res = self.run_optax(max_iter=args.max_iter)
            elif args.method == "optuna":
                res = self.run_optuna(max_iter=args.max_iter)
            elif args.method == "pso":
                res = self.run_pso(max_iter=args.max_iter, popsize=args.popsize)
            else:
                res = self.run_de(max_iter=args.max_iter, popsize=args.popsize, patience=args.early_stopping)
            res["time"] = time.time() - t0
            
            out_json = f"best_params_{model_tag}_{ts}.json"
            best_dict = {n: float(v) for n, v in zip(self.PARAM_NAMES, res["x"])}
            with open(out_json, "w") as f:
                json.dump(best_dict, f, indent=4)
            logger.info(f"Saved best parameters to {out_json}")

        self.print_results(res)
        from plot_cicr_diagnostic import pick_random_pair_syn
        
        # Always generate at least the reference plot #84/3
        plot_out = f"cicr_diagnostic_{model_tag}_{args.method}_{ts}_ref.png"
        logger.info(f"Generating diagnostic plot → {plot_out}")
        self.plot_diagnostic_from_results(pair_idx=84, syn_idx=3, output=plot_out, protocols=["10Hz_10ms", "10Hz_-10ms"], _best_x=res["x"])
        
        # Generate N additional random plots
        for i in range(args.n_plots):
            p_idx, s_idx = pick_random_pair_syn(protocol_data, ["10Hz_10ms", "10Hz_-10ms"])
            rand_out = f"cicr_diagnostic_{model_tag}_{args.method}_{ts}_rand{i+1}_p{p_idx}_s{s_idx}.png"
            logger.info(f"Generating random diagnostic plot {i+1}/{args.n_plots} → {rand_out}")
            self.plot_diagnostic_from_results(pair_idx=p_idx, syn_idx=s_idx, output=rand_out, protocols=["10Hz_10ms", "10Hz_-10ms"], _best_x=res["x"])

    def plot_diagnostic_from_results(self, pair_idx=0, syn_idx=0, output="cicr_diagnostic.png", protocols=("10Hz_10ms", "10Hz_-10ms"), _best_x=None):
        from plot_cicr_diagnostic import compute_single, plot_diagnostic
        dp = dict(self.DEFAULT_PARAMS)
        if _best_x is not None:
            for name, val in zip(self.PARAM_NAMES, _best_x): dp[name] = float(val)
        
        protocol_data = preload_all_data(protocols=list(protocols))
        results = {}
        for proto in protocols:
            pairs = protocol_data.get(proto, [])
            if not pairs: continue
            
            logger.info(f"  Computing {proto} (single synapse {syn_idx} for plot)…")
            res_single = compute_single(pairs[pair_idx], syn_idx, dp, debug_sim_fn=self.get_debug_sim_fn())
            results[proto] = res_single
            
            n_syn = pairs[pair_idx]["cai"].shape[0]
            avg = { 'cai_total': 0, 'cai_raw': 0, 'priming': 0, 'ca_er': 0, 'ca_cicr': 0, 'effcai_no': 0, 'effcai_ci': 0, 'rho_no': 0, 'rho_ci': 0 }
            for s in range(n_syn):
                r = compute_single(pairs[pair_idx], s, dp, debug_sim_fn=self.get_debug_sim_fn())
                avg['cai_total'] += r['cai_total'].max()
                avg['cai_raw'] += r['cai_raw'].max()
                avg['priming'] += r['priming'].max()
                avg['ca_er'] += r['ca_er'].max()
                avg['ca_cicr'] += r['ca_cicr'].max()
                avg['effcai_no'] += r['effcai_no'].max()
                avg['effcai_ci'] += r['effcai_ci'].max()
                avg['rho_no'] += r['rho_no'][-1]
                avg['rho_ci'] += r['rho_ci'][-1]

            def f(k): return f"{avg[k]/n_syn:.6f}"
            print(f"\n  --- {proto} (AVERAGES ACROSS {n_syn} SYNAPSES) ---")
            print(f"  cai_raw max:         {f('cai_raw')} mM")
            print(f"  cai_total max:       {f('cai_total')} mM")
            print(f"  Latent State max:    {f('priming')} au")
            
            # Smart Unit Detection
            er_val = avg['ca_er']/n_syn
            if er_val > 0.1:
                print(f"  ER metric max:       {er_val:.6f} au")
            else:
                print(f"  ca_er max:           {er_val:.6f} mM")
                
            print(f"  ca_cicr max:         {f('ca_cicr')} mM")
            print(f"  effcai max no/ci:    {f('effcai_no')} / {f('effcai_ci')}")
            print(f"  rho final no/ci:     {f('rho_no')} / {f('rho_ci')}")

        if results:
            plot_diagnostic(results, protocols=list(results.keys()), pair_idx=pair_idx, syn_idx=syn_idx, param_label=self.DESCRIPTION, output=output)