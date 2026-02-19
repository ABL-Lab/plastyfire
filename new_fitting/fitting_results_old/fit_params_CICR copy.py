#!/usr/bin/env python3
"""
Parameter fitting for Graupner-Brunel plasticity model.

Optimizes 10 parameters (gamma_d, gamma_p, a00-a31) to minimize
the sum of squared errors between in-silico and experimental EPSP ratios
across 5 STDP protocols.

Pipeline per objective evaluation:
    candidate params -> thresholds -> rho ODE (numba) -> binarize
    -> EPSP extrapolation from basis -> ratio -> loss

Usage:
    python fit_params.py --method de --max-iter 1000 --workers 4
    python fit_params.py --method cmaes --max-iter 1000
    python fit_params.py --method optuna --n-trials 5000
    python fit_params.py --method lbfgsb --n-starts 50
    python fit_params.py --method all     # run all 4 and compare
    python fit_params.py --generate-slurm  # create SLURM submission script
"""

import os
import sys
import json
import time
import pickle
import logging
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from numba import njit
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ══════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════

TAU_IND_GB = 70.0       # seconds
RHO_STAR_GB = 0.5
TAU_EFFCA_GB = 200.0    # ms, effective calcium time constant (GluSynapse.mod default)
MIN_CA_CR = 70e-6       # mM (70 nM), baseline calcium concentration

# CICR constants (from cawave.py / cawave.cfg)
K_SERCA_CICR = 0.1      # uM, SERCA half-activation (Michaelis constant)
FC_CICR = 0.83           # cytoplasmic volume fraction
FE_CICR = 0.17           # ER volume fraction
CA_AVG_CICR = 0.0017     # mM, conserved total calcium (fc*cai + fe*ca_er)

LOWER = None  # set after PARAM_BOUNDS is defined
UPPER = None

EXPERIMENTAL_TARGETS = {
    "2Hz_5ms":    0.9886,   # mrk97_01
    "5Hz_5ms":    1.0161,   # mrk97_02
    "10Hz_10ms":  1.2013,   # mrk97_07
    "10Hz_-10ms": 0.7922,   # mrk97_08
    "50Hz_10ms":  1.06,     # sjh06_02
}

FIT_PARAMS = [
    # gamma
    ("gamma_d_GB_GluSynapse", 50.0, 200.0),
    ("gamma_p_GB_GluSynapse", 150.0, 300.0),
    # basal theta_d: a00*c_pre + a01*c_post
    ("a00", 1.0, 5.0),
    ("a01", 1.0, 5.0),
    # basal theta_p: a10*c_pre + a11*c_post
    ("a10", 1.0, 5.0),
    ("a11", 1.0, 5.0),
    # apical theta_d: a20*c_pre + a21*c_post
    ("a20", 1.0, 10.0),
    ("a21", 1.0, 5.0),
    # apical theta_p: a30*c_pre + a31*c_post
    ("a30", 1.0, 10.0),
    ("a31", 1.0, 5.0),
    # CICR phenomenological parameters
    ("tau_mGluR", 100.0, 5000.0),          # mGluR leaky-bucket time constant (ms)
    ("mGluR_thresh", 0.5, 20.0),           # bucket level that activates CICR
    ("a_evt", 0.0, 0.05),                  # cpre coeff for per-synapse event-detection threshold (mM)
    ("b_evt", 0.0, 0.05),                  # cpost coeff for per-synapse event-detection threshold (mM)
    ("a_er", 0.0, 0.05),                   # cpre coeff for ER-loading threshold (mM)
    ("b_er", 0.0, 0.05),                   # cpost coeff for ER-loading threshold (mM)
    ("g_serca_cicr", 0.01, 50.0),          # SERCA pump rate (mM/ms)
    ("g_release_cicr", 0.001, 10.0),       # CICR release rate (1/ms)
]

PARAM_NAMES = [p[0] for p in FIT_PARAMS]
PARAM_BOUNDS = [(p[1], p[2]) for p in FIT_PARAMS]
LOWER = np.array([b[0] for b in PARAM_BOUNDS])
UPPER = np.array([b[1] for b in PARAM_BOUNDS])

DEFAULT_PARAMS = {
    "gamma_d_GB_GluSynapse": 101.5,
    "gamma_p_GB_GluSynapse": 216.2,
    "a00": 1.002, "a01": 1.954,
    "a10": 1.159, "a11": 2.483,
    "a20": 1.127, "a21": 2.456,
    "a30": 5.236, "a31": 1.782,
    # CICR defaults (from cawave.cfg / cawave.py)
    "tau_mGluR": 500.0,
    "mGluR_thresh": 3.0,
    "a_evt": 0.003,   # ca_event_thresh = a_evt*Cpre + b_evt*Cpost; ~0.0003 mM for Cpre~0.09
    "b_evt": 0.003,
    "a_er": 0.005,    # gamma_er = a_er*Cpre + b_er*Cpost; ~0.00045 mM for Cpre~0.09
    "b_er": 0.005,
    "g_serca_cicr": 1.9565,
    "g_release_cicr": 0.5,
}
DEFAULT_X0 = np.array([DEFAULT_PARAMS[name] for name in PARAM_NAMES])

# Paths
BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"
L5_TRACE_DIR = os.path.join(BASE_DIR, "trace_results/Chindemi_params")
L23_TRACE_DIR = os.path.join(BASE_DIR, "trace_results/L23PC_Chindemi_params")
L5_BASIS_DIR = os.path.join(BASE_DIR, "basis_results")
L23_BASIS_DIR = os.path.join(BASE_DIR, "basis_results_L23PC_L5TTPC")

PROTOCOL_PATHWAY = {
    "2Hz_5ms": "L5TTPC",
    "5Hz_5ms": "L5TTPC",
    "10Hz_10ms": "L5TTPC",
    "10Hz_-10ms": "L5TTPC",
    "50Hz_10ms": "L23PC",
}

# ══════════════════════════════════════════════════════════════════
# Module-level globals for multiprocessing-safe objective
# ══════════════════════════════════════════════════════════════════

_PROTOCOL_DATA = None
_TARGETS = None
_LAMBDA_REG = 0.0
_DEFAULT_X0 = None
_WEIGHTS = None


# ══════════════════════════════════════════════════════════════════
# Numba-accelerated rho ODE
# ══════════════════════════════════════════════════════════════════

@njit(cache=True)
def _compute_effcai_from_cai(cai_trace, t, tau_effca, min_ca):
    """
    Compute effcai analytically from cai_CR for a single synapse.

    Uses exponential integrator (exact for piecewise-constant driving):
        effcai[n+1] = effcai[n] * exp(-dt/tau) + (cai[n] - min_ca) * tau * (1 - exp(-dt/tau))
    """
    n_points = len(t)
    effcai = np.zeros(n_points, dtype=np.float64)
    for i in range(n_points - 1):
        dt = t[i + 1] - t[i]
        decay = np.exp(-dt / tau_effca)
        driving = cai_trace[i] - min_ca
        effcai[i + 1] = effcai[i] * decay + driving * tau_effca * (1.0 - decay)
    return effcai


@njit(cache=True)
def _compute_effcai_batch(cai, t, tau_effca, min_ca):
    """Compute effcai analytically from cai_CR for all synapses.

    Parameters
    ----------
    cai : (n_syn, n_time) array of cai_CR traces
    t   : (n_time,) time vector in ms

    Returns
    -------
    effcai : (n_syn, n_time) array of computed effcai_GB
    """
    n_syn = cai.shape[0]
    n_time = cai.shape[1]
    effcai = np.zeros((n_syn, n_time), dtype=np.float64)
    for s in range(n_syn):
        for i in range(n_time - 1):
            dt = t[i + 1] - t[i]
            decay = np.exp(-dt / tau_effca)
            driving = cai[s, i] - min_ca
            effcai[s, i + 1] = effcai[s, i] * decay + driving * tau_effca * (1.0 - decay)
    return effcai


@njit(cache=True)
def _apply_cicr_batch(cai, t, c_pre, c_post,
                      tau_mGluR, mGluR_thresh, a_evt, b_evt,
                      a_er, b_er, g_serca, g_release,
                      max_dt_sub=0.1):
    """Apply phenomenological CICR to calcium traces.

    Model (per synapse):
      1. Event detection: rising-edge crossing of cai > ca_event_thresh,
         where ca_event_thresh = a_evt*cpre + b_evt*cpost  (per synapse)
      2. mGluR leaky bucket: mGluR += event, decays with tau_mGluR
      3. When mGluR > mGluR_thresh -> CICR gate ON
      4. SERCA loads ER when cai_total > gamma_er = a_er*cpre + b_er*cpost
      5. CICR releases ER calcium into cytoplasm when gate ON
      6. cai_total = cai_CR + ca_cicr (accumulated CICR contribution)

    Args:
        max_dt_sub: max substep for ER dynamics (ms). 0.1 for high accuracy,
                    5.0 for fast NN training data generation.
    Returns cai_total (n_syn, n_time).
    """
    K_serca = 0.1       # uM, fixed (cawave.py)
    fc = 0.83
    fe = 0.17
    caAvg = 0.0017      # mM
    MAX_DT_SUB = max_dt_sub

    n_syn = cai.shape[0]
    n_time = cai.shape[1]
    cai_total = np.copy(cai)

    for s in range(n_syn):
        ca_event_thresh_s = a_evt * c_pre[s] + b_evt * c_post[s]
        gamma_er = a_er * c_pre[s] + b_er * c_post[s]
        ca_er = (caAvg - fc * cai[s, 0]) / fe
        if ca_er < 0.0:
            ca_er = 0.0
        mGluR_val = 0.0
        ca_cicr = 0.0
        was_above = 0

        for i in range(n_time - 1):
            dt_full = t[i + 1] - t[i]
            if dt_full <= 0.0:
                cai_total[s, i + 1] = cai[s, i + 1] + ca_cicr
                continue

            # --- event detection on raw cai (rising edge) ---
            is_above = 1 if cai[s, i] > ca_event_thresh_s else 0
            event = 1.0 if (is_above == 1 and was_above == 0) else 0.0
            was_above = is_above

            # --- mGluR leaky bucket (exact exp decay + event) ---
            mGluR_val = mGluR_val * np.exp(-dt_full / tau_mGluR) + event
            cicr_on = 1.0 if mGluR_val > mGluR_thresh else 0.0

            # --- subcycled ER / CICR dynamics ---
            n_sub = max(1, int(np.ceil(dt_full / MAX_DT_SUB)))
            dt_sub = dt_full / n_sub

            for _ in range(n_sub):
                cai_now = cai[s, i] + ca_cicr  # total cai (used for release driving force)

                # SERCA: operates only on ca_cicr, NOT on cai_now.
                # Reason: the NEURON cai_CR trace already contains its own SERCA.
                # Using cai_now would double-count the base pump and create calcium
                # from nothing when ca_cicr is clamped at 0.
                # Conservation cap (ca_cicr / dt_sub) ensures ca_cicr stays >= 0.
                ca_cicr_um = 1000.0 * ca_cicr
                if ca_cicr > gamma_er:
                    J_serca_raw = g_serca * (ca_cicr_um * ca_cicr_um) / (K_serca * K_serca + ca_cicr_um * ca_cicr_um)
                    J_serca = min(J_serca_raw, ca_cicr / dt_sub)
                else:
                    J_serca = 0.0

                # CICR release (gate ON and ER has calcium above total cytoplasmic Ca)
                if cicr_on > 0.5 and ca_er > cai_now:
                    J_release = g_release * (ca_er - cai_now)
                else:
                    J_release = 0.0

                J_net = J_release - J_serca

                # update ER and cytoplasmic CICR contribution (conserved)
                ca_er = ca_er - (fc / fe) * J_net * dt_sub
                if ca_er < 0.0:
                    ca_er = 0.0
                ca_cicr = ca_cicr + J_net * dt_sub
                if ca_cicr < 0.0:
                    ca_cicr = 0.0

            cai_total[s, i + 1] = cai[s, i + 1] + ca_cicr

    return cai_total


@njit(cache=True)
def _compute_rho_final(effcai, t, theta_d, theta_p, gamma_d, gamma_p, rho0):
    """Euler-integrate rho ODE for one synapse, return final value only."""
    n = len(t)
    rho = rho0
    inv_tau = 1.0 / (70.0 * 1000.0)
    rho_star = 0.5

    for i in range(n - 1):
        dt = t[i + 1] - t[i]
        dep = 1.0 if effcai[i] > theta_d else 0.0
        pot = 1.0 if effcai[i] > theta_p else 0.0
        drho = (-rho * (1.0 - rho) * (rho_star - rho)
                + pot * gamma_p * (1.0 - rho)
                - dep * gamma_d * rho) * inv_tau
        rho_new = rho + dt * drho
        if rho_new < 0.0:
            rho_new = 0.0
        elif rho_new > 1.0:
            rho_new = 1.0
        rho = rho_new

    return rho


@njit(cache=True)
def _compute_pair_ratio(effcai, t, c_pre, c_post, is_apical, rho0,
                        baseline_mean, singleton_means,
                        gamma_d, gamma_p,
                        a00, a01,
                        a10, a11,
                        a20, a21,
                        a30, a31):
    """
    Compute EPSP ratio for one pair.

    Threshold equations (linear):
      basal:  theta_d = a00*c_pre + a01*c_post
              theta_p = a10*c_pre + a11*c_post
      apical: theta_d = a20*c_pre + a21*c_post
              theta_p = a30*c_pre + a31*c_post

    effcai: (n_syn, n_time)
    Returns: ratio (float) or NaN if EPSP_before <= 0
    """
    n_syn = effcai.shape[0]

    # EPSP before (from initial rho)
    epsp_before = baseline_mean
    for s in range(n_syn):
        if rho0[s] >= 0.5:
            epsp_before += singleton_means[s] - baseline_mean

    # EPSP after (from rho evolved with candidate params)
    epsp_after = baseline_mean
    for s in range(n_syn):
        cp = c_pre[s]
        cq = c_post[s]
        if is_apical[s]:
            theta_d = a20 * cp + a21 * cq
            theta_p = a30 * cp + a31 * cq
        else:
            theta_d = a00 * cp + a01 * cq
            theta_p = a10 * cp + a11 * cq

        rho_final = _compute_rho_final(
            effcai[s], t, theta_d, theta_p, gamma_d, gamma_p, rho0[s]
        )

        if rho_final >= 0.5:
            epsp_after += singleton_means[s] - baseline_mean

    if epsp_before <= 0.0:
        return np.nan

    return epsp_after / epsp_before


# ══════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════

def _load_pkl(pkl_path):
    """Extract cai_CR, load rho0 & synprops.  effcai computed on-the-fly."""
    with open(pkl_path, "rb") as f:
        data = pickle.load(f)

    t = np.asarray(data["t"], dtype=np.float64)

    # cai_CR -> (n_syn, n_time)
    cai = np.asarray(data["cai_CR"], dtype=np.float64)
    if cai.ndim == 1:
        cai = cai.reshape(1, -1)
    elif cai.shape[0] == len(t) and cai.shape[1] != len(t):
        cai = cai.T

    n_syn = cai.shape[0]
    cai_contig = np.ascontiguousarray(cai)
    t_contig = np.ascontiguousarray(t)

    # rho_GB -> initial rho from trace start
    rho0 = np.zeros(n_syn, dtype=np.float64)
    if "rho_GB" in data:
        rho_gb = np.asarray(data["rho_GB"], dtype=np.float64)
        if rho_gb.ndim == 1:
            rho_gb = rho_gb.reshape(1, -1)
        elif rho_gb.shape[0] == len(t) and rho_gb.shape[1] != len(t):
            rho_gb = rho_gb.T
        rho0 = rho_gb[:, 0].copy()

    # synprops
    synprops = data.get("synprop", {})
    c_pre = np.asarray(synprops.get("Cpre", np.zeros(n_syn)), dtype=np.float64)
    c_post = np.asarray(synprops.get("Cpost", np.zeros(n_syn)), dtype=np.float64)
    loc_list = synprops.get("loc", ["basal"] * n_syn)
    is_apical = np.array([loc == "apical" for loc in loc_list], dtype=np.bool_)

    return {
        "cai": cai_contig,
        "t": t_contig,
        "c_pre": np.ascontiguousarray(c_pre),
        "c_post": np.ascontiguousarray(c_post),
        "is_apical": is_apical,
        "rho0": np.ascontiguousarray(rho0),
    }


def _load_basis(pre_gid, post_gid, basis_dir):
    """Load baseline_mean and singleton_means from basis CSV."""
    csv_path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    if not os.path.exists(csv_path):
        return None

    df = pd.read_csv(csv_path)
    configs = df["config"].apply(lambda x: [int(i) for i in x.split(",")])
    n_syn = len(configs.iloc[0])

    baseline_row = df[configs.apply(lambda x: sum(x) == 0)]
    if baseline_row.empty:
        return None
    baseline_mean = float(baseline_row["mean"].values[0])

    singleton_means = np.zeros(n_syn, dtype=np.float64)
    for i in range(n_syn):
        row = df[configs.apply(lambda x: sum(x) == 1 and x[i] == 1)]
        if row.empty:
            return None
        singleton_means[i] = row["mean"].values[0]

    return {"baseline_mean": baseline_mean, "singleton_means": singleton_means, "n_syn": n_syn}


def preload_all_data():
    """Preload all trace + basis data, organized by protocol."""
    protocol_data = {proto: [] for proto in EXPERIMENTAL_TARGETS}

    # ── L5TTPC pairs (4 protocols) ──
    l5_protocols = [p for p, pw in PROTOCOL_PATHWAY.items() if pw == "L5TTPC"]

    if Path(L5_TRACE_DIR).exists():
        l5_pair_dirs = sorted(
            d for d in Path(L5_TRACE_DIR).iterdir() if d.is_dir()
        )
        logger.info(f"Found {len(l5_pair_dirs)} L5TTPC trace directories")

        for pair_dir in l5_pair_dirs:
            parts = pair_dir.name.split("-")
            if len(parts) != 2:
                continue
            pre_gid, post_gid = int(parts[0]), int(parts[1])

            basis = _load_basis(pre_gid, post_gid, L5_BASIS_DIR)
            if basis is None:
                continue

            for proto in l5_protocols:
                pkl_path = pair_dir / proto / "simulation_traces.pkl"
                if not pkl_path.exists():
                    continue
                try:
                    pd_item = _load_pkl(str(pkl_path))
                    if pd_item["cai"].shape[0] != basis["n_syn"]:
                        logger.warning(
                            f"Synapse mismatch {pre_gid}-{post_gid} {proto}: "
                            f"traces={pd_item['cai'].shape[0]} basis={basis['n_syn']}"
                        )
                        continue
                    pd_item["baseline_mean"] = basis["baseline_mean"]
                    pd_item["singleton_means"] = basis["singleton_means"]
                    protocol_data[proto].append(pd_item)
                except Exception as e:
                    logger.warning(f"Error loading {pkl_path}: {e}")

    # ── L23PC pairs (1 protocol: 50Hz_10ms) ──
    if Path(L23_TRACE_DIR).exists():
        l23_pair_dirs = sorted(
            d for d in Path(L23_TRACE_DIR).iterdir() if d.is_dir()
        )
        logger.info(f"Found {len(l23_pair_dirs)} L23PC trace directories")

        for pair_dir in l23_pair_dirs:
            parts = pair_dir.name.split("-")
            if len(parts) != 2:
                continue
            pre_gid, post_gid = int(parts[0]), int(parts[1])

            basis = _load_basis(pre_gid, post_gid, L23_BASIS_DIR)
            if basis is None:
                continue

            pkl_path = pair_dir / "50Hz_10ms" / "simulation_traces.pkl"
            if not pkl_path.exists():
                continue
            try:
                pd_item = _load_pkl(str(pkl_path))
                if pd_item["cai"].shape[0] != basis["n_syn"]:
                    logger.warning(f"Synapse mismatch L23 {pre_gid}-{post_gid}")
                    continue
                pd_item["baseline_mean"] = basis["baseline_mean"]
                pd_item["singleton_means"] = basis["singleton_means"]
                protocol_data["50Hz_10ms"].append(pd_item)
            except Exception as e:
                logger.warning(f"Error loading L23 {pkl_path}: {e}")

    # Summary
    total = 0
    for proto, pairs in protocol_data.items():
        logger.info(f"  {proto}: {len(pairs)} pairs")
        total += len(pairs)
    logger.info(f"  Total pair-protocol combos: {total}")

    return protocol_data


# ══════════════════════════════════════════════════════════════════
# Objective function (module-level for pickling)
# ══════════════════════════════════════════════════════════════════

def objective(x):
    """
    SSE between predicted and experimental EPSP ratios.

    Uses module-level _PROTOCOL_DATA, _TARGETS, _LAMBDA_REG, _DEFAULT_X0.

    Parameter vector layout (17 params):
      x[0:2]   gamma_d, gamma_p
      x[2:4]   a00, a01       (basal theta_d)
      x[4:6]   a10, a11       (basal theta_p)
      x[6:8]   a20, a21       (apical theta_d)
      x[8:10]  a30, a31       (apical theta_p)
      x[10:17] CICR params
    """
    gamma_d, gamma_p = x[0], x[1]
    a00, a01 = x[2], x[3]
    a10, a11 = x[4], x[5]
    a20, a21 = x[6], x[7]
    a30, a31 = x[8], x[9]

    # CICR parameters
    tau_mGluR = x[10]
    mGluR_thresh = x[11]
    a_evt = x[12]
    b_evt = x[13]
    a_er = x[14]
    b_er = x[15]
    g_serca_cicr = x[16]
    g_release_cicr = x[17]

    total_loss = 0.0

    for p_idx, (proto, target_ratio) in enumerate(_TARGETS.items()):
        ratios = []
        for pd_item in _PROTOCOL_DATA[proto]:
            # Apply CICR to raw cai -> cai_total
            cai_total = _apply_cicr_batch(
                pd_item["cai"], pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                tau_mGluR, mGluR_thresh, a_evt, b_evt,
                a_er, b_er, g_serca_cicr, g_release_cicr,
            )
            # Compute effcai from CICR-modified calcium
            effcai = _compute_effcai_batch(
                cai_total, pd_item["t"], TAU_EFFCA_GB, MIN_CA_CR,
            )
            ratio = _compute_pair_ratio(
                effcai, pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                pd_item["is_apical"], pd_item["rho0"],
                pd_item["baseline_mean"], pd_item["singleton_means"],
                gamma_d, gamma_p,
                a00, a01,
                a10, a11,
                a20, a21,
                a30, a31,
            )
            if not np.isnan(ratio):
                ratios.append(ratio)

        if ratios:
            mean_ratio = np.mean(ratios)
            w = _WEIGHTS[p_idx] if _WEIGHTS is not None else 1.0
            total_loss += w * (mean_ratio - target_ratio) ** 2

    if _LAMBDA_REG > 0 and _DEFAULT_X0 is not None:
        total_loss += _LAMBDA_REG * np.sum((x - _DEFAULT_X0) ** 2)

    return total_loss


def _init_objective(protocol_data, targets, lambda_reg=0.0, default_x0=None, weights=None):
    """Set module-level globals for the objective function."""
    global _PROTOCOL_DATA, _TARGETS, _LAMBDA_REG, _DEFAULT_X0, _WEIGHTS
    _PROTOCOL_DATA = protocol_data
    _TARGETS = targets
    _LAMBDA_REG = lambda_reg
    _DEFAULT_X0 = default_x0
    _WEIGHTS = weights


# ══════════════════════════════════════════════════════════════════
# Optimizers
# ══════════════════════════════════════════════════════════════════

def run_de(max_iter=1000, workers=1, seed=42, popsize=15, **kw):
    """Differential Evolution (scipy)."""
    from scipy.optimize import differential_evolution

    logger.info(f"Running Differential Evolution "
                f"(maxiter={max_iter}, popsize={popsize}, workers={workers})")

    callback_history = []
    _start = time.time()

    def callback(xk, convergence):
        loss = objective(xk)
        elapsed = time.time() - _start
        callback_history.append({"iter": len(callback_history), "loss": float(loss)})
        logger.info(f"  DE gen {len(callback_history):>4d}: loss={loss:.10f}  "
                     f"conv={convergence:.6f}  [{elapsed/60:.1f}min]")

    result = differential_evolution(
        objective, PARAM_BOUNDS,
        maxiter=max_iter,
        popsize=popsize,
        tol=1e-10,
        seed=seed,
        workers=workers if workers > 1 else 1,
        disp=True,
        polish=True,
        callback=callback,
        updating="deferred" if workers > 1 else "immediate",
    )

    return {
        "method": "differential_evolution",
        "x": result.x.tolist(),
        "fun": float(result.fun),
        "nfev": result.nfev,
        "nit": result.nit,
        "success": result.success,
        "message": result.message,
        "history": callback_history,
    }


def run_cmaes(max_iter=1000, seed=42, **kw):
    """CMA-ES optimization."""
    import cma

    logger.info(f"Running CMA-ES (maxiter={max_iter})")

    x0 = DEFAULT_X0.copy()
    lower = np.array([b[0] for b in PARAM_BOUNDS])
    upper = np.array([b[1] for b in PARAM_BOUNDS])
    sigma0 = 0.3 * np.mean(upper - lower)

    opts = {
        "bounds": [lower.tolist(), upper.tolist()],
        "maxiter": max_iter,
        "tolfun": 1e-12,
        "seed": seed,
        "verb_disp": 100,
    }

    es = cma.CMAEvolutionStrategy(x0.tolist(), sigma0, opts)
    es.optimize(objective)

    res = es.result
    return {
        "method": "cma-es",
        "x": list(res.xbest),
        "fun": float(res.fbest),
        "nfev": res.evaluations,
        "iterations": res.iterations,
    }


def run_optuna(n_trials=5000, seed=42, **kw):
    """Optuna (TPE) Bayesian optimization."""
    import optuna

    logger.info(f"Running Optuna TPE (n_trials={n_trials})")
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def optuna_obj(trial):
        x = np.array([
            trial.suggest_float(name, lo, hi)
            for name, lo, hi in FIT_PARAMS
        ])
        return objective(x)

    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="minimize", sampler=sampler)

    # Seed with default params
    study.enqueue_trial({name: DEFAULT_PARAMS[name] for name in PARAM_NAMES})

    study.optimize(optuna_obj, n_trials=n_trials, show_progress_bar=True)

    best = study.best_trial
    x_best = [best.params[name] for name in PARAM_NAMES]

    return {
        "method": "optuna",
        "x": x_best,
        "fun": float(best.value),
        "nfev": len(study.trials),
    }


def run_lbfgsb(n_starts=50, seed=42, **kw):
    """Multi-start L-BFGS-B."""
    from scipy.optimize import minimize

    logger.info(f"Running multi-start L-BFGS-B (n_starts={n_starts})")

    rng = np.random.default_rng(seed)
    lower = np.array([b[0] for b in PARAM_BOUNDS])
    upper = np.array([b[1] for b in PARAM_BOUNDS])

    best_result = None

    for i in range(n_starts):
        x0 = DEFAULT_X0.copy() if i == 0 else rng.uniform(lower, upper)

        result = minimize(
            objective, x0, method="L-BFGS-B", bounds=PARAM_BOUNDS,
            options={"maxiter": 500, "ftol": 1e-14},
        )

        if best_result is None or result.fun < best_result.fun:
            best_result = result
            logger.info(f"  Start {i+1}/{n_starts}: loss={result.fun:.10f} * (new best)")

    return {
        "method": "L-BFGS-B",
        "x": best_result.x.tolist(),
        "fun": float(best_result.fun),
        "nfev": best_result.nfev,
        "success": best_result.success,
    }


# ══════════════════════════════════════════════════════════════════
# Neural network surrogate fitting
# ══════════════════════════════════════════════════════════════════

ALL_PROTO_NAMES = list(EXPERIMENTAL_TARGETS.keys())
ALL_TARGETS = np.array(list(EXPERIMENTAL_TARGETS.values()))


def _update_proto_globals(targets):
    """Update ALL_PROTO_NAMES/ALL_TARGETS to match filtered targets."""
    global ALL_PROTO_NAMES, ALL_TARGETS
    ALL_PROTO_NAMES = list(targets.keys())
    ALL_TARGETS = np.array(list(targets.values()))


def _compute_ratios_cicr(x):
    """Forward model: params (17,) -> mean EPSP ratios (5,) via CICR pipeline.

    Runs the full CICR -> effcai -> rho -> ratio pipeline for all protocols.
    Uses module-level _PROTOCOL_DATA.
    Uses relaxed substep (MAX_DT_SUB=5.0) for speed — acceptable for NN training.
    """
    gamma_d, gamma_p = x[0], x[1]
    a00, a01 = x[2], x[3]
    a10, a11 = x[4], x[5]
    a20, a21 = x[6], x[7]
    a30, a31 = x[8], x[9]
    tau_mGluR = x[10]
    mGluR_thresh = x[11]
    a_evt = x[12]
    b_evt = x[13]
    a_er = x[14]
    b_er = x[15]
    g_serca_cicr = x[16]
    g_release_cicr = x[17]

    ratios = np.full(len(ALL_PROTO_NAMES), np.nan)
    for p_idx, proto in enumerate(ALL_PROTO_NAMES):
        if proto not in _PROTOCOL_DATA:
            continue
        pair_ratios = []
        for pd_item in _PROTOCOL_DATA[proto]:
            cai_total = _apply_cicr_batch(
                pd_item["cai"], pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                tau_mGluR, mGluR_thresh, a_evt, b_evt,
                a_er, b_er, g_serca_cicr, g_release_cicr,
                5.0,  # relaxed substep for NN training speed
            )
            effcai = _compute_effcai_batch(
                cai_total, pd_item["t"], TAU_EFFCA_GB, MIN_CA_CR,
            )
            r = _compute_pair_ratio(
                effcai, pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                pd_item["is_apical"], pd_item["rho0"],
                pd_item["baseline_mean"], pd_item["singleton_means"],
                gamma_d, gamma_p,
                a00, a01,
                a10, a11,
                a20, a21,
                a30, a31,
            )
            if not np.isnan(r):
                pair_ratios.append(r)
        if pair_ratios:
            ratios[p_idx] = np.mean(pair_ratios)
    return ratios


def generate_training_data_cicr(n_samples, workers=1, seed=42):
    """Generate (params, ratios) pairs via Sobol quasi-random sampling."""
    from scipy.stats.qmc import Sobol

    logger.info(f"Sampling {n_samples} parameter vectors (Sobol, {len(PARAM_NAMES)}D)...")

    m = int(np.ceil(np.log2(n_samples)))
    n_sobol = 2 ** m

    sampler = Sobol(d=len(PARAM_NAMES), scramble=True, seed=seed)
    X_unit = sampler.random(n_sobol)[:n_samples]
    X = LOWER + X_unit * (UPPER - LOWER)
    X[0] = DEFAULT_X0  # always include default

    logger.info(f"Running CICR forward model for {n_samples} samples with {workers} workers...")
    t0 = time.time()

    if workers > 1:
        with Pool(workers) as pool:
            Y_list = list(pool.imap(_compute_ratios_cicr, X, chunksize=64))
    else:
        Y_list = []
        for i, xi in enumerate(X):
            Y_list.append(_compute_ratios_cicr(xi))
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


def train_surrogate_cicr(X, Y, seed=42):
    """Train sklearn MLPRegressor surrogate: params (25) -> ratios (5)."""
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.model_selection import train_test_split

    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    X_train, X_val, Y_train, Y_val = train_test_split(
        X_scaled, Y, test_size=0.15, random_state=seed,
    )

    logger.info(f"Training MLP surrogate: {len(X_train)} train, {len(X_val)} val, "
                f"{X.shape[1]} inputs -> {Y.shape[1]} outputs")

    model = MLPRegressor(
        hidden_layer_sizes=(256, 256, 128),
        activation="relu",
        solver="adam",
        alpha=1e-4,
        batch_size=256,
        learning_rate="adaptive",
        learning_rate_init=1e-3,
        max_iter=1000,
        tol=1e-5,
        early_stopping=True,
        validation_fraction=0.15,
        n_iter_no_change=20,
        random_state=seed,
        verbose=True,
    )

    model.fit(X_train, Y_train)

    Y_pred_val = model.predict(X_val)
    mae_per_proto = np.mean(np.abs(Y_pred_val - Y_val), axis=0)
    mse_val = np.mean((Y_pred_val - Y_val) ** 2)

    logger.info(f"  Converged in {model.n_iter_} iterations")
    logger.info(f"  Validation MSE: {mse_val:.8f}")
    logger.info(f"  Validation MAE per protocol:")
    for i, proto in enumerate(ALL_PROTO_NAMES):
        logger.info(f"    {proto}: {mae_per_proto[i]:.6f}")

    return model, scaler


def optimize_surrogate_cicr(model, scaler, n_starts=500, top_k=20, seed=42, workers=1):
    """Multi-start L-BFGS-B over NN surrogate. Returns top_k diverse candidates."""
    import heapq

    logger.info(f"Optimizing over surrogate ({n_starts} L-BFGS-B starts, "
                f"keeping top {top_k}, workers={workers})...")
    t0 = time.time()

    rng = np.random.default_rng(seed)
    x0_list = [DEFAULT_X0.copy()] + [rng.uniform(LOWER, UPPER) for _ in range(n_starts - 1)]
    targets = ALL_TARGETS
    weights = np.ones(len(ALL_PROTO_NAMES))

    if workers > 1:
        from joblib import Parallel, delayed
        from scipy.optimize import minimize as _minimize

        def _run_one(x0):
            def nn_obj(x):
                x_scaled = scaler.transform(x.reshape(1, -1))
                pred = np.atleast_1d(model.predict(x_scaled)[0])
                return float(np.sum(weights * (pred - targets) ** 2))
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
        from scipy.optimize import minimize
        top_candidates = []
        best_loss = np.inf

        def nn_objective(x):
            x_scaled = scaler.transform(x.reshape(1, -1))
            pred = np.atleast_1d(model.predict(x_scaled)[0])
            return float(np.sum(weights * (pred - targets) ** 2))

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

    candidates = sorted([(-loss, x) for loss, x in top_candidates])
    best_params = candidates[0][1]
    best_surr_loss = candidates[0][0]

    x_scaled = scaler.transform(best_params.reshape(1, -1))
    nn_preds = np.atleast_1d(model.predict(x_scaled)[0])

    logger.info(f"\n  Best surrogate SSE: {best_surr_loss:.12f}")
    logger.info(f"  Top {top_k} surrogate SSE range: [{candidates[0][0]:.8f}, {candidates[-1][0]:.8f}]")
    for i, proto in enumerate(ALL_PROTO_NAMES):
        logger.info(f"    {proto}: nn={nn_preds[i]:.6f} target={targets[i]:.6f} "
                     f"err={nn_preds[i]-targets[i]:+.6f}")

    return best_params, best_surr_loss, nn_preds, [x for _, x in candidates]


def refine_on_real_model_cicr(candidates, default_loss, workers=1,
                              refine=False, max_time_per_candidate=120):
    """Evaluate top-K NN candidates on real CICR model, optionally refine with Nelder-Mead."""
    from multiprocessing import Pool

    logger.info(f"\nEvaluating {len(candidates)} NN candidates on real CICR model "
                f"(workers={workers})...")

    all_xs = [DEFAULT_X0.copy()] + [c.copy() for c in candidates]
    labels = ["default"] + [f"NN-{i+1}" for i in range(len(candidates))]

    # Parallel evaluation on real CICR model
    if workers > 1:
        with Pool(workers) as pool:
            losses = pool.map(objective, all_xs)
    else:
        losses = [objective(x) for x in all_xs]

    scored = []
    for loss, label, x0 in zip(losses, labels, all_xs):
        scored.append((loss, label, x0))
        logger.info(f"  {label}: real SSE = {loss:.8f}")

    scored.sort(key=lambda t: t[0])
    logger.info(f"\n  Best 3 by direct eval:")
    for loss, label, _ in scored[:3]:
        logger.info(f"    {label}: {loss:.8f}")

    best_params = scored[0][2].copy()
    best_loss = scored[0][0]

    if not refine:
        logger.info(f"\n  Best real SSE: {best_loss:.10f} "
                    f"(default was {default_loss:.10f})")
        return best_params, best_loss

    # Optional Nelder-Mead refinement
    from scipy.optimize import minimize

    n_refine = min(3, len(scored))
    logger.info(f"\nRefining top {n_refine} candidates with Nelder-Mead "
                f"(max {max_time_per_candidate}s each)...")

    for rank in range(n_refine):
        init_loss, label, x0 = scored[rank]
        t0 = time.time()
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
            final_loss = _best_so_far[0]
            final_x = x0
            nfev = _eval_count[0]

        elapsed = time.time() - t0
        logger.info(f"  {label}: {init_loss:.8f} -> {final_loss:.8f} "
                     f"({nfev} evals, {elapsed:.0f}s)")

        if final_loss < best_loss:
            best_loss = final_loss
            best_params = final_x
            logger.info(f"    ^ new best real SSE: {best_loss:.10f}")

    logger.info(f"\n  Best real SSE after refinement: {best_loss:.10f} "
                f"(default was {default_loss:.10f})")
    return best_params, best_loss


def run_nn(n_samples=50000, workers=1, n_starts=500, seed=42,
           data_path=None, model_path=None, save_dir="fitting_results", **kw):
    """NN surrogate pipeline: generate data -> train -> optimize -> refine."""
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    default_loss = objective(DEFAULT_X0)
    logger.info(f"Default params real loss: {default_loss:.10f}")

    # Warmup the forward model
    _ = _compute_ratios_cicr(DEFAULT_X0)

    # ── Generate or load training data ──
    if data_path:
        logger.info(f"Loading training data from {data_path}")
        npz = np.load(data_path)
        X, Y = npz["X"], npz["Y"]
        valid = ~np.any(np.isnan(Y), axis=1)
        X, Y = X[valid], Y[valid]
        logger.info(f"Loaded {len(X)} valid samples")
    else:
        X, Y = generate_training_data_cicr(n_samples, workers, seed)
        dp = save_dir / "training_data_cicr.npz"
        np.savez(dp, X=X, Y=Y)
        logger.info(f"Saved training data to {dp}")

    # ── Train or load surrogate ──
    if model_path:
        logger.info(f"Loading model from {model_path}")
        with open(model_path, "rb") as f:
            saved = pickle.load(f)
        model, scaler = saved["model"], saved["scaler"]
    else:
        model, scaler = train_surrogate_cicr(X, Y, seed=seed)
        mp = save_dir / "surrogate_cicr.pkl"
        with open(mp, "wb") as f:
            pickle.dump({"model": model, "scaler": scaler,
                         "protocols": ALL_PROTO_NAMES}, f)
        logger.info(f"Saved model to {mp}")

    # ── Optimize over surrogate ──
    t0 = time.time()
    best_params, surr_loss, nn_preds, top_candidates = optimize_surrogate_cicr(
        model, scaler, n_starts=n_starts, seed=seed, workers=workers,
    )
    opt_time = time.time() - t0

    # ── Refine on real CICR model ──
    best_params, real_loss = refine_on_real_model_cicr(top_candidates, default_loss, workers=workers)
    real_eval = evaluate_params(best_params)

    print(f"\n{'='*70}")
    print("  Neural Network Surrogate Results (CICR)")
    print(f"{'='*70}")
    print(f"  Training samples:  {len(X)}")
    print(f"  Surrogate SSE:     {surr_loss:.10f}")
    print(f"  Real SSE:          {real_loss:.10f}")
    print(f"  Default SSE:       {default_loss:.10f}")
    print(f"  Surrogate opt:     {opt_time:.1f}s")

    print(f"\n  Best parameters:")
    for i, name in enumerate(PARAM_NAMES):
        d = DEFAULT_PARAMS[name]
        b = best_params[i]
        pct = (b - d) / d * 100 if d != 0 else 0
        print(f"    {name:30s} = {b:10.4f}  (default: {d:.4f}, {pct:+.1f}%)")

    print(f"\n  {'Protocol':<15s} {'NN':>8s} {'Real':>8s} {'Exper':>8s} {'RealErr':>9s}")
    print(f"  {'-'*48}")
    for i, proto in enumerate(ALL_PROTO_NAMES):
        r = real_eval[proto]
        print(f"  {proto:<15s} {nn_preds[i]:8.4f} {r['predicted']:8.4f} "
              f"{r['experimental']:8.4f} {r['error']:+9.4f}")
    print(f"{'='*70}\n")

    return {
        "method": "neural_network_surrogate",
        "x": best_params.tolist(),
        "fun": float(real_loss),
        "fun_surrogate": float(surr_loss),
        "fun_default": float(default_loss),
        "nn_predictions": nn_preds.tolist(),
        "per_protocol": real_eval,
        "nfev": len(X),
        "n_training_samples": len(X),
        "opt_time_seconds": opt_time,
    }


OPTIMIZERS = {
    "de": run_de,
    "cmaes": run_cmaes,
    "optuna": run_optuna,
    "lbfgsb": run_lbfgsb,
    "nn": run_nn,
}


# ══════════════════════════════════════════════════════════════════
# Analysis & output
# ══════════════════════════════════════════════════════════════════

def evaluate_params(x):
    """Return per-protocol predictions for a given parameter vector."""
    gamma_d, gamma_p = x[0], x[1]
    a00, a01 = x[2], x[3]
    a10, a11 = x[4], x[5]
    a20, a21 = x[6], x[7]
    a30, a31 = x[8], x[9]

    # CICR parameters
    tau_mGluR = x[10]
    mGluR_thresh = x[11]
    a_evt = x[12]
    b_evt = x[13]
    a_er = x[14]
    b_er = x[15]
    g_serca_cicr = x[16]
    g_release_cicr = x[17]

    results = {}
    for proto, target in _TARGETS.items():
        ratios = []
        for pd_item in _PROTOCOL_DATA[proto]:
            cai_total = _apply_cicr_batch(
                pd_item["cai"], pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                tau_mGluR, mGluR_thresh, a_evt, b_evt,
                a_er, b_er, g_serca_cicr, g_release_cicr,
            )
            effcai = _compute_effcai_batch(
                cai_total, pd_item["t"], TAU_EFFCA_GB, MIN_CA_CR,
            )
            ratio = _compute_pair_ratio(
                effcai, pd_item["t"],
                pd_item["c_pre"], pd_item["c_post"],
                pd_item["is_apical"], pd_item["rho0"],
                pd_item["baseline_mean"], pd_item["singleton_means"],
                gamma_d, gamma_p,
                a00, a01,
                a10, a11,
                a20, a21,
                a30, a31,
            )
            if not np.isnan(ratio):
                ratios.append(ratio)

        mean_r = float(np.mean(ratios)) if ratios else float("nan")
        results[proto] = {
            "predicted": mean_r,
            "experimental": float(target),
            "error": float(mean_r - target),
            "n_pairs": len(ratios),
        }
    return results


def print_results(opt_result, eval_results):
    """Pretty-print optimization results."""
    print(f"\n{'='*70}")
    print(f"  Method:    {opt_result['method']}")
    print(f"  Loss(SSE): {opt_result['fun']:.10f}")
    print(f"  Evals:     {opt_result.get('nfev', 'N/A')}")
    elapsed = opt_result.get("elapsed_seconds", 0)
    print(f"  Time:      {elapsed/60:.1f} min")
    print(f"{'='*70}")

    print("\n  Best parameters:")
    for i, name in enumerate(PARAM_NAMES):
        default = DEFAULT_PARAMS[name]
        best = opt_result["x"][i]
        pct = (best - default) / default * 100 if default != 0 else 0
        print(f"    {name:30s} = {best:10.4f}  (default: {default:.4f}, {pct:+.1f}%)")

    print(f"\n  {'Protocol':<15s} {'Predicted':>10s} {'Experiment':>11s} {'Error':>10s} {'N':>5s}")
    print(f"  {'-'*51}")
    for proto in eval_results:
        r = eval_results[proto]
        print(f"  {proto:<15s} {r['predicted']:10.4f} {r['experimental']:11.4f} "
              f"{r['error']:+10.4f} {r['n_pairs']:5d}")
    print(f"{'='*70}\n")


# ══════════════════════════════════════════════════════════════════
# SLURM script generation
# ══════════════════════════════════════════════════════════════════

def generate_slurm_script(args):
    """Generate SLURM submission script to run all 4 methods in parallel."""
    script_path = os.path.join(BASE_DIR, "new_fitting", "submit_fitting.sh")
    output_dir = os.path.join(BASE_DIR, "new_fitting", "fitting_results")

    script = f"""#!/bin/bash
#SBATCH --job-name=fit_params
#SBATCH --account=ctb-emuller
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output={output_dir}/fit_%j.out
#SBATCH --error={output_dir}/fit_%j.err

source {BASE_DIR}/setupenv.sh
mkdir -p {output_dir}

SCRIPT="{BASE_DIR}/new_fitting/fit_params.py"

echo "Starting parameter fitting at $(date)"
echo "============================================"

# Run all 4 optimizers sequentially (each uses JIT-compiled ODE)
echo "\\n--- Differential Evolution ---"
python $SCRIPT --method de --max-iter 1000 --workers 8 \\
    --output {output_dir}/result_de.json

echo "\\n--- CMA-ES ---"
python $SCRIPT --method cmaes --max-iter 1000 \\
    --output {output_dir}/result_cmaes.json

echo "\\n--- Optuna ---"
python $SCRIPT --method optuna --n-trials 5000 \\
    --output {output_dir}/result_optuna.json

echo "\\n--- L-BFGS-B ---"
python $SCRIPT --method lbfgsb --n-starts 100 \\
    --output {output_dir}/result_lbfgsb.json

echo "\\nAll methods completed at $(date)"
"""
    with open(script_path, "w") as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    print(f"SLURM script written to: {script_path}")
    print(f"Submit with: sbatch {script_path}")
    return script_path


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Fit Graupner-Brunel plasticity parameters to experimental EPSP ratios",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--method", choices=list(OPTIMIZERS.keys()) + ["all"],
                        help="Optimization method (or 'all' to run all 5)")
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--n-trials", type=int, default=5000, help="Optuna trials")
    parser.add_argument("--n-starts", type=int, default=50, help="L-BFGS-B / NN surrogate restarts")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers (DE / NN)")
    parser.add_argument("--popsize", type=int, default=15, help="DE population multiplier")
    parser.add_argument("--lambda-reg", type=float, default=0.0,
                        help="L2 regularization toward default params")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default=None, help="Output JSON path")
    parser.add_argument("--generate-slurm", action="store_true",
                        help="Generate SLURM script and exit")
    parser.add_argument("--eval-only", action="store_true",
                        help="Only evaluate default params (no optimization)")
    # NN surrogate args
    parser.add_argument("--n-samples", type=int, default=50000,
                        help="Number of Sobol samples for NN training data")
    parser.add_argument("--data", type=str, default=None,
                        help="Load pre-generated training data .npz (skip generation)")
    parser.add_argument("--nn-model", type=str, default=None,
                        help="Load pre-trained NN surrogate .pkl (skip training)")
    parser.add_argument("--protocols", nargs="+", default=None,
                        choices=list(EXPERIMENTAL_TARGETS.keys()),
                        help="Subset of protocols to fit (default: all 5)")

    args = parser.parse_args()

    if args.generate_slurm:
        generate_slurm_script(args)
        return

    if not args.method and not args.eval_only:
        parser.print_help()
        return

    # 1. Preload data
    logger.info("Preloading trace and basis data...")
    t0 = time.time()
    protocol_data = preload_all_data()
    load_time = time.time() - t0
    logger.info(f"Data loaded in {load_time:.1f}s")

    # 2. Warmup numba JIT
    logger.info("Warming up numba JIT...")
    _dummy_cai = np.zeros((1, 100), dtype=np.float64)
    _dummy_t = np.linspace(0, 100, 100, dtype=np.float64)
    _compute_effcai_from_cai(_dummy_cai[0], _dummy_t, TAU_EFFCA_GB, MIN_CA_CR)
    _dummy_e = _compute_effcai_batch(_dummy_cai, _dummy_t, TAU_EFFCA_GB, MIN_CA_CR)
    _compute_rho_final(_dummy_e[0], _dummy_t, 0.1, 0.2, 100.0, 200.0, 0.0)
    _compute_pair_ratio(
        _dummy_e, _dummy_t,
        np.zeros(1), np.zeros(1), np.array([False]), np.zeros(1),
        3.0, np.array([4.0]),
        100.0, 200.0,
        1.0, 2.0,   # basal theta_d: a00, a01
        1.0, 2.0,   # basal theta_p: a10, a11
        1.0, 2.0,   # apical theta_d: a20, a21
        5.0, 2.0,   # apical theta_p: a30, a31
    )
    _apply_cicr_batch(
        _dummy_cai, _dummy_t, np.zeros(1), np.zeros(1),
        500.0, 3.0, 0.003, 0.003, 0.005, 0.005, 1.9565, 0.5,
    )
    logger.info("JIT ready")

    # 3. Initialize objective globals (with optional protocol filter)
    targets = EXPERIMENTAL_TARGETS
    if args.protocols:
        targets = {p: v for p, v in EXPERIMENTAL_TARGETS.items() if p in args.protocols}
        logger.info(f"Fitting {len(targets)} protocols: {list(targets.keys())}")
    _update_proto_globals(targets)
    _init_objective(protocol_data, targets, args.lambda_reg, DEFAULT_X0)

    # Sanity check: default params
    default_loss = objective(DEFAULT_X0)
    logger.info(f"Default params loss: {default_loss:.10f}")

    default_eval = evaluate_params(DEFAULT_X0)
    print("\n  Default parameters evaluation:")
    print(f"  {'Protocol':<15s} {'Predicted':>10s} {'Experiment':>11s} {'Error':>10s} {'N':>5s}")
    print(f"  {'-'*51}")
    for proto in targets:
        r = default_eval[proto]
        print(f"  {proto:<15s} {r['predicted']:10.4f} {r['experimental']:11.4f} "
              f"{r['error']:+10.4f} {r['n_pairs']:5d}")
    print()

    if args.eval_only:
        return

    # 4. Run optimizer(s)
    methods = list(OPTIMIZERS.keys()) if args.method == "all" else [args.method]
    all_results = {}

    for method in methods:
        logger.info(f"\n{'='*50}")
        logger.info(f"Starting: {method}")
        logger.info(f"{'='*50}")

        t0 = time.time()
        kwargs = dict(
            max_iter=args.max_iter,
            n_trials=args.n_trials,
            n_starts=args.n_starts,
            workers=args.workers,
            popsize=args.popsize,
            seed=args.seed,
        )
        if method == "nn":
            kwargs.update(
                n_samples=args.n_samples,
                data_path=args.data,
                model_path=args.nn_model,
                save_dir=str(Path(args.output).parent) if args.output else "fitting_results",
            )
        opt_result = OPTIMIZERS[method](**kwargs)
        opt_result["elapsed_seconds"] = time.time() - t0

        best_x = np.array(opt_result["x"])
        eval_result = evaluate_params(best_x)
        print_results(opt_result, eval_result)

        all_results[method] = {
            "optimization": opt_result,
            "per_protocol": eval_result,
        }

    # 5. Compare (if multiple methods)
    if len(methods) > 1:
        print(f"\n{'='*70}")
        print("  COMPARISON SUMMARY")
        print(f"{'='*70}")
        print(f"  {'Method':<25s} {'Loss':>12s} {'Time (min)':>12s} {'Evals':>8s}")
        print(f"  {'-'*57}")
        for method in methods:
            r = all_results[method]["optimization"]
            print(f"  {r['method']:<25s} {r['fun']:12.10f} "
                  f"{r.get('elapsed_seconds',0)/60:12.1f} {r.get('nfev','?'):>8}")
        print(f"{'='*70}\n")

    # 6. Save
    output_path = args.output
    if output_path is None and len(methods) == 1:
        pass  # no auto-save
    elif output_path is None and len(methods) > 1:
        output_path = os.path.join(BASE_DIR, "new_fitting", "fitting_results", "all_results.json")

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        save_data = {
            "param_names": PARAM_NAMES,
            "default_params": DEFAULT_PARAMS,
            "experimental_targets": EXPERIMENTAL_TARGETS,
            "results": {},
        }
        for method, res in all_results.items():
            # Remove non-serializable history if present
            opt = dict(res["optimization"])
            opt.pop("history", None)
            save_data["results"][method] = {
                "optimization": opt,
                "per_protocol": res["per_protocol"],
            }
        with open(output_path, "w") as f:
            json.dump(save_data, f, indent=2)
        logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
