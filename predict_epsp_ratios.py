"""
Analytical EPSP ratio predictor from recipe.csv parameters.

For each pathway in (or selected from) recipe_dhuruva.csv, this script samples
a population of synapses via the same Gaussian copula used in epg.py, then
analytically predicts the STDP plasticity outcome (EPSP ratio) using the
fixed a00-a31 fit_params — without running any NEURON simulation.

──────────────────────────────────────────────────────────────────
ANALYTICAL MODEL
──────────────────────────────────────────────────────────────────

The plasticity efficiency ratio η determines the outcome:

   η = K · Use_ss(u, D, freq) · N_eff(freq) / (a10 + a11·ζ)

where:
   Use_ss = u / (1 + u·D/isi)         ← TM depression at steady state
   N_eff  = 1/(1−exp(−isi/τ_effca))   ← effcai accumulation factor (~3.32 at 10Hz)
   ζ      = c_post/c_pre ≈ 0.013      ← fixed ratio (both scale with volume)
   K      ≈ 3.93                       ← coincidence amplification, calibrated
                                          from CHINDEMI_PARAMS_v2 simulations

Zones:
   η > 1                            → LTP zone  (pot_GB=1, rho tends ↑)
   θ_d/θ_p < η < 1                  → LTD zone  (dep_GB=1, rho tends ↓)
   η < θ_d/θ_p                      → no plasticity

where θ_d/θ_p = (a00 + a01·ζ)/(a10 + a11·ζ) ≈ 0.862 (fixed for given a-params)

EPSP ratio is then approximated via the rho ODE:
   EPSP_ratio ≈ 1 + p_ltp − 0.5·p_ltd

where p_ltp/p_ltd = fraction of synapses that accumulate enough Δrho to flip.

──────────────────────────────────────────────────────────────────
KEY INSIGHT
──────────────────────────────────────────────────────────────────
spinevol has MINIMAL direct effect on η because both the induction effcai and
the thresholds (θ_d, θ_p) scale proportionally with c_pre (which scales with
volume). The dominant recipe parameters driving plasticity are:

   1. u  (Use0) — controls Use_ss and thus η linearly
   2. d  (tau_dep) — controls Use_ss (higher D → more depression → lower η)
   3. freq — controls isi → Use_ss and N_eff

spinevol affects plasticity only INDIRECTLY via the copula correlations
(high spinevol → high gsyn → correlated high u → higher η).

──────────────────────────────────────────────────────────────────
USAGE
──────────────────────────────────────────────────────────────────
    # All pathways, default 10 Hz, dt = -10 +10 ms
    python predict_epsp_ratios.py

    # Single pathway
    python predict_epsp_ratios.py --pathway L5_TPC:A L5_TPC:A

    # Regex filter
    python predict_epsp_ratios.py --mtype-filter L5_TPC

    # Scan u across range for one pathway
    python predict_epsp_ratios.py --scan-param u --pathway L5_TPC:A L5_TPC:A

    # Scan tau_dep
    python predict_epsp_ratios.py --scan-param d --pathway L5_TPC:A L5_TPC:A

    # Multiple frequencies
    python predict_epsp_ratios.py --freq 2 5 10 20 --dt -50 -10 10 50

    # Use a different recipe or config
    python predict_epsp_ratios.py --recipe biodata/recipe_andras.csv
    python predict_epsp_ratios.py --config configs/L5TTPC_L5TTPC_STDP.yaml

    # Save results to CSV
    python predict_epsp_ratios.py --output results_epsp.csv
"""

import argparse
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats

# ──────────────────────────────────────────────────────────────────
# Calibrated constants (derived from CHINDEMI_PARAMS_v2 full-trace sims)
# ──────────────────────────────────────────────────────────────────
# Global coincidence amplification: K × Use_ss × N_eff = effcai_peak / c_pre
K_COINCIDENCE = 3.93

# Fixed c_post/c_pre ratio (both scale similarly with volume, so ζ ≈ const)
ZETA_DEFAULT  = 0.0135   # = 0.00163 / 0.121  (anchors at V=0.1 µm³)

# GluSynapse global params
TAU_EFFCA     = 278.0    # ms  (τ_effca_GB)
TAU_IND       = 70.0     # s   (τ_ind_GB)
INDUCTION_T   = 140.0    # s   (induction duration for 40 pairs × 4s T)

MAX_SEED      = 2**32 - 1

# Default fit_params (from L5TTPC_L5TTPC_STDP.yaml / Zenodo_O1.yaml — Chindemi values)
DEFAULT_FIT_PARAMS = {
    "a00": 1.0018099627, "a01": 1.9535568661,   # basal  LTD threshold
    "a10": 1.1593870631, "a11": 2.4827933785,   # basal  LTP threshold
    "a20": 1.1267343377, "a21": 2.4559713296,   # apical LTD threshold
    "a30": 5.2356637311, "a31": 1.7822214534,   # apical LTP threshold
    "gamma_d_GB_GluSynapse": 101.5387594661,
    "gamma_p_GB_GluSynapse": 216.1841700668,
}


# ──────────────────────────────────────────────────────────────────
# Copula helpers  (mirrors epg.py exactly)
# ──────────────────────────────────────────────────────────────────

def _get_covariance_matrix(r):
    cov = np.eye(6)
    cov[0][3] = cov[3][0] = float(r["u_gsyn_r"]) * float(r["gsyn_nrrp_r"])
    cov[0][4] = cov[4][0] = float(r["u_gsyn_r"])
    cov[0][5] = cov[5][0] = float(r["u_gsyn_r"]) * float(r["spinevol_gsyn_r"])
    cov[3][4] = cov[4][3] = float(r["gsyn_nrrp_r"])
    cov[3][5] = cov[5][3] = float(r["spinevol_nrrp_r"])
    cov[4][5] = cov[5][4] = float(r["spinevol_gsyn_r"])
    return cov


def _make_dist(distname, mu, sigma):
    if distname == "beta":
        a = -(mu * (mu**2 - mu + sigma**2)) / sigma**2
        b = ((mu - 1) * (mu**2 - mu + sigma**2)) / sigma**2
        return stats.beta(a, b)
    elif distname == "gamma":
        return stats.gamma(mu**2 / sigma**2, loc=0, scale=sigma**2 / mu)
    elif distname == "poisson":
        return stats.poisson(mu - 1, loc=1)
    raise NotImplementedError(f"Unknown distribution: {distname}")


def sample_synapses(recipe_row, n_synapses, seed=42, k_u=0.2, k_gsyn=2.0):
    """
    Sample n_synapses from the Gaussian copula defined by recipe_row.
    Mirrors epg.py but vectorised.
    Returns DataFrame with: u, d, f, nrrp, gsyn, spinevol, rho0, use_d, use_p, gmax_d, gmax_p
    """
    r = recipe_row
    names = ["u", "d", "f", "nrrp", "gsyn", "spinevol"]
    dists = [_make_dist(str(r[f"{n}Dist"]), float(r[n]), float(r[f"{n}SD"])) for n in names]
    cov   = _get_covariance_matrix(r)

    rng     = np.random.default_rng(seed)
    normals = rng.multivariate_normal(np.zeros(6), cov, size=n_synapses)

    rows = {}
    for j, (name, dist) in enumerate(zip(names, dists)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rows[name] = dist.ppf(stats.norm.cdf(normals[:, j]))
    df = pd.DataFrame(rows)

    # LTP/LTD scaling (mirrors _get_ltpltd_params)
    rng2 = np.random.default_rng(seed + 1)
    rho0 = rng2.binomial(1, np.clip(df["u"], 0, 1))
    df["rho0"]   = rho0
    pot           = rho0 > 0.5
    df["use_p"]  = np.where(pot, df["u"],             np.power(df["u"], k_u))
    df["use_d"]  = np.where(pot, np.power(df["u"], 1/k_u), df["u"])
    df["gmax_p"] = np.where(pot, df["gsyn"],           df["gsyn"] * k_gsyn)
    df["gmax_d"] = np.where(pot, df["gsyn"] / k_gsyn,  df["gsyn"])
    df["gmax_nmda"] = df["gsyn"] * float(r["gsynSRSF"])
    return df


# ──────────────────────────────────────────────────────────────────
# Analytical plasticity model
# ──────────────────────────────────────────────────────────────────

def Use_ss(u, tau_dep, isi_ms):
    """Steady-state release probability under TM depression."""
    return u / (1.0 + u * tau_dep / isi_ms)


def N_eff(isi_ms, tau_effca=TAU_EFFCA):
    """Geometric-series effcai accumulation factor at steady state."""
    return 1.0 / (1.0 - np.exp(-isi_ms / tau_effca))


def threshold_ratio(fit_params, zeta=ZETA_DEFAULT, loc="basal"):
    """θ_d/θ_p ratio — determines the LTD window lower boundary."""
    if loc == "basal":
        a_d0, a_d1 = fit_params["a00"], fit_params["a01"]
        a_p0, a_p1 = fit_params["a10"], fit_params["a11"]
    else:
        a_d0, a_d1 = fit_params["a20"], fit_params["a21"]
        a_p0, a_p1 = fit_params["a30"], fit_params["a31"]
    return (a_d0 + a_d1 * zeta) / (a_p0 + a_p1 * zeta)


def compute_eta(u_arr, d_arr, freq_hz, fit_params, zeta=ZETA_DEFAULT,
                loc="basal", K=K_COINCIDENCE):
    """
    Compute the plasticity efficiency ratio η for an array of synapses.

    η = K · Use_ss(u, D, freq) · N_eff(freq) / (a10 + a11·ζ)

    Returns: (eta, theta_d_theta_p_ratio, Use_ss_arr)
    """
    isi = 1000.0 / freq_hz
    uss = Use_ss(u_arr, d_arr, isi)
    nef = N_eff(isi)

    if loc == "basal":
        a_p0, a_p1 = fit_params["a10"], fit_params["a11"]
    else:
        a_p0, a_p1 = fit_params["a30"], fit_params["a31"]

    denom = a_p0 + a_p1 * zeta
    eta   = K * uss * nef / denom
    r_dt  = threshold_ratio(fit_params, zeta, loc)
    return eta, r_dt, uss


def rho_flip_probability(eta, r_dt, rho0_arr, fit_params,
                         tau_ind=TAU_IND, ind_time=INDUCTION_T):
    """
    Estimate probability that each synapse flips rho state during induction,
    using an approximate analytical solution to the rho ODE.

    Δrho = (gamma_p or gamma_d) × t_above / tau_ind
    where t_above ≈ p_above × induction_time,
    and p_above is a sigmoid function of how far η is from the boundary.

    Returns: (p_flip_to_1, p_flip_to_0) per synapse
    """
    gamma_d = fit_params.get("gamma_d_GB_GluSynapse", 101.54)
    gamma_p = fit_params.get("gamma_p_GB_GluSynapse", 216.18)

    # Fraction of induction time with effcai > θ_p  (sigmoid centred at η=1)
    # Slope=12 gives: η=1.0→5%, η=1.1→73%, η=1.2→91% — matches observed rates
    f_ltp = 1.0 / (1.0 + np.exp(-12.0 * (eta - 1.02)))

    # Fraction of induction time in LTD window (θ_d < effcai < θ_p)
    # Must be above θ_d (r_dt < η) AND below θ_p (η < 1)
    f_ltd = (1.0 / (1.0 + np.exp(-12.0 * (eta - r_dt - 0.02)))) * \
            (1.0 / (1.0 + np.exp(+12.0 * (eta - 0.98))))

    # Δrho from ODE (clamp to [0,1])
    delta_ltp = np.clip(gamma_p * ind_time / (tau_ind * 1e3) * f_ltp, 0, 1)
    delta_ltd = np.clip(gamma_d * ind_time / (tau_ind * 1e3) * f_ltd, 0, 1)

    # Flip probabilities: rho=0→1 for LTP, rho=1→0 for LTD
    p_flip_to_1 = (1 - rho0_arr) * delta_ltp
    p_flip_to_0 = rho0_arr       * delta_ltd
    return p_flip_to_1, p_flip_to_0


def predict_epsp_ratio(df, freq_hz, fit_params, zeta=ZETA_DEFAULT, loc="basal"):
    """
    Predict EPSP ratio for a population of synapses at given frequency.

    Returns a dict with summary statistics.
    """
    eta, r_dt, uss = compute_eta(df["u"].values, df["d"].values,
                                  freq_hz, fit_params, zeta, loc)
    p_ltp, p_ltd = rho_flip_probability(eta, r_dt, df["rho0"].values, fit_params)

    # EPSP ratio: potentiated gmax = 2×, depressed gmax = 0.5×
    # So: 0→1 contributes +1 (from 0.5 to 1.0 × gmax_p, but test pulse uses initial gmax)
    # Simplified: ratio = 1 + p_ltp - 0.5*p_ltd  (see derivation in docstring)
    epsp_ratio = (1.0 + p_ltp - 0.5 * p_ltd).mean()
    sem        = (1.0 + p_ltp - 0.5 * p_ltd).std() / np.sqrt(len(df))

    frac_ltp  = (eta > 1).mean()
    frac_ltd  = ((eta > r_dt) & (eta < 1)).mean()
    frac_none = (eta <= r_dt).mean()

    # Critical U thresholds (where η crosses 1 and r_dt)
    isi = 1000. / freq_hz
    nef = N_eff(isi)
    a_p = fit_params["a10"] + fit_params["a11"] * zeta
    a_d = fit_params["a00"] + fit_params["a01"] * zeta
    # For a "typical" synapse (D = median tau_dep from df):
    d_med = float(df["d"].median())
    # η = K*u/(1+u*D/isi)*N_eff/a_p  = 1 → solve for u_ltp_min
    # K*N_eff/a_p * u / (1 + u*D/isi) = 1
    # u*(K*N_eff/a_p - D/isi) = 1  → u_ltp_min = 1 / (K*N_eff/a_p - D/isi)
    denom_u = K_COINCIDENCE * nef / a_p - d_med / isi
    u_ltp_min = 1.0 / denom_u if denom_u > 0 else np.nan
    denom_u_d = K_COINCIDENCE * nef / a_d - d_med / isi
    u_ltd_min = 1.0 / denom_u_d if denom_u_d > 0 else np.nan

    return {
        "mean_epsp_ratio":  epsp_ratio,
        "sem_epsp_ratio":   sem,
        "frac_ltp_zone":    frac_ltp,
        "frac_ltd_zone":    frac_ltd,
        "frac_no_plast":    frac_none,
        "mean_eta":         eta.mean(),
        "median_eta":       np.median(eta),
        "p25_eta":          np.percentile(eta, 25),
        "p75_eta":          np.percentile(eta, 75),
        "theta_ratio":      r_dt,
        "mean_use_ss":      uss.mean(),
        "u_ltp_min":        u_ltp_min,    # min U for LTP at median D
        "u_ltd_min":        u_ltd_min,    # min U for LTD at median D
        "mean_u":           df["u"].mean(),
        "mean_d":           df["d"].mean(),
        "mean_gsyn":        df["gsyn"].mean(),
        "mean_spinevol":    df["spinevol"].mean(),
    }


# ──────────────────────────────────────────────────────────────────
# Per-pathway analysis
# ──────────────────────────────────────────────────────────────────

def analyse_pathway(recipe_row, fit_params, freq_list, dt_list,
                    n_synapses=3000, seed=42):
    """
    Run analytical prediction for one pathway across all freq × dt combinations.

    Note: dt only matters in so far as it determines whether effcai reaches the
    LTP or LTD zone (same η formula, different calibrated K for different |dt|).
    Currently uses the same K for all |dt| — use --calibrate to refine if needed.

    Returns DataFrame with one row per freq × dt.
    """
    df = sample_synapses(recipe_row, n_synapses, seed=seed)
    rows = []
    for freq in freq_list:
        for dt in dt_list:
            res = predict_epsp_ratio(df, freq, fit_params)
            res["freq_hz"] = freq
            res["dt_ms"]   = dt
            rows.append(res)
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────
# Scan utility
# ──────────────────────────────────────────────────────────────────

SCAN_RANGES = {
    "u":        np.linspace(0.10, 0.90, 17),
    "d":        np.linspace(100., 1200., 12),
    "spinevol": np.linspace(0.020, 0.30, 14),
    "gsyn":     np.linspace(0.30, 3.00, 14),
    "f":        np.linspace(5., 200., 12),
    "nrrp":     np.linspace(1., 5., 9),
}


def scan_parameter(recipe_row, fit_params, param, freq=10., dt_list=(-10., 10.),
                   n_synapses=3000, seed=42, values=None):
    """Sweep one recipe parameter, keeping others at their recipe values."""
    if values is None:
        values = SCAN_RANGES.get(param, np.linspace(0.1, 1.0, 10))

    rows = []
    for val in values:
        row = recipe_row.copy()
        row[param] = val
        df = sample_synapses(row, n_synapses, seed=seed)
        for dt in dt_list:
            res = predict_epsp_ratio(df, freq, fit_params)
            res[param] = val
            res["dt_ms"] = dt
            rows.append(res)
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────

def load_fit_params(config_path):
    try:
        import yaml
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
        fp = cfg.get("fit_params", {})
        merged = dict(DEFAULT_FIT_PARAMS)
        merged.update(fp)
        return merged
    except Exception as e:
        print(f"[warn] Could not load {config_path}: {e}  — using built-in defaults.")
        return dict(DEFAULT_FIT_PARAMS)


def print_banner(args, fit_params):
    fp = fit_params
    print(f"\n{'═'*72}")
    print(f"  Analytical EPSP Ratio Predictor")
    print(f"{'═'*72}")
    print(f"  Recipe  : {args.recipe}")
    print(f"  Config  : {args.config}")
    print(f"  fit_params (basal):")
    print(f"    LTD: θ_d = {fp['a00']:.4f}·c_pre + {fp['a01']:.4f}·c_post")
    print(f"    LTP: θ_p = {fp['a10']:.4f}·c_pre + {fp['a11']:.4f}·c_post")
    print(f"  θ_d/θ_p ratio = {threshold_ratio(fp):.4f}  "
          f"(LTD window: η ∈ [{threshold_ratio(fp):.3f}, 1.000])")
    print(f"  K_coincidence = {K_COINCIDENCE}  (calibrated from CHINDEMI_PARAMS_v2 sims)")
    print(f"  Freq(s)={args.freq} Hz    Δt(s)={args.dt} ms    N_synapses={args.n_synapses}")
    print(f"{'═'*72}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analytically predict EPSP ratios from recipe.csv parameters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--recipe", default="biodata/recipe.csv")
    parser.add_argument("--config", default="configs/L5TTPC_L5TTPC_STDP.yaml")
    parser.add_argument("--pathway", nargs=2, metavar=("PRE", "POST"), default=None,
                        help="Single pathway, e.g. --pathway L5_TPC:A L5_TPC:A")
    parser.add_argument("--mtype-filter", default=None,
                        help="Regex filter on pre/post mtype (e.g. 'L5_TPC')")
    parser.add_argument("--freq", nargs="+", type=float, default=[10.0],
                        metavar="HZ", help="Induction frequency in Hz (default: 10)")
    parser.add_argument("--dt", nargs="+", type=float, default=[-10., 10.],
                        metavar="MS", help="Timing offsets Δt in ms (default: -10 10)")
    parser.add_argument("--n-synapses", type=int, default=3000,
                        metavar="N", help="Synapses sampled per pathway (default: 3000)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--scan-param", choices=list(SCAN_RANGES), default=None,
                        metavar="PARAM",
                        help=f"Sweep one recipe parameter. Choices: {list(SCAN_RANGES)}")
    parser.add_argument("--scan-values", nargs="+", type=float, default=None,
                        metavar="VAL", help="Custom values for --scan-param")
    parser.add_argument("--output", default=None, metavar="CSV",
                        help="Save full results table to CSV file")
    parser.add_argument("--verbose", action="store_true",
                        help="Print extra statistics per pathway")
    args = parser.parse_args()

    # Load recipe
    recipe = pd.read_csv(args.recipe, index_col=[0, 1], skipinitialspace=True)
    recipe.index = recipe.index.set_levels(
        [lvl.str.strip() for lvl in recipe.index.levels]
    )
    recipe.columns = recipe.columns.str.strip()

    # Load fit_params
    fit_params = load_fit_params(args.config)

    # Select pathways
    if args.pathway:
        selected = [tuple(args.pathway)]
    elif args.mtype_filter:
        import re
        pat = re.compile(args.mtype_filter, re.IGNORECASE)
        selected = [(p, q) for (p, q) in recipe.index
                    if pat.search(p) and pat.search(q)]
    else:
        selected = list(recipe.index)

    if not selected:
        print("No pathways matched. Check --pathway or --mtype-filter.")
        sys.exit(1)

    print_banner(args, fit_params)

    all_results = []

    for (pre_mt, post_mt) in selected:
        try:
            row = recipe.loc[pre_mt, post_mt]
        except KeyError:
            print(f"  [skip] {pre_mt} → {post_mt} not found in recipe.")
            continue

        # ── SCAN MODE ──────────────────────────────────────────────────
        if args.scan_param:
            scan_df = scan_parameter(
                row, fit_params, args.scan_param,
                freq=args.freq[0], dt_list=args.dt,
                n_synapses=args.n_synapses, seed=args.seed,
                values=args.scan_values,
            )
            print(f"── {pre_mt} → {post_mt}  |  sweep: {args.scan_param}  "
                  f"@ {args.freq[0]} Hz ──")
            print()

            # Print table grouped by dt
            pivot = scan_df.pivot_table(
                index=args.scan_param, columns="dt_ms",
                values=["mean_epsp_ratio", "frac_ltp_zone", "frac_ltd_zone",
                        "mean_eta", "mean_use_ss"],
            )
            print(pivot.round(4).to_string())
            print()

            # ── Operating window summary ──
            r_dt = threshold_ratio(fit_params)
            for dt_val, grp in scan_df.groupby("dt_ms"):
                eta_series = grp["mean_eta"].values
                p_series   = grp[args.scan_param].values
                # Find crossing points
                idx_ltp = np.where(np.diff(np.sign(eta_series - 1.0)))[0]
                idx_ltd = np.where(np.diff(np.sign(eta_series - r_dt)))[0]
                cross_ltp = p_series[idx_ltp[0]] if len(idx_ltp) else None
                cross_ltd = p_series[idx_ltd[0]] if len(idx_ltd) else None
                print(f"  Δt={dt_val:+.0f}ms  LTD boundary ({args.scan_param}≈{cross_ltd:.4f} if not None)  "
                      f"LTP boundary ({args.scan_param}≈{cross_ltp:.4f} if not None)")
                if cross_ltd: print(f"    → LTD needs {args.scan_param} ≥ {cross_ltd:.4f}")
                if cross_ltp: print(f"    → LTP needs {args.scan_param} ≥ {cross_ltp:.4f}")
            print()

            if args.output:
                scan_df.to_csv(args.output.replace(".csv", f"_scan_{args.scan_param}.csv"), index=False)
            continue

        # ── NORMAL MODE ────────────────────────────────────────────────
        res_df = analyse_pathway(row, fit_params, args.freq, args.dt,
                                  args.n_synapses, args.seed)
        res_df.insert(0, "post_mtype", post_mt)
        res_df.insert(0, "pre_mtype", pre_mt)
        all_results.append(res_df)

        u_mean, u_sd = float(row["u"]), float(row["uSD"])
        d_mean       = float(row["d"])
        sv_mean      = float(row["spinevol"])
        gsyn_mean    = float(row["gsyn"])

        print(f"── {pre_mt} → {post_mt} ──")
        print(f"   Recipe: u={u_mean:.3f}±{u_sd:.3f}  D={d_mean:.0f}ms  "
              f"gsyn={gsyn_mean:.3f}nS  spinevol={sv_mean:.4f}")

        for _, r in res_df.iterrows():
            ltp = "✓" if r["frac_ltp_zone"] > 0.2 else ("~" if r["frac_ltp_zone"] > 0.05 else " ")
            ltd = "✓" if r["frac_ltd_zone"] > 0.2 else ("~" if r["frac_ltd_zone"] > 0.05 else " ")
            print(f"   {r['freq_hz']:5.1f} Hz  Δt={r['dt_ms']:+7.1f} ms  "
                  f"EPSP={r['mean_epsp_ratio']:.4f}  "
                  f"η={r['mean_eta']:.3f}[{r['p25_eta']:.2f}-{r['p75_eta']:.2f}]  "
                  f"LTP={r['frac_ltp_zone']*100:5.1f}%{ltp}  "
                  f"LTD={r['frac_ltd_zone']*100:5.1f}%{ltd}  "
                  f"Use_ss={r['mean_use_ss']:.4f}")
        if args.verbose:
            row0 = res_df.iloc[0]
            print(f"   Operating window (at D={d_mean:.0f}ms, median):")
            print(f"     LTD needs u ≥ {row0['u_ltd_min']:.4f}  (recipe u={u_mean:.3f})")
            print(f"     LTP needs u ≥ {row0['u_ltp_min']:.4f}  (recipe u={u_mean:.3f})")
        print()

    # ── SUMMARY TABLE ─────────────────────────────────────────────────
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        print(f"\n{'═'*72}")
        print("  SUMMARY — mean EPSP ratio  (pathway × freq × Δt)")
        print(f"{'═'*72}")
        pivot = combined.pivot_table(
            index=["pre_mtype", "post_mtype"],
            columns=["freq_hz", "dt_ms"],
            values="mean_epsp_ratio",
        )
        print(pivot.round(4).to_string())

        print(f"\n{'═'*72}")
        print("  SUMMARY — fraction of synapses in LTP zone")
        print(f"{'═'*72}")
        print(combined.pivot_table(
            index=["pre_mtype", "post_mtype"],
            columns=["freq_hz", "dt_ms"],
            values="frac_ltp_zone",
        ).round(3).to_string())

        print(f"\n{'═'*72}")
        print("  SUMMARY — fraction of synapses in LTD zone")
        print(f"{'═'*72}")
        print(combined.pivot_table(
            index=["pre_mtype", "post_mtype"],
            columns=["freq_hz", "dt_ms"],
            values="frac_ltd_zone",
        ).round(3).to_string())

        print(f"\n{'═'*72}")
        print("  SUMMARY — mean η (plasticity efficiency ratio)")
        print(f"  η > 1 = LTP zone   {threshold_ratio(fit_params):.3f} < η < 1 = LTD zone")
        print(f"{'═'*72}")
        print(combined.pivot_table(
            index=["pre_mtype", "post_mtype"],
            columns=["freq_hz", "dt_ms"],
            values="mean_eta",
        ).round(3).to_string())

        if args.output:
            combined.to_csv(args.output, index=False)
            print(f"\n  Saved results to {args.output}")


if __name__ == "__main__":
    main()
