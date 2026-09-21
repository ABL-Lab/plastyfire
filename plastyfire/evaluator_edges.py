"""
BluePyOpt evaluator for the Graupner & Brunel model using dhuruva_modified_edges.h5.

Replicates evaluator.py but drives induction simulations via the prefire+edges pipeline
(simulator_edges.py / pairrunner_edges_fit.py) instead of the recipe-based EPG pipeline.

Key differences from the original Evaluator:
  - Simulations are dispatched to pairrunner_edges_fit.py (subprocess), not pairrunner.py
  - Result pkl files are named simulation_edges_{hash}.pkl
  - rho arrays come as flat float lists (initial_rho / final_rho) rather than a time-series matrix
  - tau_effca_GB_GluSynapse is fixed at FITTED_TAU and not included in the optimised params

EPSP ratio computation — basis technique:
  Instead of running bluecellulab test-pulse simulations, the EPSP ratio is computed
  analytically from a pre-generated basis CSV (basis_{pre}_{post}.csv).

  The basis CSV is produced once per pair by run_basis_pair_edges.py and contains the
  mean EPSP for:
    - all-depressed config: all rho=0
    - each singleton: one rho_i=1, rest=0
    - all-potentiated config: all rho=1  (validation only)

  Because gmax_AMPA = gmax_d + rho*(gmax_p - gmax_d), and EPSPs sum linearly across
  synapses, the EPSP for any continuous rho vector is:

      EPSP(ρ) = e₀ + Σᵢ [ ρᵢ × (eᵢ_singleton - e₀) ]

  where ρᵢ ∈ [0,1] are the continuous rho values from the induction simulation.
  No ephys pkl files are needed.
"""

import csv
import hashlib
import json
import logging
import multiprocessing as mp
import os
import pickle
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

from plastyfire.evaluator import (
    Evaluator,
    FITTED_TAU,
    CONFIGS_DIR,
)
import plastyfire.simulator_edges as sim_mod

logger = logging.getLogger(__name__)

# JSON file listing individuals to tag and log in detail when evaluated.
# Format: {"label": [param0, param1, ...], ...}
# e.g. {"chindemi": [101.5, 216.2, 1.002, 1.954, 1.159, 2.483, 1.127, 2.456, 5.236, 1.782]}
TAGGED_INDIVIDUALS_FILE = "tagged_individuals.json"

# Tolerance for matching a tagged individual to param_values
_TAG_MATCH_ATOL = 1e-6


def _load_tagged_individuals():
    """Load tagged individuals from JSON file, return dict {label: list_of_values}."""
    if not os.path.isfile(TAGGED_INDIVIDUALS_FILE):
        return {}
    try:
        with open(TAGGED_INDIVIDUALS_FILE) as f:
            return json.load(f)
    except Exception as e:
        logger.warning("Could not load %s: %s", TAGGED_INDIVIDUALS_FILE, e)
        return {}


def _match_tag(param_values, tagged):
    """Return label if param_values match a tagged individual, else None."""
    for label, values in tagged.items():
        if len(values) == len(param_values) and np.allclose(param_values, values, atol=_TAG_MATCH_ATOL):
            return label
    return None


# Parameters optimised in the edges fitting loop (tau_effca is fixed)
EDGES_FIT_PARAM_NAMES = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00", "a01",
    "a10", "a11",
    "a20", "a21",
    "a30", "a31",
]

# Reduced parameter set used when apical and basal share one set of a-params
# (modelfitter_edges.py --same-apical-basal). The GA then searches 6 values
# instead of 10, and a20/a21/a30/a31 are mirrored from a00/a01/a10/a11 at the
# single point where the vector is turned into a parameter dict.
EDGES_FIT_PARAM_NAMES_TIED = [
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00", "a01",
    "a10", "a11",
]

# Set by modelfitter_edges.py before the run starts. Module-level because the
# multiprocessing workers import this module fresh and must agree on the
# vector layout; it is read, never mutated, once the run begins.
SAME_APICAL_BASAL = False


def expand_tied_params(param_values):
    """Map a 6-value tied vector to the full 10-value parameter vector.

    apical (a2x/a3x) is set equal to basal (a0x/a1x):
        a20 = a00, a21 = a01, a30 = a10, a31 = a11

    A vector that is already full length is returned unchanged, so callers can
    apply this unconditionally.
    """
    if not SAME_APICAL_BASAL or len(param_values) != len(EDGES_FIT_PARAM_NAMES_TIED):
        return list(param_values)
    gd, gp, a00, a01, a10, a11 = param_values
    return [gd, gp, a00, a01, a10, a11, a00, a01, a10, a11]


def active_param_names():
    """Names matching the vector length the GA is currently searching."""
    return EDGES_FIT_PARAM_NAMES_TIED if SAME_APICAL_BASAL else EDGES_FIT_PARAM_NAMES

# Direction constraints: True = LTP (mean EPSP ratio should be > 1),
#                        False = LTD (mean EPSP ratio should be < 1).
# Penalty weight per unit of violation (in the same units as the weighted error).
PROTOCOL_DIRECTION = {
    "mrk97_07": True,   # dt = +10 ms → LTP → ratio must exceed 1
    "mrk97_08": False,  # dt = -10 ms → LTD → ratio must be below 1
}
DIRECTION_PENALTY_SCALE = 10.0


# ---------------------------------------------------------------------------
# Module-level helpers (must be top-level for multiprocessing pickling)
# ---------------------------------------------------------------------------

def _load_basis(basis_dir, pre_gid, post_gid):
    """Load the basis CSV for a pair, returning a DataFrame."""
    csv_path = os.path.join(basis_dir, f"basis_{pre_gid}_{post_gid}.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Basis CSV not found: {csv_path}\n"
            f"Run run_basis_pair_edges.py --pre-gid {pre_gid} --post-gid {post_gid} first."
        )
    return pd.read_csv(csv_path)


def compute_epsp_from_basis(basis_df, rho_vec):
    """
    Compute EPSP mean and std for a rho vector using linear superposition with
    binary rho thresholding (rho >= 0.5 → 1, else 0) to match bluecellulab BCL.

    Mean:  EPSP(ρ) = e₀ + Σᵢ [ ρᵢ_bin × (eᵢ_singleton - e₀) ]
    Std:   Var(EPSP) = (1-k)² × Var(e₀) + Σᵢ_active Var(eᵢ_singleton)
           where k = number of potentiated synapses.

    :param basis_df: DataFrame from basis_{pre}_{post}.csv (must have 'mean' and 'std' cols)
    :param rho_vec: list/array of floats in [0, 1], one per synapse
    :returns: (epsp_mean, epsp_std) tuple
    """
    n = len(rho_vec)
    rho_bin = [1 if r >= 0.5 else 0 for r in rho_vec]
    all_zeros_key = ",".join(["0"] * n)

    row0 = basis_df.loc[basis_df["config"] == all_zeros_key]
    if row0.empty:
        raise ValueError(f"Basis CSV missing all-zeros config: {all_zeros_key!r}")
    e0     = float(row0["mean"].values[0])
    var0   = float(row0["std"].values[0]) ** 2

    epsp_mean = e0
    variance  = 0.0
    k = 0  # number of potentiated synapses

    for i, rho_i in enumerate(rho_bin):
        if rho_i == 0:
            continue
        singleton = ["0"] * n
        singleton[i] = "1"
        key = ",".join(singleton)
        row_i = basis_df.loc[basis_df["config"] == key]
        if row_i.empty:
            raise ValueError(f"Basis CSV missing singleton config: {key!r}")
        ei_mean = float(row_i["mean"].values[0])
        ei_std  = float(row_i["std"].values[0])
        epsp_mean += (ei_mean - e0)
        variance  += ei_std ** 2
        k += 1

    # Depressed baseline contribution to variance
    variance += (1 - k) ** 2 * var0
    epsp_std = float(np.sqrt(max(0.0, variance)))

    return epsp_mean, epsp_std


def compute_epsp_ratio_basis_batch(param_values, sim_dict, workdir, basis_dir):
    """
    Read simulation_edges_{hash}.pkl and compute the EPSP ratio using the basis technique.

    Matches bluecellulab BCL computation:
      - rho is binarized (>= 0.5 → 1, else 0) before EPSP evaluation
      - delta-method bias correction: ratio = (after/before) * (1 + CV²), CV² capped at 0.25

    Returns (protocol_id, epsp_ratio, pre_gid, post_gid) on success or (None, None, None, None) on failure.
    """
    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
    pkl_file = os.path.join(workdir, f"simulation_edges_{param_hash}.pkl")

    if not os.path.exists(pkl_file):
        logger.error("Edges simulation pkl not found: %s", pkl_file)
        return None, None, None, None

    try:
        with open(pkl_file, "rb") as f:
            raw = pickle.load(f)
    except (pickle.UnpicklingError, EOFError, KeyError, ValueError) as e:
        logger.warning("Corrupted edges pkl %s: %s — deleting", pkl_file, e)
        os.remove(pkl_file)
        return None, None, None, None

    # Use continuous rho values directly (no 0.5 threshold)
    initial_rho = list(raw["initial_rho"])
    final_rho   = list(raw["final_rho"])

    pre_gid, post_gid = workdir.split("/")[-2].split("-")

    try:
        basis_df = _load_basis(basis_dir, pre_gid, post_gid)
    except FileNotFoundError as e:
        logger.error("%s", e)
        return None, None, None, None

    try:
        epsp_before_mean, epsp_before_std = compute_epsp_from_basis(basis_df, initial_rho)
        epsp_after_mean,  _               = compute_epsp_from_basis(basis_df, final_rho)
    except ValueError as e:
        logger.error("Basis computation failed for %s: %s", sim_dict["protocol_id"], e)
        return None, None, None, None

    if epsp_before_mean == 0:
        logger.error("EPSP before induction is zero for %s", sim_dict["protocol_id"])
        return None, None, None, None

    # Delta-method bias correction: (after/before) * (1 + CV²), CV² capped at 0.25
    # Matches bluecellulab BCL computation (plot_compare_ndamus_bluecellulab.py)
    cv2 = min((epsp_before_std / epsp_before_mean) ** 2, 0.25) if epsp_before_mean > 0 else 0.0
    epsp_ratio = (epsp_after_mean / epsp_before_mean) * (1.0 + cv2)
    logger.info("EPSP ratio = %.4f (before=%.4f±%.4f, after=%.4f, cv2=%.4f) for protocol %s",
                epsp_ratio, epsp_before_mean, epsp_before_std, epsp_after_mean, cv2,
                sim_dict["protocol_id"])
    return sim_dict["protocol_id"], epsp_ratio, pre_gid, post_gid


def run_simulation_worker_edges(args):
    """
    Pool worker: run one induction simulation via pairrunner_edges_fit.py (subprocess).

    args tuple: (param_values, sim_dict, param_hash, edges_h5, cpre_cpost_cache, basis_dir, circuit_config)

    Using subprocess (not multiprocessing.Process) means this function may safely be
    called from within a Pool worker (daemon process) — the OS-level fork+exec of the
    subprocess is not restricted by the daemon flag.
    """
    param_values, sim_dict, param_hash, edges_h5, cpre_cpost_cache, basis_dir, circuit_config = args

    # Build per-parameter CLI arguments. Under --same-apical-basal the GA vector
    # carries only the 6 free values; mirror them into the apical slots here so
    # pairrunner_edges_fit.py still receives all 10 it requires.
    param_dict = dict(zip(EDGES_FIT_PARAM_NAMES, expand_tied_params(param_values)))
    param_args = [f"--{name}={value}" for name, value in param_dict.items()]
    param_args.append(f"--tau_effca_GB_GluSynapse={FITTED_TAU}")
    param_args.append(f"--fastforward={sim_dict['fastforward']}")
    param_args.append(f"--edges-h5={edges_h5}")
    param_args.append(f"--param_hash={param_hash}")
    if cpre_cpost_cache:
        param_args.append(f"--cpre-cpost-cache={cpre_cpost_cache}")
    if circuit_config:
        param_args.append(f"--circuit-config={circuit_config}")

    script_path = os.path.join(os.path.dirname(__file__), "pairrunner_edges_fit.py")
    cmd = [sys.executable, script_path] + param_args

    logger.info("Running edges simulation for %s in %s", sim_dict["protocol_id"], sim_dict["workdir"])
    try:
        result = subprocess.run(
            cmd, cwd=sim_dict["workdir"], capture_output=True, text=True, timeout=1800
        )
    except subprocess.TimeoutExpired:
        logger.error("Simulation timed out for %s", sim_dict["protocol_id"])
        return None

    if result.returncode != 0:
        logger.error("Simulation failed for %s:\n%s", sim_dict["protocol_id"], result.stderr[-2000:])
        return None

    pkl_file = os.path.join(sim_dict["workdir"], f"simulation_edges_{param_hash}.pkl")
    if not os.path.exists(pkl_file):
        logger.error("PKL file not created for %s", sim_dict["protocol_id"])
        return None

    return compute_epsp_ratio_basis_batch(param_values, sim_dict, sim_dict["workdir"], basis_dir)


# ---------------------------------------------------------------------------
# EvaluatorEdges class
# ---------------------------------------------------------------------------

class EvaluatorEdges(Evaluator):
    """
    Graupner & Brunel plasticity model evaluator using dhuruva_modified_edges.h5.

    Inherits all setup logic (__init__, objectives, all_sims, scoring) from Evaluator
    and overrides only the simulation dispatch + result parsing to use the edges pipeline.
    """

    def __init__(
        self,
        fit_params,
        invitro_db,
        seed,
        sample_size,
        ipp_id,
        work_dir=None,
        max_jobs=900,
        use_multiprocessing=True,
        max_workers=None,
        fitness_schedule_file=None,
        edges_h5=None,
        cpre_cpost_cache=None,
        basis_dir=None,
        circuit_config=None,
    ):
        """
        :param edges_h5: path to dhuruva_modified_edges.h5 (default: sim_mod.EDGES_H5_DEFAULT)
        :param cpre_cpost_cache: optional precomputed c_pre/c_post cache pkl
        :param basis_dir: directory containing basis_{pre}_{post}.csv files produced by
            run_basis_pair_edges.py (e.g. ``basis_results_edges_mini/``). Required.
        :param circuit_config: optional path to a circuit config JSON to override the
            'network' field in each prefire_simulation_config.json (e.g. for new ion channels).
        All other parameters mirror Evaluator.__init__.
        """
        super().__init__(
            fit_params=fit_params,
            invitro_db=invitro_db,
            seed=seed,
            sample_size=sample_size,
            ipp_id=ipp_id,
            work_dir=work_dir,
            max_jobs=max_jobs,
            use_multiprocessing=use_multiprocessing,
            max_workers=max_workers,
            fitness_schedule_file=fitness_schedule_file,
        )
        self.edges_h5 = edges_h5 or sim_mod.EDGES_H5_DEFAULT
        # Always store as absolute path: the cache pkl is passed as a CLI arg to
        # pairrunner_edges_fit.py which runs with cwd=workdir (a different directory).
        self.cpre_cpost_cache = os.path.abspath(cpre_cpost_cache) if cpre_cpost_cache else None
        if basis_dir is None:
            raise ValueError("basis_dir is required — point it at the directory containing basis_{pre}_{post}.csv files")
        self.basis_dir = os.path.abspath(basis_dir)
        self.circuit_config = os.path.abspath(circuit_config) if circuit_config else None
        logger.info("EvaluatorEdges: edges_h5  = %s", self.edges_h5)
        logger.info("EvaluatorEdges: basis_dir = %s", self.basis_dir)
        if self.cpre_cpost_cache:
            logger.info("EvaluatorEdges: cpre_cpost_cache = %s", self.cpre_cpost_cache)
        if self.circuit_config:
            logger.info("EvaluatorEdges: circuit_config   = %s", self.circuit_config)

    # ------------------------------------------------------------------
    # Override: direction-violation penalty
    # ------------------------------------------------------------------

    def _process_evaluation_results(self, results, param_values):
        """
        Add a direction-violation penalty on top of the parent's error.

        After calling the parent (which computes |mean_invitro - mean_insilico| / sem × weight),
        we add DIRECTION_PENALTY_SCALE × violation for each protocol whose simulated mean
        goes the wrong way:
          - mrk97_07 (LTP, dt=+10ms): mean EPSP ratio must be > 1
          - mrk97_08 (LTD, dt=-10ms): mean EPSP ratio must be < 1
        If the direction is correct, no penalty is added.
        """
        error = super()._process_evaluation_results(results, param_values)

        if not results:
            return error

        res_db = pd.DataFrame(results, columns=["protocol_id", "epsp_ratio"])
        means = res_db.groupby("protocol_id")["epsp_ratio"].mean()

        for i, obj in enumerate(self.objectives):
            expected_ltp = PROTOCOL_DIRECTION.get(obj.name)
            if expected_ltp is None:
                continue
            mean = means.get(obj.name)
            if mean is None:
                continue
            if expected_ltp and mean < 1.0:
                penalty = (1.0 - mean) * DIRECTION_PENALTY_SCALE
                logger.info(
                    "Direction penalty for %s: LTP mean=%.4f < 1.0, penalty=%.4f",
                    obj.name, mean, penalty,
                )
                error[i] += penalty
            elif not expected_ltp and mean > 1.0:
                penalty = (mean - 1.0) * DIRECTION_PENALTY_SCALE
                logger.info(
                    "Direction penalty for %s: LTD mean=%.4f > 1.0, penalty=%.4f",
                    obj.name, mean, penalty,
                )
                error[i] += penalty

        return error

    # ------------------------------------------------------------------
    # Override: multiprocessing-based evaluation using the edges pipeline
    # ------------------------------------------------------------------

    def evaluate_with_multiprocessing(self, param_values):
        """Evaluate one individual using the edges-based prefire simulation pipeline."""
        try:
            logger.info("EvaluatorEdges: evaluating individual %s", param_values)

            # Check if this individual is tagged before hitting the cache
            tagged = _load_tagged_individuals()
            tag_label = _match_tag(param_values, tagged)
            if tag_label:
                logger.info("EvaluatorEdges: tagged individual detected — label='%s'", tag_label)

            # Cache check (reuse parent's format) — skip cache for tagged individuals
            # so we always get fresh per-pair detail for logging.
            cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
            pklf_name = os.path.join(".cache", f"{cachekey}.pkl")
            if not tag_label and os.path.isfile(pklf_name):
                with open(pklf_name, "rb") as f:
                    cache_data = pickle.load(f)
                np.testing.assert_array_equal(param_values, cache_data["individual"])
                if "resdb" in cache_data:
                    res_db = cache_data["resdb"]
                    current_ids = [obj.name for obj in self.objectives]
                    filtered = res_db[res_db["protocol_id"].isin(current_ids)]
                    if len(filtered) == len(self.objectives):
                        results = [(r["protocol_id"], r["epsp_ratio"]) for _, r in filtered.iterrows()]
                        return self._process_evaluation_results(results, param_values)
                elif len(cache_data.get("error", [])) == len(self.objectives):
                    return cache_data["error"]

            param_hash = cachekey[:12]
            existing_results, existing_detailed, jobs_to_run = [], [], []

            for sim_dict in self.all_sims:
                pkl_file = os.path.join(sim_dict["workdir"], f"simulation_edges_{param_hash}.pkl")
                file_key = f"{pkl_file}_{param_hash}"

                if file_key in self._processed_files:
                    continue

                if os.path.exists(pkl_file):
                    self._processed_files.add(file_key)
                    try:
                        r = compute_epsp_ratio_basis_batch(
                            param_values, sim_dict, sim_dict["workdir"], self.basis_dir
                        )
                        if r is not None and r[1] is not None:
                            existing_results.append((r[0], r[1]))
                            existing_detailed.append(r)  # (protocol_id, epsp_ratio, pre_gid, post_gid)
                            continue
                    except Exception as e:
                        logger.info("Failed to load existing edges pkl %s: %s", pkl_file, e)

                jobs_to_run.append((
                    param_values, sim_dict, param_hash,
                    self.edges_h5, self.cpre_cpost_cache, self.basis_dir,
                    self.circuit_config,
                ))

            results = list(existing_results)
            detailed_results = list(existing_detailed)  # (protocol_id, epsp_ratio, pre_gid, post_gid)
            if jobs_to_run:
                # run_simulation_worker_edges uses subprocess.run (OS-level fork+exec),
                # which is safe from daemon workers. Use ThreadPoolExecutor so this works
                # whether called from a daemon Pool worker or the main process — threads
                # have no daemon-process child restriction.
                n_workers = self.max_workers or min(len(jobs_to_run), mp.cpu_count())
                logger.info(
                    "Running %d edges simulations with %d threads",
                    len(jobs_to_run), n_workers
                )
                with ThreadPoolExecutor(max_workers=n_workers) as tex:
                    futures = {tex.submit(run_simulation_worker_edges, job): job for job in jobs_to_run}
                    for fut in as_completed(futures):
                        r = fut.result()
                        if r is not None and len(r) >= 2 and r[0] is not None and r[1] is not None:
                            results.append((r[0], r[1]))
                            detailed_results.append(r)  # (protocol_id, epsp_ratio, pre_gid, post_gid)

            # Write tagged CSV if this is a labelled individual
            if tag_label and detailed_results:
                self._write_tagged_csv(tag_label, detailed_results, param_values)

            error = self._process_evaluation_results(results, param_values)
            logger.info("EvaluatorEdges evaluation complete")
            return error

        except Exception:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))

    def _write_tagged_csv(self, tag_label, detailed_results, param_values):
        """
        Write per-pair EPSP ratios for a tagged individual to a CSV in logs/.

        CSV columns: tag, generation, pre_gid, post_gid, protocol_id, epsp_ratio
        Also logs a summary table to the logger so it's visible in the run log.

        :param tag_label: str label from tagged_individuals.json
        :param detailed_results: list of (protocol_id, epsp_ratio, pre_gid, post_gid)
        :param param_values: parameter vector for this individual
        """
        os.makedirs("logs", exist_ok=True)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        csv_path = os.path.join("logs", f"tagged_{tag_label}_gen{self.current_generation}_{timestamp}.csv")

        # Log the full 10 params even when only 6 are free, so tagged CSVs stay
        # comparable across tied and untied runs.
        param_values = expand_tied_params(param_values)
        param_names = EDGES_FIT_PARAM_NAMES
        rows = []
        for entry in detailed_results:
            protocol_id, epsp_ratio, pre_gid, post_gid = entry
            rows.append({
                "tag": tag_label,
                "generation": self.current_generation,
                "pre_gid": pre_gid,
                "post_gid": post_gid,
                "protocol_id": protocol_id,
                "epsp_ratio": epsp_ratio,
            })

        fieldnames = ["tag", "generation", "pre_gid", "post_gid", "protocol_id", "epsp_ratio"]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        logger.info(
            "=== TAGGED INDIVIDUAL '%s' (gen %d) — per-pair EPSP ratios written to: %s ===",
            tag_label, self.current_generation, csv_path,
        )
        logger.info("  Parameters: %s", dict(zip(param_names, param_values)))

        # Summary by protocol
        df = pd.DataFrame(rows)
        summary = (
            df.groupby("protocol_id")["epsp_ratio"]
            .agg(["mean", "sem", "count"])
            .rename(columns={"mean": "mean_epsp_ratio", "sem": "sem_epsp_ratio", "count": "n_pairs"})
            .reset_index()
        )
        logger.info("  Per-protocol summary:\n%s", summary.to_string(index=False))
