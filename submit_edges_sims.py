#!/usr/bin/env python3
"""
Submit bluecellulab STDP induction simulations using dhuruva_modified_edges.h5.

theta_d_GB / theta_p_GB are injected from the precomputed edges HDF5 file;
all other synapse parameters (rho0_GB, Use_d/p, gmax_d/p) are loaded automatically
by bluecellulab from dhuruva_circuit_config.json.  No EPG / c_pre / c_post step.

Results are written as bluecellulab_results/rho.h5 (SONATA format) inside each workdir.

Usage:
    python submit_edges_sims.py --params chindemi --freq 10Hz --fastforward 280000
    python submit_edges_sims.py --params chindemi --freq 10Hz --fastforward 280000 --dry-run
    python submit_edges_sims.py --params chindemi --freq 10Hz --fastforward 280000 --skip-existing
    python submit_edges_sims.py --params chindemi --freq 10Hz --fastforward 280000 --execution-mode cpu --workers 4
"""

import argparse
import glob
import logging
import multiprocessing
import multiprocessing.pool
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────

PLASTYFIRE_ROOT  = os.path.dirname(os.path.abspath(__file__))
PAIRRUNNER       = os.path.join(PLASTYFIRE_ROOT, "plastyfire", "pairrunner_edges.py")
DEFAULT_RESULTS  = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
DEFAULT_EDGES_H5 = "/project/ctb-emuller/dhuruva/plastyfire/data/dhuruva_modified_edges.h5"
DEFAULT_TRACES_DIR = "/project/ctb-emuller/dhuruva/plastyfire/bluecellulab_results"
ACCOUNT          = "ctb-emuller"
WALLTIME         = "01:00:00"

# ── Parameter presets (global HOC params only — thresholds come from edges.h5) ──

PARAM_PRESETS = {
    "chindemi": {
        "tau_effca_GB_GluSynapse": 278.3177658387,
        "gamma_d_GB_GluSynapse":   101.5387594661,
        "gamma_p_GB_GluSynapse":   216.1841700668,
    },
}

# ── sbatch template ──────────────────────────────────────────────────────────

SBATCH_TMPL = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --time={walltime}
#SBATCH --chdir={workdir}
#SBATCH --output={log_file}

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR NEURODAMUS_DIR HOC_LIBRARY_PATH
hash -r

module --force purge
module load StdEnv/2023 scipy-stack/2024a gcc/12.3 openmpi/4.1.5 \\
    hdf5-mpi/1.14.2 cmake/3.31.0 mpi4py/4.0.3 python/3.11.5

export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages
export PYTHONPATH=$PYTHONPATH:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python3.11/site-packages
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
export LD_LIBRARY_PATH=$JPKGDIR/lib
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:{plastyfire_root}

python {pairrunner} {runner_args}
"""


# ── Helpers ──────────────────────────────────────────────────────────────────

def find_workdirs(results_dir, freq):
    """
    Find all simulation workdirs matching the given frequency string.
    Expected structure:
        <results_dir>/fitting/*/seed*/*_STDP/simulations/<pre>-<post>/<freq>_<dt>ms/
    Requires prefire_simulation_config.json + prefire_prespikes.h5.
    """
    pattern = os.path.join(
        results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*", f"{freq}_*ms",
    )
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs
            if os.path.isfile(os.path.join(d, "prefire_simulation_config.json"))
            and os.path.isfile(os.path.join(d, "prefire_prespikes.h5"))]


def build_runner_args(params_name, fastforward, edges_h5, traces_output_dir):
    parts = [f"--params {params_name}"]
    if fastforward is not None:
        parts.append(f"--fastforward {fastforward:.1f}")
    if edges_h5 != DEFAULT_EDGES_H5:
        parts.append(f"--edges-h5 {edges_h5}")
    if traces_output_dir is not None:
        parts.append(f"--traces-output-dir {traces_output_dir}")
    return " ".join(parts)


def write_sbatch(workdir, params_name, fastforward, edges_h5, traces_output_dir):
    pair  = os.path.basename(os.path.dirname(workdir))
    freq  = os.path.basename(workdir)
    name  = f"edges_{freq}_{pair}"
    content = SBATCH_TMPL.format(
        job_name=name,
        account=ACCOUNT,
        walltime=WALLTIME,
        workdir=workdir,
        log_file=f"{name}.log",
        pairrunner=PAIRRUNNER,
        runner_args=build_runner_args(params_name, fastforward, edges_h5, traces_output_dir),
        plastyfire_root=PLASTYFIRE_ROOT,
    )
    script = os.path.join(workdir, "edges_sim.batch")
    with open(script, "w") as f:
        f.write(content)
    os.chmod(script, 0o755)
    return script


class _NoDaemonProcess(multiprocessing.Process):
    """Process subclass whose daemon flag is permanently False so it can spawn children."""
    @property
    def daemon(self):
        return False
    @daemon.setter
    def daemon(self, value):
        pass  # pool sets w.daemon=True after creation; silently ignore it


class _NoDaemonPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super().Process(*args, **kwds)
        proc.__class__ = _NoDaemonProcess
        return proc


def _cpu_worker(args):
    """Top-level function (picklable) for multiprocessing.Pool."""
    workdir, fit_params, edges_h5, fastforward, traces_output_dir = args
    import plastyfire.simulator_edges as sim_mod
    import pickle, os, logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )
    log = logging.getLogger(__name__)
    log.info("CPU worker starting: %s", workdir)
    try:
        results = sim_mod.runconnectedpair_prefire_from_edges(
            workdir,
            fit_params=fit_params or None,
            edges_h5=edges_h5,
            fastforward=fastforward,
            traces_output_dir=traces_output_dir,
        )
        out_path = os.path.join(workdir, "simulation_edges.pkl")
        with open(out_path, "wb") as f:
            pickle.dump(results, f, protocol=-1)
        log.info("Done: %s → %s", workdir, out_path)
        return (workdir, True, None)
    except Exception as e:
        log.error("Failed %s: %s", workdir, e)
        return (workdir, False, str(e))


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Submit edges-based STDP prefire sims (slurm or local CPU)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--params", choices=list(PARAM_PRESETS.keys()), default="chindemi",
        help="Named parameter preset for global HOC vars (default: chindemi)",
    )
    parser.add_argument(
        "--freq", default="10Hz",
        help="Induction frequency prefix to filter workdirs (default: 10Hz)",
    )
    parser.add_argument(
        "--fastforward", type=float, default=None,
        help="Fast-forward point in ms (e.g. 280000)",
    )
    parser.add_argument(
        "--results-dir", default=DEFAULT_RESULTS,
        help=f"Root results directory (default: {DEFAULT_RESULTS})",
    )
    parser.add_argument(
        "--edges-h5", default=DEFAULT_EDGES_H5,
        help=f"Path to dhuruva_modified_edges.h5 (default: {DEFAULT_EDGES_H5})",
    )
    parser.add_argument(
        "--execution-mode", choices=["slurm", "cpu"], default="slurm",
        help="slurm: submit one sbatch per workdir; cpu: run locally in parallel (default: slurm)",
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel workers for --execution-mode cpu (default: 4)",
    )
    parser.add_argument(
        "--traces-output-dir", default=DEFAULT_TRACES_DIR,
        help=(
            "Directory where simulation_traces.pkl files are written in the structure "
            "<traces-output-dir>/<pair_name>/<protocol>/simulation_traces.pkl "
            f"(default: {DEFAULT_TRACES_DIR}). "
            "This is the layout plastyfitting/cicr_common.py expects. "
            "Pass 'none' to fall back to writing inside bluecellulab_results/ per workdir."
        ),
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip workdirs that already have both bluecellulab_results/rho.h5 "
             "and simulation_traces.pkl in the traces output dir",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Write batch scripts but print sbatch commands instead of running them "
             "(slurm mode); or list workdirs without running (cpu mode)",
    )
    args = parser.parse_args()

    traces_output_dir = None if args.traces_output_dir.lower() == "none" else args.traces_output_dir

    workdirs = find_workdirs(args.results_dir, args.freq)
    if not workdirs:
        logger.error("No workdirs found for freq=%s in %s", args.freq, args.results_dir)
        sys.exit(1)

    logger.info("Found %d workdirs for %s", len(workdirs), args.freq)

    if args.skip_existing:
        before = len(workdirs)
        def _is_done(d):
            has_rho = os.path.isfile(os.path.join(d, "bluecellulab_results", "rho.h5"))
            if traces_output_dir is not None:
                pair = os.path.basename(os.path.dirname(d))
                proto = os.path.basename(d)
                has_traces = os.path.isfile(
                    os.path.join(traces_output_dir, pair, proto, "simulation_traces.pkl")
                )
            else:
                has_traces = os.path.isfile(
                    os.path.join(d, "bluecellulab_results", "simulation_traces.pkl")
                )
            return has_rho and has_traces
        workdirs = [d for d in workdirs if not _is_done(d)]
        logger.info("  %d already done → %d remaining", before - len(workdirs), len(workdirs))

    if not workdirs:
        logger.info("All simulations already completed.")
        return

    if args.execution_mode == "slurm":
        submitted = failed = 0
        for workdir in workdirs:
            script = write_sbatch(workdir, args.params, args.fastforward, args.edges_h5, traces_output_dir)
            if args.dry_run:
                print(f"sbatch {script}")
                submitted += 1
            else:
                result = subprocess.run(["sbatch", script], capture_output=True, text=True)
                if result.returncode == 0:
                    logger.info("%s  (%s)", result.stdout.strip(), os.path.basename(workdir))
                    submitted += 1
                else:
                    logger.error("FAILED %s: %s", os.path.basename(workdir), result.stderr.strip())
                    failed += 1
        logger.info("Done: %d submitted, %d failed", submitted, failed)

    else:  # cpu
        fit_params = dict(PARAM_PRESETS[args.params]) if args.params else {}
        worker_args = [
            (workdir, fit_params, args.edges_h5, args.fastforward, traces_output_dir)
            for workdir in workdirs
        ]

        if args.dry_run:
            for workdir in workdirs:
                print(workdir)
            return

        n_workers = min(args.workers, len(workdirs))
        logger.info("Running %d sims with %d workers", len(workdirs), n_workers)

        succeeded = failed = 0
        with _NoDaemonPool(processes=n_workers) as pool:
            for workdir, ok, err in pool.imap_unordered(_cpu_worker, worker_args):
                if ok:
                    logger.info("OK  %s", os.path.basename(workdir))
                    succeeded += 1
                else:
                    logger.error("ERR %s: %s", os.path.basename(workdir), err)
                    failed += 1

        logger.info("Done: %d succeeded, %d failed", succeeded, failed)


if __name__ == "__main__":
    main()
