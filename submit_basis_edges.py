#!/usr/bin/env python3
"""
Submit EPSP basis simulations using dhuruva_modified_edges.h5.

Finds all unique pairs in the refitting_results tree, then submits one job per
pair via sbatch (default) or runs them locally in parallel (--execution-mode cpu).

Each job calls run_basis_pair_edges.py, which:
  - Reads gmax_d/p_AMPA and Use_d/p from bluecellulab (auto-loaded from edges.h5)
  - Runs 10 test-pulse trials per rho configuration (all-0, singletons, all-1)
  - Saves basis_{pre}_{post}.csv to --output-dir

Usage:
    python submit_basis_edges.py
    python submit_basis_edges.py --dry-run
    python submit_basis_edges.py --skip-existing
    python submit_basis_edges.py --execution-mode cpu --workers 4
    python submit_basis_edges.py --output-dir /path/to/basis_results_edges
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

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PLASTYFIRE_ROOT  = os.path.dirname(os.path.abspath(__file__))
PAIR_RUNNER      = os.path.join(PLASTYFIRE_ROOT, "run_basis_pair_edges.py")
DEFAULT_RESULTS  = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
DEFAULT_OUT_DIR  = os.path.join(PLASTYFIRE_ROOT, "basis_results_edges_mini")
ACCOUNT          = "ctb-emuller"
WALLTIME         = "06:00:00"
NUM_TRIALS       = 5
CPUS_PER_TASK    = 12   # one per trial config; 10 trials + overhead

SBATCH_TMPL = """\
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --cpus-per-task={cpus}
#SBATCH --ntasks=1
#SBATCH --mem=24g
#SBATCH --time={walltime}
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
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin

python {runner} \\
    --pre-gid {pre_gid} \\
    --post-gid {post_gid} \\
    --sim-config {sim_config} \\
    --output-csv {output_csv} \\
    --num-trials {num_trials} \\
    --workers {workers}
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_pairs(results_dir):
    """
    Return {pair_str: sim_config_path} for all unique pairs that have
    simulation_config.json + prespikes.h5 in at least one frequency subdir.
    Uses the first matching subdir found for each pair.
    """
    pattern = os.path.join(
        results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*",
    )
    pair_dirs = sorted(glob.glob(pattern))

    pairs = {}
    for pair_dir in pair_dirs:
        pair = os.path.basename(pair_dir)
        if pair in pairs:
            continue
        # Find any freq subdir with the required files
        for sub in sorted(os.listdir(pair_dir)):
            sub_path = os.path.join(pair_dir, sub)
            cfg  = os.path.join(sub_path, "simulation_config.json")
            spk  = os.path.join(sub_path, "prespikes.h5")
            if os.path.isfile(cfg) and os.path.isfile(spk):
                pairs[pair] = cfg
                break
    return pairs


def write_sbatch(pair, pre_gid, post_gid, sim_config, output_csv,
                 out_dir, num_trials, walltime):
    name    = f"basis_{pair}"
    log     = os.path.join(out_dir, "logs", f"{name}.log")
    os.makedirs(os.path.dirname(log), exist_ok=True)
    content = SBATCH_TMPL.format(
        job_name=name,
        account=ACCOUNT,
        cpus=CPUS_PER_TASK,
        walltime=walltime,
        log_file=log,
        runner=PAIR_RUNNER,
        pre_gid=pre_gid,
        post_gid=post_gid,
        sim_config=sim_config,
        output_csv=output_csv,
        num_trials=num_trials,
        workers=CPUS_PER_TASK,
    )
    script = os.path.join(out_dir, "scripts", f"{name}.batch")
    os.makedirs(os.path.dirname(script), exist_ok=True)
    with open(script, "w") as f:
        f.write(content)
    os.chmod(script, 0o755)
    return script


class _NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False
    @daemon.setter
    def daemon(self, value):
        pass


class _NoDaemonPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super().Process(*args, **kwds)
        proc.__class__ = _NoDaemonProcess
        return proc


def _cpu_worker(args):
    pair, pre_gid, post_gid, sim_config, output_csv, num_trials, workers = args
    import subprocess, logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger(__name__)
    log.info("Starting %s", pair)
    cmd = [
        sys.executable, PAIR_RUNNER,
        "--pre-gid",    str(pre_gid),
        "--post-gid",   str(post_gid),
        "--sim-config", sim_config,
        "--output-csv", output_csv,
        "--num-trials", str(num_trials),
        "--workers",    str(workers),
    ]
    try:
        subprocess.check_call(cmd)
        log.info("Done %s", pair)
        return (pair, True, None)
    except subprocess.CalledProcessError as e:
        log.error("Failed %s: %s", pair, e)
        return (pair, False, str(e))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit EPSP basis simulations using dhuruva_modified_edges.h5",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS,
                        help=f"Root refitting results directory (default: {DEFAULT_RESULTS})")
    parser.add_argument("--output-dir",  default=DEFAULT_OUT_DIR,
                        help=f"Directory for basis CSV output (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--num-trials",  type=int, default=NUM_TRIALS,
                        help=f"Test-pulse trials per rho config (default: {NUM_TRIALS})")
    parser.add_argument("--walltime",    default=WALLTIME,
                        help=f"Slurm wall time (default: {WALLTIME})")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm",
                        help="slurm: one sbatch per pair; cpu: local parallel (default: slurm)")
    parser.add_argument("--workers",     type=int, default=4,
                        help="Parallel workers for --execution-mode cpu (default: 4)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip pairs whose output CSV already exists")
    parser.add_argument("--dry-run",     action="store_true",
                        help="Print commands/workdirs without submitting")
    args = parser.parse_args()

    pairs = find_pairs(args.results_dir)
    if not pairs:
        logger.error("No pairs found in %s", args.results_dir)
        sys.exit(1)

    logger.info("Found %d unique pairs", len(pairs))
    os.makedirs(args.output_dir, exist_ok=True)

    if args.skip_existing:
        before = len(pairs)
        pairs  = {p: v for p, v in pairs.items()
                  if not os.path.isfile(
                      os.path.join(args.output_dir, f"basis_{p.replace('-', '_')}.csv")
                  )}
        logger.info("  %d already done → %d remaining", before - len(pairs), len(pairs))

    if not pairs:
        logger.info("All basis results already exist.")
        return

    if args.execution_mode == "slurm":
        submitted = failed = 0
        for pair, sim_config in sorted(pairs.items()):
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            output_csv = os.path.join(args.output_dir, f"basis_{pre_gid}_{post_gid}.csv")
            script = write_sbatch(
                pair, pre_gid, post_gid, sim_config, output_csv,
                args.output_dir, args.num_trials, args.walltime,
            )
            if args.dry_run:
                print(f"sbatch {script}")
                submitted += 1
            else:
                result = subprocess.run(["sbatch", script], capture_output=True, text=True)
                if result.returncode == 0:
                    logger.info("%s  (%s)", result.stdout.strip(), pair)
                    submitted += 1
                else:
                    logger.error("FAILED %s: %s", pair, result.stderr.strip())
                    failed += 1
        logger.info("Done: %d submitted, %d failed", submitted, failed)

    else:  # cpu
        worker_args = []
        for pair, sim_config in sorted(pairs.items()):
            pre_gid, post_gid = [int(x) for x in pair.split("-")]
            output_csv = os.path.join(args.output_dir, f"basis_{pre_gid}_{post_gid}.csv")
            worker_args.append((
                pair, pre_gid, post_gid, sim_config, output_csv,
                args.num_trials, CPUS_PER_TASK,
            ))

        if args.dry_run:
            for wa in worker_args:
                print(wa[0], wa[3])
            return

        n_workers = min(args.workers, len(worker_args))
        logger.info("Running %d pairs with %d pool workers", len(worker_args), n_workers)

        succeeded = failed = 0
        with _NoDaemonPool(processes=n_workers) as pool:
            for pair, ok, err in pool.imap_unordered(_cpu_worker, worker_args):
                if ok:
                    logger.info("OK  %s", pair)
                    succeeded += 1
                else:
                    logger.error("ERR %s: %s", pair, err)
                    failed += 1

        logger.info("Done: %d succeeded, %d failed", succeeded, failed)


if __name__ == "__main__":
    main()
