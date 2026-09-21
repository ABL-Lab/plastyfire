#!/usr/bin/env python3
"""
Generate prefire_simulation_config_defit2.json and neurodamus_sbatch_defit2.sh
for every workdir, pointing neurodamus at the ion-channel circuit config and the
DE fit #2 best params.  Writes to out_defit2/ so neither the baseline out/ nor
the earlier (cancelled) DE fit #1 out_defit/ results are touched.

Changes vs the baseline prefire_simulation_config.json:
  - "network"  → dhuruva_modified_ion_channels_circuit_config.json
  - output_dir → {workdir}/out_defit2
  - GluSynapse.gamma_d_GB → 77.7558   (DE fit #2)
  - GluSynapse.gamma_p_GB → 299.9121  (DE fit #2)
  - GluSynapse.tau_effca_GB stays 278.3177658387

theta_d / theta_p are read from dhuruva_modified_edges.h5, so the DE fit #2
a-params MUST be injected first (see check_thresholds below, which refuses to
submit otherwise):

    python compute_thresholds_from_cache.py \\
        --cache cpre_cpost_cache/ion_channels_tau278.pkl \\
        --output-dir threshold_results_de_fit2/ \\
        --a00 1.003498 --a01 2.902478 --a10 1.644558 --a11 2.764812 \\
        --a20 1.003498 --a21 2.902478 --a30 1.644558 --a31 2.764812 \\
        --edges-h5 data/dhuruva_modified_edges.h5
    python inject_thresholds.py --results-dir threshold_results_de_fit2/ \\
        --edges-h5 data/dhuruva_modified_edges.h5 --force

Usage:
    python submit_neurodamus_defit.py --dry-run
    python submit_neurodamus_defit.py
    python submit_neurodamus_defit.py --skip-existing
"""

import argparse
import copy
import glob
import json
import logging
import os
import subprocess
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RESULTS = os.path.join(PLASTYFIRE_ROOT, "refitting_results")
IC_CIRCUIT_CONFIG = os.path.join(PLASTYFIRE_ROOT, "data", "dhuruva_modified_ion_channels_circuit_config.json")
ACCOUNT = "ctb-emuller"
# Suffix distinguishing this run's config / sbatch / output dir / job names from
# both the baseline and the cancelled DE fit #1 run (tag "defit").
TAG = "defit2"

# DE fit #2 from analytical_method/fit.py, refit at 1 ms with the physiological
# gamma bounds gamma_d in [50, 200] and gamma_p in [150, 300]. Apical tied to
# basal. Confirmed in bluecellulab as hash b8c7ff3ecf0a over 658 pairs:
#   dt=-10 0.8065 | +5 1.2376 | +10 1.1549   (in vitro 0.7922 / 1.2038 / 1.2013)
# Supersedes DE fit #1 (gamma_p 483.4465), whose neurodamus jobs were cancelled.
#
# Neurodamus does NOT have bluecellulab's runtime a-param path: it reads
# theta_d / theta_p straight out of dhuruva_modified_edges.h5. So the DE-fit
# thresholds MUST be injected with compute_thresholds_from_cache.py +
# inject_thresholds.py before this runs, or it silently uses whatever thetas
# happen to be in the edge file.
DE_FIT_PARAMS = {
    "gamma_d_GB":    77.7558,
    "gamma_p_GB":   299.9121,
    "tau_effca_GB": 278.3177658387,
}

DEES_SPECIAL = "/project/ctb-emuller/dhuruva/DEES_cell_packages/x86_64/special"
DEES_LIB     = "/project/ctb-emuller/dhuruva/DEES_cell_packages/x86_64"

SBATCH_TMPL = """\
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --account={account}
#SBATCH --job-name={job_name}
#SBATCH --output={log_file}

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH JPKGDIR PYTHONPATH PATH HOC_LIBRARY_PATH NEURODAMUS_DIR
hash -r

module --force purge
module load StdEnv/2023 scipy-stack/2025a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.4 \\
    cmake/3.31.0 cuda/12.9 mpi4py/4.0.3 pytest/8.2.2 rust/1.91.0 boost/1.85.0
module load python/3.11.5

# Python packages from the newer build (uses libsonata for edges, avoids broken hdf5-mpi h5py)
export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc
# DEES special was compiled against build-20.03.2025-py311/libnrniv.so — put it first.
# Module-set LD_LIBRARY_PATH is preserved at the end (provides hdf5-mpi lib dir etc.)
DEES_NRNIV=/project/def-emuller/opt/jupyterhub-pkgs/build-20.03.2025-py311/lib
export LD_LIBRARY_PATH={dees_lib}:$DEES_NRNIV:$LD_LIBRARY_PATH
# SonataReport is not compiled into DEES special — load it from a minimal lib
# (standard libnrnmech.so conflicts because it re-defines CaDynamics_DC0 already in DEES special)
export NRNMECH_LIB_PATH=/project/ctb-emuller/dhuruva/sonata_mods/x86_64/libnrnmech.so

BASE_DIR=`pwd`
SIM_CFG=$BASE_DIR/{sim_cfg_name}

echo "=== libnrniv.so resolution ==="
ldd {dees_special} | grep nrniv
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "Starting ion-channels neurodamus run (DEES special) ..."
echo "Config: $SIM_CFG"

srun {dees_special} -mpi \\
    -python $JPKGDIR/bin/neurodamus_init.py \\
    --configFile=$SIM_CFG --lb-mode=WholeCell --verbose
"""


def find_workdirs(results_dir, freq):
    pattern = os.path.join(
        results_dir, "fitting", "*", "seed*", "*_STDP",
        "simulations", "*-*", f"{freq}_*ms",
    )
    dirs = sorted(glob.glob(pattern))
    return [d for d in dirs
            if os.path.isfile(os.path.join(d, "prefire_simulation_config.json"))
            and os.path.isfile(os.path.join(d, "prefire_prespikes.h5"))]


def make_ic_config(workdir, src_cfg):
    """Return a modified copy of src_cfg for the ion-channel run."""
    cfg = copy.deepcopy(src_cfg)
    cfg["network"] = IC_CIRCUIT_CONFIG
    # Use absolute path so neurodamus doesn't double-resolve relative to SimulationConfigDir
    cfg["output"]["output_dir"] = os.path.join(os.path.abspath(workdir), f"out_{TAG}")
    # DEES special has no CoreNEURON-compiled mechanisms — run under plain NEURON
    cfg.pop("target_simulator", None)
    mech = cfg["conditions"]["mechanisms"]["GluSynapse"]
    mech["gamma_d_GB"] = DE_FIT_PARAMS["gamma_d_GB"]
    mech["gamma_p_GB"] = DE_FIT_PARAMS["gamma_p_GB"]
    mech["tau_effca_GB"] = DE_FIT_PARAMS["tau_effca_GB"]
    return cfg


def write_ic_files(workdir):
    """Write prefire_simulation_config_{TAG}.json and neurodamus_sbatch_{TAG}.sh."""
    workdir = os.path.abspath(workdir)
    src_path = os.path.join(workdir, "prefire_simulation_config.json")
    with open(src_path) as f:
        src_cfg = json.load(f)

    cfg_ic = make_ic_config(workdir, src_cfg)
    cfg_path = os.path.join(workdir, f"prefire_simulation_config_{TAG}.json")
    with open(cfg_path, "w") as f:
        json.dump(cfg_ic, f, indent=4)

    pair  = os.path.basename(os.path.dirname(workdir))
    proto = os.path.basename(workdir)
    job_name = f"nd{TAG}_{proto}_{pair}"
    log_file = os.path.join(workdir, f"{job_name}.log")
    sbatch = SBATCH_TMPL.format(
        account=ACCOUNT,
        job_name=job_name,
        log_file=log_file,
        dees_special=DEES_SPECIAL,
        dees_lib=DEES_LIB,
        sim_cfg_name=os.path.basename(cfg_path),
    )
    sbatch_path = os.path.join(workdir, f"neurodamus_sbatch_{TAG}.sh")
    with open(sbatch_path, "w") as f:
        f.write(sbatch)
    os.chmod(sbatch_path, 0o755)

    return cfg_path, sbatch_path


A_DE_FIT = dict(a00=1.003498, a01=2.902478, a10=1.644558, a11=2.764812,
                a20=1.003498, a21=2.902478, a30=1.644558, a31=2.764812)
CACHE_PKL = os.path.join(PLASTYFIRE_ROOT, "cpre_cpost_cache", "ion_channels_tau278.pkl")
EDGES_H5  = os.path.join(PLASTYFIRE_ROOT, "data", "dhuruva_modified_edges.h5")
EDGE_POP  = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"


def check_thresholds():
    """Verify edges.h5 theta matches the DE-fit a-params before launching.

    Neurodamus reads theta from the edge file, so a stale injection produces a
    plausible-looking run against the wrong parameters. Fail loudly instead.
    """
    import pickle
    import h5py
    import numpy as np
    with open(CACHE_PKL, "rb") as f:
        cache = pickle.load(f)
    cp, cq = {}, {}
    for v in cache.values():
        for sid, x in v["c_pre"].items():  cp[int(sid)] = x
        for sid, x in v["c_post"].items(): cq[int(sid)] = x
    ids = sorted(cp)
    with h5py.File(EDGES_H5, "r") as f:
        g = f[f"edges/{EDGE_POP}/0"]
        td = g["theta_d"][ids]
        tp = g["theta_p"][ids]
    exp_d = np.array([A_DE_FIT["a00"]*cp[i] + A_DE_FIT["a01"]*cq[i] for i in ids])
    exp_p = np.array([A_DE_FIT["a10"]*cp[i] + A_DE_FIT["a11"]*cq[i] for i in ids])
    ok = (np.abs(td-exp_d) < 1e-6) & (np.abs(tp-exp_p) < 1e-6)
    logger.info("edges.h5 theta matches DE-fit a-params on %d/%d synapses",
                ok.sum(), len(ids))
    if ok.mean() < 0.999:
        logger.error(
            "edges.h5 theta does NOT match the DE-fit a-params (%d/%d). "
            "Run compute_thresholds_from_cache.py + inject_thresholds.py with "
            "a00=%.6f a01=%.6f a10=%.6f a11=%.6f (apical tied) first, "
            "or pass --skip-threshold-check to override.",
            ok.sum(), len(ids), A_DE_FIT["a00"], A_DE_FIT["a01"],
            A_DE_FIT["a10"], A_DE_FIT["a11"])
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS)
    parser.add_argument("--freq", default="10Hz",
                        help="Frequency prefix (default: 10Hz)")
    parser.add_argument("--skip-existing", action="store_true",
                        help=f"Skip workdirs that already have out_{TAG}/rho.h5")
    parser.add_argument("--skip-threshold-check", action="store_true",
                        help="do not verify edges.h5 theta against the DE-fit a-params")
    parser.add_argument("--dry-run", action="store_true",
                        help="Write config/sbatch files but do not sbatch")
    args = parser.parse_args()

    if not args.skip_threshold_check:
        check_thresholds()

    workdirs = find_workdirs(args.results_dir, args.freq)
    if not workdirs:
        logger.error("No workdirs found for freq=%s in %s", args.freq, args.results_dir)
        sys.exit(1)
    logger.info("Found %d workdirs for %s", len(workdirs), args.freq)

    if args.skip_existing:
        before = len(workdirs)
        workdirs = [d for d in workdirs
                    if not os.path.isfile(os.path.join(d, f"out_{TAG}", "rho.h5"))]
        logger.info("  %d already done → %d remaining", before - len(workdirs), len(workdirs))

    if not workdirs:
        logger.info("All done.")
        return

    submitted = failed = 0
    for workdir in workdirs:
        try:
            cfg_path, sbatch_path = write_ic_files(workdir)
        except Exception as e:
            logger.error("Config write failed for %s: %s", workdir, e)
            failed += 1
            continue

        if args.dry_run:
            print(f"sbatch --chdir={workdir} {sbatch_path}")
            submitted += 1
        else:
            result = subprocess.run(
                ["sbatch", f"--chdir={workdir}", sbatch_path],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                logger.info("%s  (%s)", result.stdout.strip(), os.path.basename(workdir))
                submitted += 1
            else:
                logger.error("FAILED %s: %s", os.path.basename(workdir), result.stderr.strip())
                failed += 1

    action = "Would submit" if args.dry_run else "Submitted"
    logger.info("%s %d jobs, %d failed", action, submitted, failed)
    if args.dry_run:
        logger.info("Re-run without --dry-run to actually submit.")


if __name__ == "__main__":
    main()
