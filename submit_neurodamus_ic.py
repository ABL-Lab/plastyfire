#!/usr/bin/env python3
"""
Generate prefire_simulation_config_ic.json and neurodamus_sbatch_ic.sh for every
workdir, pointing neurodamus at the new ion-channel circuit config and gen-4 best
params.  Writes to out_ic/ instead of out/ so baseline results are untouched.

Changes vs the baseline prefire_simulation_config.json:
  - "network"  → dhuruva_modified_ion_channels_circuit_config.json
  - output_dir → {workdir}/out_ic
  - GluSynapse.gamma_d_GB → 101.5  (gen-4 best)
  - GluSynapse.gamma_p_GB → 199.773931  (gen-4 best)
  - GluSynapse.tau_effca_GB stays 278.3177658387

theta_d / theta_p are read from dhuruva_modified_edges.h5 which already holds
the gen-4 a-param values injected by inject_thresholds.py.

Usage:
    python submit_neurodamus_ic.py --dry-run
    python submit_neurodamus_ic.py
    python submit_neurodamus_ic.py --skip-existing
    python submit_neurodamus_ic.py --freq 10 --dry-run
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

# Gen-4 ion-channel optimizer best
GEN4_PARAMS = {
    "tau_effca_GB": 278.3177658387,
    "gamma_d_GB":   101.5,
    "gamma_p_GB":   199.773931,
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
SIM_CFG=$BASE_DIR/prefire_simulation_config_ic.json

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
    cfg["output"]["output_dir"] = os.path.join(os.path.abspath(workdir), "out_ic")
    # DEES special has no CoreNEURON-compiled mechanisms — run under plain NEURON
    cfg.pop("target_simulator", None)
    mech = cfg["conditions"]["mechanisms"]["GluSynapse"]
    mech["gamma_d_GB"] = GEN4_PARAMS["gamma_d_GB"]
    mech["gamma_p_GB"] = GEN4_PARAMS["gamma_p_GB"]
    mech["tau_effca_GB"] = GEN4_PARAMS["tau_effca_GB"]
    return cfg


def write_ic_files(workdir):
    """Write prefire_simulation_config_ic.json and neurodamus_sbatch_ic.sh."""
    workdir = os.path.abspath(workdir)
    src_path = os.path.join(workdir, "prefire_simulation_config.json")
    with open(src_path) as f:
        src_cfg = json.load(f)

    cfg_ic = make_ic_config(workdir, src_cfg)
    cfg_path = os.path.join(workdir, "prefire_simulation_config_ic.json")
    with open(cfg_path, "w") as f:
        json.dump(cfg_ic, f, indent=4)

    pair  = os.path.basename(os.path.dirname(workdir))
    proto = os.path.basename(workdir)
    job_name = f"ndic_{proto}_{pair}"
    log_file = os.path.join(workdir, f"{job_name}.log")
    sbatch = SBATCH_TMPL.format(
        account=ACCOUNT,
        job_name=job_name,
        log_file=log_file,
        dees_special=DEES_SPECIAL,
        dees_lib=DEES_LIB,
    )
    sbatch_path = os.path.join(workdir, "neurodamus_sbatch_ic.sh")
    with open(sbatch_path, "w") as f:
        f.write(sbatch)
    os.chmod(sbatch_path, 0o755)

    return cfg_path, sbatch_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS)
    parser.add_argument("--freq", default="10Hz",
                        help="Frequency prefix (default: 10Hz)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip workdirs that already have out_ic/rho.h5")
    parser.add_argument("--dry-run", action="store_true",
                        help="Write config/sbatch files but do not sbatch")
    args = parser.parse_args()

    workdirs = find_workdirs(args.results_dir, args.freq)
    if not workdirs:
        logger.error("No workdirs found for freq=%s in %s", args.freq, args.results_dir)
        sys.exit(1)
    logger.info("Found %d workdirs for %s", len(workdirs), args.freq)

    if args.skip_existing:
        before = len(workdirs)
        workdirs = [d for d in workdirs
                    if not os.path.isfile(os.path.join(d, "out_ic", "rho.h5"))]
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
