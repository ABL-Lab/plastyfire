#!/usr/bin/env python3
"""Neurodamus validation of a fitted param set (fit.py --save-json) on any emodel.

Per workdir writes prefire_simulation_config_<tag>.json (network -> the given
circuit config, output -> out_<tag>/, GluSynapse gammas/tau from the JSON), then
submits ONE slurm array over all workdirs and prints ARRAY_JOB_ID=<id>.
Compare afterwards with compare_nd_bcl_prefire.py --tag <tag>.

theta_d / theta_p come from dhuruva_modified_edges.h5 in neurodamus (no runtime
a-param path), so they must be injected from the SAME cache first:
    compute_thresholds_from_cache.py --cache <cache> --a00 .. --a31 ..
    inject_thresholds.py --force
check_thresholds() refuses to submit otherwise.

Runs with neurodamus 4.1 (build-22.07.2026-py311) and a special compiled from the
DEES mods (other_mods/modified_mechanisms, identical to what bluecellulab loads)
plus neurodamus' ALU/SonataReport/current-source mods:
    DEES_cell_packages/nd_build_22.07.2026/x86_64/special

    python submit_neurodamus_fit.py --tag delta-cooker \\
        --params-json fit_results/delta-cooker.json \\
        --cache cpre_cpost_cache/sabrina_n120_delta.pkl \\
        --circuit-config data/dhuruva_delta_circuit_config.json
"""
import argparse, copy, glob, json, os, pickle, subprocess, sys

import numpy as np

PLASTYFIRE_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SIMS = os.path.join(PLASTYFIRE_ROOT, "refitting_results/fitting/n120/seed20262009/"
                                             "Sabrina_L5TTPC_L5TTPC_STDP/simulations")
EDGES_H5 = os.path.join(PLASTYFIRE_ROOT, "data", "dhuruva_modified_edges.h5")
EDGE_POP = "S1nonbarrel_neurons__S1nonbarrel_neurons__chemical"
ACCOUNT = "rrg-emuller"
ND_BUILD = "/project/rrg-emuller/dhuruva/DEES_cell_packages/nd_build_22.07.2026/x86_64"
A_KEYS = ("a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31")

ARRAY_TMPL = """\
#!/bin/bash
#SBATCH --job-name=nd_{tag}
#SBATCH --account={account}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=1-{n}
#SBATCH --output={logdir}/nd_{tag}_%a.log

WD=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {listfile})
cd $WD

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH LIBRARY_PATH JPKGDIR HOC_LIBRARY_PATH NEURODAMUS_DIR
hash -r
module --force purge
module load StdEnv/2023 scipy-stack/2025a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.4 \\
    cmake/3.31.0 cuda/12.9 mpi4py/4.0.3 pytest/8.2.2 rust/1.91.0 boost/1.85.0
module load python/3.11.5

export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-22.07.2026-py311
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/:$JPKGDIR/lib/python
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export PATH=$PATH:/opt/software/slurm/bin/:$JPKGDIR/bin
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc:$JPKGDIR/share/nrn/lib/hoc
export LD_LIBRARY_PATH={nd_build}:$JPKGDIR/lib:$LD_LIBRARY_PATH
# NEURON's compiled-in prefix is the builder's home (unreadable) -> point it at JPKGDIR
export NRNHOME=$JPKGDIR

echo "workdir: $WD"
srun {nd_build}/special -mpi -python $JPKGDIR/bin/neurodamus_init.py \\
    --configFile=$WD/prefire_simulation_config_{tag}.json --lb-mode=WholeCell --verbose
"""


def find_workdirs(sims, freq):
    return [d for d in sorted(glob.glob(os.path.join(sims, "*-*", f"{freq}_*ms")))
            if os.path.isfile(os.path.join(d, "prefire_simulation_config.json"))
            and os.path.isfile(os.path.join(d, "prefire_prespikes.h5"))]


def load_params(path):
    p = json.load(open(path))
    p = p.get("params", p)
    missing = [k for k in A_KEYS + ("gamma_d_GB_GluSynapse", "gamma_p_GB_GluSynapse") if k not in p]
    if missing:
        sys.exit(f"{path}: missing {missing}")
    p.setdefault("tau_effca_GB_GluSynapse", 278.3177658387)
    return p


def check_thresholds(params, cache_pkl):
    """edges.h5 theta must equal a-params x the cache's c_pre/c_post."""
    import h5py
    cache = pickle.load(open(cache_pkl, "rb"))
    cp, cq = {}, {}
    for v in cache.values():
        cp.update({int(s): x for s, x in v["c_pre"].items()})
        cq.update({int(s): x for s, x in v["c_post"].items()})
    ids = sorted(cp)
    with h5py.File(EDGES_H5, "r") as f:
        g = f[f"edges/{EDGE_POP}/0"]
        td, tp = g["theta_d"][ids], g["theta_p"][ids]
        sec = g["afferent_section_type"][ids]
    apical = sec == 3
    a = params
    c_pre, c_post = np.array([cp[i] for i in ids]), np.array([cq[i] for i in ids])
    exp_d = np.where(apical, a["a20"] * c_pre + a["a21"] * c_post, a["a00"] * c_pre + a["a01"] * c_post)
    exp_p = np.where(apical, a["a30"] * c_pre + a["a31"] * c_post, a["a10"] * c_pre + a["a11"] * c_post)
    ok = (np.abs(td - exp_d) < 1e-4) & (np.abs(tp - exp_p) < 1e-4)
    print(f"edges.h5 theta matches params on {ok.sum()}/{len(ids)} synapses")
    if ok.mean() < 0.999:
        sys.exit("edges.h5 thetas do not match -- run compute_thresholds_from_cache.py + "
                 "inject_thresholds.py with these a-params and this cache first")


def write_config(workdir, tag, params, circuit):
    cfg = json.load(open(os.path.join(workdir, "prefire_simulation_config.json")))
    cfg = copy.deepcopy(cfg)
    cfg["network"] = circuit
    cfg["output"]["output_dir"] = os.path.join(workdir, f"out_{tag}")
    cfg.pop("target_simulator", None)          # DEES special has no CoreNEURON mechs
    for m in cfg["conditions"].get("modifications", []):
        if m.get("type") == "configure_all_sections":   # libsonata 22.07 wants CamelCase
            m["type"] = "ConfigureAllSections"
    mech = cfg["conditions"]["mechanisms"]["GluSynapse"]
    mech["gamma_d_GB"] = params["gamma_d_GB_GluSynapse"]
    mech["gamma_p_GB"] = params["gamma_p_GB_GluSynapse"]
    mech["tau_effca_GB"] = params["tau_effca_GB_GluSynapse"]
    with open(os.path.join(workdir, f"prefire_simulation_config_{tag}.json"), "w") as f:
        json.dump(cfg, f, indent=4)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--params-json", required=True)
    ap.add_argument("--cache", required=True, help="cpre/cpost cache the thetas were computed from")
    ap.add_argument("--circuit-config", required=True, help="picks the emodel")
    ap.add_argument("--results-dir", default=DEFAULT_SIMS, help=".../simulations dir")
    ap.add_argument("--freq", default="10Hz")
    ap.add_argument("--skip-existing", action="store_true", help="skip workdirs with out_<tag>/rho.h5")
    ap.add_argument("--skip-threshold-check", action="store_true")
    ap.add_argument("--limit", type=int, default=None, help="first N workdirs only (pilot)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    params = load_params(a.params_json)
    circuit = os.path.abspath(a.circuit_config)
    if not a.skip_threshold_check:
        check_thresholds(params, a.cache)

    wds = [os.path.abspath(w) for w in find_workdirs(a.results_dir, a.freq)]
    if a.skip_existing:
        wds = [w for w in wds if not os.path.isfile(os.path.join(w, f"out_{a.tag}", "rho.h5"))]
    if a.limit:
        wds = wds[:a.limit]
    if not wds:
        print("nothing to run"); print("ARRAY_JOB_ID="); return
    for w in wds:
        write_config(w, a.tag, params, circuit)

    logdir = os.path.join(PLASTYFIRE_ROOT, "logs", f"nd_{a.tag}")
    os.makedirs(logdir, exist_ok=True)
    listfile = os.path.join(logdir, "workdirs.txt")
    open(listfile, "w").write("\n".join(wds) + "\n")
    script = os.path.join(logdir, "nd_array.sh")
    open(script, "w").write(ARRAY_TMPL.format(tag=a.tag, account=ACCOUNT, n=len(wds), logdir=logdir,
                                              listfile=listfile, nd_build=ND_BUILD))
    print(f"{len(wds)} workdirs -> {script}")
    if a.dry_run:
        print("ARRAY_JOB_ID="); return
    # when called from inside a job (the pipeline), the inherited SLURM_* vars
    # (e.g. SLURM_TRES_PER_TASK=cpu=60) make the array's srun abort immediately
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    r = subprocess.run(["sbatch", "--parsable", script], capture_output=True, text=True, env=env)
    if r.returncode:
        sys.exit(f"sbatch failed: {r.stderr.strip()}")
    print(f"ARRAY_JOB_ID={r.stdout.strip().split(';')[0]}")


if __name__ == "__main__":
    main()
