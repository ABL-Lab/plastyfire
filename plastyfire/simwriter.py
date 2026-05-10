"""
Writes files for model optimization (more involved)
and simple ones for generalization (AKA finding thresholds based on the optimized parameters)
last modified: András Ecker 06.2024
"""

import os
import argparse
import h5py
import json
import pickle
import pathlib
import shutil
import hashlib
import warnings
from tqdm import tqdm
import numpy as np
import pandas as pd
import threading
from concurrent.futures import ProcessPoolExecutor
from bluepysnap import Circuit
from conntility.connectivity import ConnectivityMatrix

from plastyfire.config import OptConfig, Config
from plastyfire.simulator import spike_threshold_finder

MIN2MS = 60 * 1000.
OPT_CPU_TIME = 2.  # heuristics: it takes ~2x compute time (w/ CVode w/ reporting w/o fastforward) as biological time
CPU_TIME = 1.  # heuristics: c_pre and c_post for a single connection takes <1 minute to simulate/calculate
FIGS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/figures/plastyfire"


def check_geom_constraint(conn_mat, pre_mtype, post_gid, max_dist):
    """Check if cell has any presynaptic partners within `max_dists`"""
    nrn = conn_mat.vertices
    post_mtype = nrn.loc[nrn["node_ids"] == post_gid, "mtype"].to_numpy()[0]
    coords = nrn.loc[nrn["node_ids"] == post_gid, ["ss_flat_x", "depth", "ss_flat_y"]]
    dists = (nrn.loc[nrn["mtype"].isin(pre_mtype), ["ss_flat_x", "depth", "ss_flat_y"]] - coords.to_numpy()).abs()
    if post_mtype in pre_mtype:
        dists.drop(coords.index, inplace=True)
    idx = dists.loc[(dists["ss_flat_x"] < max_dist[0]) &
                    (dists["depth"] < max_dist[1]) &
                    (dists["ss_flat_y"] < max_dist[2])].index.to_numpy()
    valid_gids = nrn.loc[idx, "node_ids"].to_numpy()
    sub_mat = conn_mat.submatrix(valid_gids, sub_gids_post=[post_gid])
    if sub_mat.size:
        return valid_gids[sub_mat.tocoo().row]
    else:
        return None


def check_electrical_constraint(sim_config, gid, stim_config, save_dir):
    """Check if the cell can fire correctly at every stim. frequency
    (and save params. of current injection that makes it fire)"""
    pklf_name = os.path.join(save_dir, "%i.pkl" % gid)
    if os.path.isfile(pklf_name):  # check if these sims were already run...
        return True
    results = {}
    for freq in stim_config["freq"]:
        simres = spike_threshold_finder(sim_config, gid, stim_config["nspikes"],
                                        freq, stim_config["width"], stim_config["offset"], stim_config["amp_min"],
                                        stim_config["amp_max"], stim_config["amp_lev"], fixhp=True)
        if simres is None:
            return False
        else:
            results[freq] = simres
    # Store results for manual validation and simulation setup
    with open(pklf_name, "wb") as f:
        pickle.dump(results, f, -1)
    return True


def save_spikes(h5f_name, prefix, spike_times, spiking_gids):
    """Save spikes to SONATA format"""
    assert (spiking_gids.shape == spike_times.shape)
    with h5py.File(h5f_name, "w") as h5f:
        grp = h5f.require_group("spikes/%s" % prefix)
        grp.create_dataset("timestamps", data=spike_times)
        grp["timestamps"].attrs["units"] = "ms"
        grp.create_dataset("node_ids", data=spiking_gids, dtype=int)


def get_cpu_time(n_afferents):
    """CPU time heuristics: 5 min setup and stimulus calculation + 10 mins. extra just to make sure +
    gid specific simulation time, based on the number of its afferent gids"""
    cpu_time_sec = (15 + n_afferents * CPU_TIME) * 60
    h, m = np.divmod(cpu_time_sec, 3600)
    m, s = np.divmod(m, 60)
    cpu_time_str = "%.2i:%.2i:%.2i" % (h, m, s)
    qos = "#SBATCH --qos=longjob" if h >= 24 else ""
    return cpu_time_str, qos


def plot_evolution(logbook, fig_name):
    """Saves figure with the evolution of fitting error"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_subplot(1, 1, 1)
    gen = np.array([data["gen"] for data in logbook])
    mean = np.array([data["avg"] for data in logbook])
    std = np.array([data["std"] for data in logbook])
    ax.plot(gen, mean, "k-", linewidth=2, label="pop. mean")
    ax.fill_between(gen, mean - std, mean + std, color="lightgray", label="pop. std")
    ax.plot(gen, np.array([data["min"] for data in logbook]), "r-", linewidth=2, label="pop. min")
    ax.legend(frameon=False)
    ax.set_xlabel("Generation")
    ax.set_xlim([1, gen[-1]])
    ax.set_ylabel("Error")
    fig.savefig(fig_name, dpi=100, bbox_inches="tight")
    plt.close()


def plot_epsp_ratios(res_db, fig_name):
    """Quick and dirty plot of EPSP ratios"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 6.5))
    ax = fig.add_subplot(2, 3, 1)
    epsp_ratios = res_db.loc[res_db["protocol_id"] == "mrk97_01", "epsp_ratio"].to_numpy()
    ax.hist(epsp_ratios, bins=10, range=(0, 2.5), color="lightgray")
    ax.axvline(0.98, color="red", label="in vitro: 0.98")
    ax.axvline(np.mean(epsp_ratios), color="black", label="in silico: %.2f" % np.mean(epsp_ratios))
    ax.legend(frameon=False)
    ax.set_ylabel("Count")
    ax.set_title("(Mrk97) f: 2Hz, dt: +5ms")
    ax2 = fig.add_subplot(2, 3, 2)
    epsp_ratios = res_db.loc[res_db["protocol_id"] == "mrk97_02", "epsp_ratio"].to_numpy()
    ax2.hist(epsp_ratios, bins=10, range=(0, 2.5), color="lightgray")
    ax2.axvline(1.01, color="red", label="in vitro: 1.01")
    ax2.axvline(np.mean(epsp_ratios), color="black", label="in silico: %.2f" % np.mean(epsp_ratios))
    ax2.legend(frameon=False)
    ax2.set_title("(Mrk97) f: 5Hz, dt: +5ms")
    ax3 = fig.add_subplot(2, 3, 3)
    epsp_ratios = res_db.loc[res_db["protocol_id"] == "sjh06_02", "epsp_ratio"].to_numpy()
    ax3.hist(epsp_ratios, bins=10, range=(0, 2.5), color="lightgray")
    ax3.axvline(1.06, color="red", label="in vitro: 1.06")
    ax3.axvline(np.mean(epsp_ratios), color="black", label="in silico: %.2f" % np.mean(epsp_ratios))
    ax3.legend(frameon=False)
    ax3.set_xlabel("EPSP ratio")
    ax3.set_title("(Sjh06) f: 50Hz, dt: +10ms")
    ax4 = fig.add_subplot(2, 3, 4)
    epsp_ratios = res_db.loc[res_db["protocol_id"] == "mrk97_08", "epsp_ratio"].to_numpy()
    ax4.hist(epsp_ratios, bins=10, range=(0, 2.5), color="lightgray")
    ax4.axvline(0.79, color="red", label="in vitro: 0.79")
    ax4.axvline(np.mean(epsp_ratios), color="black", label="in silico: %.2f" % np.mean(epsp_ratios))
    ax4.legend(frameon=False)
    ax4.set_xlabel("EPSP ratio")
    ax4.set_ylabel("Count")
    ax4.set_title("(Mrk97) f: 10Hz, dt: -10ms")
    ax5 = fig.add_subplot(2, 3, 5)
    epsp_ratios = res_db.loc[res_db["protocol_id"] == "mrk97_07", "epsp_ratio"].to_numpy()
    ax5.hist(epsp_ratios, bins=10, range=(0, 2.5), color="lightgray")
    ax5.axvline(1.20, color="red", label="in vitro: 1.20")
    ax5.axvline(np.mean(epsp_ratios), color="black", label="in silico: %.2f" % np.mean(epsp_ratios))
    ax5.legend(frameon=False)
    ax5.set_xlabel("EPSP ratio")
    ax5.set_title("(Mrk97) f: 10Hz, dt: +10ms")
    fig.tight_layout()
    fig.savefig(fig_name, dpi=100, bbox_inches="tight")
    plt.close()


class OptSimWriter(OptConfig):
    """Class to setup single cell simulations for the optimization of model parameters"""
    
    def _process_post_gid(self, args):
        """Worker function to process a single post_gid for parallel execution"""
        post_gid, i, conn_mat, pre_mtype, max_dist, sim_config_path, stim_config, save_dir, seed = args
        pre_gids = check_geom_constraint(conn_mat, pre_mtype, post_gid, max_dist)
        if pre_gids is not None:
            if check_electrical_constraint(sim_config_path, post_gid, stim_config, save_dir):
                np.random.seed(seed + i)
                return (np.random.choice(pre_gids, 1)[0], post_gid)
        return None
    
    def find_pairs(self):
        """Finds connected pairs of gids based on constraints specified in the config"""
        save_dir = os.path.join(self.out_dir, "single_cells")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # Write simple simulation config (for checking electrical constraints)
        sim_config = {"run": {"dt": 0.025, "tstop": self.T, "random_seed": self.seed},
                      "network": self.circuit_config,
                      "node_sets_file": self.node_set,
                      "node_set": self.target,
                      "output": {"output_dir": os.path.join(self.sims_dir, "out")}}
        sim_config_path = os.path.join(save_dir, "simulation_config.json")
        with open(sim_config_path, "w", encoding="utf-8") as f:
            json.dump(sim_config, f, indent=4)

        c = Circuit(self.circuit_config)
        # Get connectivity matrix and flatmap locations (used for distance based filtering) with `conntility`
        load_cfg = {"loading": {"base_target": self.target,
                                "properties": ["mtype", "x", "y", "z",
                                               "ss_flat_x", "ss_flat_y", "depth"]},
                    "filtering": [{"column": "mtype",
                                   "values": np.unique(self.pre_mtype + self.post_mtype).tolist()}]}
        conn_mat = ConnectivityMatrix.from_bluepy(c, load_cfg, connectome=self.edge_pop)
        nrn = conn_mat.vertices
        post_gids = nrn.loc[nrn["mtype"].isin(self.post_mtype), "node_ids"].to_numpy()
        np.random.seed(self.seed)
        np.random.shuffle(post_gids)
        # Find pairs (parallelized)
        pairs = []
        pbar = tqdm(total=self.npairs, desc="Finding pairs")
        
        # Prepare arguments for parallel processing
        args_list = [(post_gid, i, conn_mat, self.pre_mtype, self.max_dist, 
                     sim_config_path, self.config["stimulus"], save_dir, self.seed) 
                    for i, post_gid in enumerate(post_gids)]
        
        # Process in parallel
        with ProcessPoolExecutor(max_workers=min(30, self.npairs)) as pool:
            futures = [pool.submit(self._process_post_gid, args) for args in args_list]
            try:
                for future in futures:
                    if len(pairs) >= self.npairs:
                        break
                    try:
                        result = future.result(timeout=300)  # 5 minute timeout per future
                        if result is not None:
                            pairs.append(result)
                            pbar.update(1)
                            if len(pairs) >= self.npairs:
                                break
                    except TimeoutError:
                        print(f"Warning: Operation timed out, skipping...")
                        continue
            finally:
                # Cancel remaining futures when enough pairs found
                for future in futures:
                    if not future.done():
                        future.cancel()
        
        pbar.close()
        if len(pairs) < self.npairs:
            warnings.warn("Not enough pairs found")
        return pairs

    def write_batch_sript(self, f_name, templ, cpu_time):
        """Writes single cell batch script (simulation.batch, for pairrunner.py / bluecellulab)"""
        workdir = os.path.dirname(f_name)
        tmp = os.path.split(workdir)
        name = "%s_%s" % (tmp[1], os.path.split(tmp[0])[1])
        param_args = "--fastforward=%.1f" % self.fastforward if self.fastforward is not None else ""
        unique_id = abs(hash(workdir)) % 100000
        with open(f_name, "w+", encoding="latin1") as f:
            f.write(templ.format(name=name, cpu_time=cpu_time, qos="#SBATCH --chdir=%s" % workdir, log=name,
                                 env=self.env, run=self.run, param_args=param_args, workdir=workdir,
                                 fastforward=0.0, param_hash_arg="",
                                 individual_id_arg=f" --individual_id={unique_id}",
                                 generation_arg=" --generation=0"))

    def write_neurodamus_sbatch(self, workdir, cpu_time):
        """Writes neurodamus_sbatch.sh to directly run the pair simulation with neurodamus/CORENEURON"""
        account = self.config.get("simulator", {}).get("account", "ctb-emuller")
        param_args = "--fastforward=%.1f" % self.fastforward if self.fastforward is not None else ""
        content = """\
#!/bin/bash
#SBATCH --nodes=1                                   # number of nodes
#SBATCH --ntasks-per-node=1                       # tasks per node: MPI procs = nodes*ntasks-per-node
#SBATCH --mem=1G                                   # memory; default unit is megabytes
#SBATCH --time={cpu_time}                          # time (DD-HH:MM)
#SBATCH --account={account}

unset PATH PYTHONPATH LD_LIBRARY_PATH CMAKE_PREFIX_PATH
unset LIBRARY_PATH
unset JPKGDIR
unset PYTHONPATH
unset PATH
unset HOC_LIBRARY_PATH
unset NEURODAMUS_DIR

hash -r

module --force purge
module load StdEnv/2023 scipy-stack/2025a gcc/12.3 openmpi/4.1.5 hdf5-mpi/1.14.4 cmake/3.31.0 cuda/12.9 mpi4py/4.0.3 pytest/8.2.2 rust/1.91.0 boost/1.85.0
module load python/3.11.5

export JPKGDIR=/project/def-emuller/opt/jupyterhub-pkgs/build-07.04.2026-py311
# Bashrc setup
export PYTHONPATH=/cvmfs/soft.computecanada.ca/easybuild/python/site-packages:/cvmfs/soft.computecanada.ca/custom/python/site-packages
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/
export PATH=$PATH:/opt/software/slurm/bin/
export PATH=$PATH:$JPKGDIR/bin
# For NEURON
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python
# for neurodamus (mapping.py)
export PYTHONPATH=$PYTHONPATH:$JPKGDIR/lib/python3.11/site-packages/neurodamus/core/hoc
export HOC_LIBRARY_PATH=$JPKGDIR/lib/python3.11/site-packages/neurodamus/data/hoc:$JPKGDIR/share/neurodamus_neocortex/hoc

export NEURODAMUS_DIR=$JPKGDIR
NEURODAMUS_INIT_PY=$NEURODAMUS_DIR/bin/neurodamus_init.py
NEURODAMUS_PARAMS="{param_args}"

BASE_DIR=`pwd`
SIM_CFG=$BASE_DIR/prefire_simulation_config.json

echo Starting ...
echo Python: `which python`

python - << 'EOF'
import neuron
print("neuron.__version__", neuron.__version__)
print("neuron.__file__", neuron.__file__)
EOF

echo JPKGDIR $JPKGDIR
echo PYTHONPATH $PYTHONPATH
echo PATH $PATH
echo NEURODAMUS_DIR $NEURODAMUS_DIR
echo HOC_LIBRARY_PATH $HOC_LIBRARY_PATH
echo LD_LIBRARY_PATH $LD_LIBRARY_PATH
echo LIBRARY_PATH $LIBRARY_PATH
export CORENRN_DEBUG=1
export NEURON_LOG_LEVEL=DEBUG

echo "=== BUILD ==="
srun $NEURODAMUS_DIR/bin/special -mpi --debug \\
    -python $NEURODAMUS_INIT_PY \\
    --configFile=$SIM_CFG $NEURODAMUS_PARAMS --verbose
""".format(cpu_time=cpu_time, account=account, param_args=param_args)
        sbatch_path = os.path.join(workdir, "neurodamus_sbatch.sh")
        with open(sbatch_path, "w", encoding="utf-8") as f:
            f.write(content)
        os.chmod(sbatch_path, 0o755)

    def write_sim_files(self, pairs):
        """Writes pair, frequency, and dt specific `simulation_config.json` used by `bluecellulab`
        and batch scripts to launch single cell sims"""
        # Copy config file to sims dir (just to make sure that one can more or less know what was run...)
        basedir = os.path.split(os.path.split(self.out_dir)[0])[0]
        if not os.path.exists(basedir):
            os.makedirs(basedir)
        shutil.copyfile(self._config_path, os.path.join(basedir, "%s.yaml" % self.label))
        # Generate presynaptic spike times (for C01 and C02)
        n_spikes_before = int(self.C01_duration * MIN2MS / self.C01_T)
        n_spikes_after = int(self.C02_duration * MIN2MS / self.C02_T)
        before_pre_spikes = [self.offset + i * self.C01_T for i in range(n_spikes_before)]
        after_pre_spikes = [self.offset + n_spikes_before * self.C01_T + self.nreps * self.T +
                            i * self.C02_T for i in range(n_spikes_after)]
        before_duration = n_spikes_before * self.C01_T
        pairing_duration = self.nreps * self.T
        t_stop = before_duration + pairing_duration + n_spikes_after * self.C02_T
        prefire_t_stop = self.offset + pairing_duration + 1000.0
        # CPU time heuristics: estimated sim time + 2 hour buffer
        h, m = np.divmod(OPT_CPU_TIME * t_stop / 1000. + 7200, 3600)
        m, s = np.divmod(m, 60)
        cpu_time = "%.2i:%.2i:%.2i" % (h, m, s)
        with open(os.path.join("templates", "simulation.batch.tmpl"), "r") as f:
            templ = f.read()

        all_sims = []
        for freq in self.freq:
            for dt in self.dt:
                for pre_gid, post_gid in pairs:
                    workdir = os.path.join(self.out_dir, "%i-%i" % (pre_gid, post_gid),
                                           "%iHz_%ims" % (int(freq), int(dt)))
                    if not os.path.exists(workdir):
                        os.makedirs(workdir)
                    # Write node set with idx of pre- and postsynaptic neurons
                    jsonf_name = os.path.join(workdir, "node_sets.json")
                    node_sets = {"precell": {"node_id": [int(pre_gid)], "population": self.node_pop},
                                 "postcell": {"node_id": [int(post_gid)], "population": self.node_pop},
                                 "paircells": {"node_id": [int(pre_gid), int(post_gid)], "population": self.node_pop}}
                    with open(jsonf_name, "w", encoding="utf-8") as f:
                        json.dump(node_sets, f, indent=4)
                    try:  # load pulse amplitude and compute spike delays at the given frequency (independent of dt...)
                        with open(os.path.join(self.out_dir, "single_cells", "%i.pkl" % post_gid), "rb") as f:
                            simres = pickle.load(f)
                        amplitude = simres[freq]["amp"]
                        spike_delay = simres[freq]["t_spikes"] - simres[freq]["t_stimuli"]
                    except IOError:  # fallback to default current pulse
                        warnings.warn("Cannot read stimulus amplitude from cache, using 1 nA")
                        spike_delay = self.nspikes * [self.width / 2.]
                        amplitude = 1.0  # nA
                    # Generate postsynaptic stimulus (only one period, as the stimulus will be periodic)
                    isi = 1000.0 / freq  # Inter Spike Interval (ms)
                    post_spikes = np.array([self.offset + before_duration + i * isi for i in range(self.nspikes)])
                    inputs = {"pulse%i" % i: {"input_type": "current_clamp", "module": "pulse",
                                              "node_set": "postcell",
                                              "delay": post_spike, "duration": pairing_duration, "amp_start": amplitude,
                                              "width": self.width, "frequency": 1000. / self.T}
                              for i, post_spike in enumerate(post_spikes)}
                    
                    prefire_post_spikes = np.array([self.offset + i * isi for i in range(self.nspikes)])
                    prefire_inputs = {"pulse%i" % i: {"input_type": "current_clamp", "module": "pulse",
                                              "node_set": "postcell",
                                              "delay": post_spike, "duration": pairing_duration, "amp_start": amplitude,
                                              "width": self.width, "frequency": 1000. / self.T}
                              for i, post_spike in enumerate(prefire_post_spikes)}

                    # Generate (full) presynaptic spike train used as spike replay stimulus
                    pre_spikes = [self.offset + before_duration + i * isi - dt + spike_delay[i] + j * self.T
                                  for j in range(self.nreps) for i in range(self.nspikes)]
                    pre_spikes = np.array(before_pre_spikes + pre_spikes + after_pre_spikes)
                    h5f_name = os.path.join(workdir, "prespikes.h5")
                    save_spikes(h5f_name, self.node_pop, pre_spikes, pre_gid * np.ones(len(pre_spikes), dtype=int))
                    
                    prefire_pre_spikes = [self.offset + i * isi - dt + spike_delay[i] + j * self.T
                                  for j in range(self.nreps) for i in range(self.nspikes)]
                    prefire_pre_spikes = np.array(prefire_pre_spikes)
                    h5f_prefire_name = os.path.join(workdir, "prefire_prespikes.h5")
                    save_spikes(h5f_prefire_name, self.node_pop, prefire_pre_spikes, pre_gid * np.ones(len(prefire_pre_spikes), dtype=int))
                    
                    # Add spike replay to inputs for neurodamus (needed for correct synapse creation).
                    # node_set is the TARGET (postsynaptic) cell; source is the population for replay matching.
                    inputs["prespikes"] = {"input_type": "spikes", "module": "synapse_replay", "node_set": "postcell",
                                           "delay": 0., "duration": t_stop, "source": self.node_pop,
                                           "spike_file": h5f_name}
                    glusynapse_conditions = {
                        "cao_CR": 2.0,
                        "tau_effca_GB": 278.3177658387,
                        "gamma_d_GB": 101.5387594661,
                        "gamma_p_GB": 216.1841700668,
                        "init_depleted": True,
                        "minis_single_vesicle": False
                    }
                    no_sk_e2_modification = {
                        "name": "no_SK_E2",
                        "node_set": "postcell",
                        "type": "ConfigureAllSections",
                        "section_configure": "%s.gSK_E2bar_SK_E2 = 0"
                    }
                    connection_overrides = [
                        {"name": "plasticity", "source": "paircells", "target": "paircells",
                         "modoverride": "GluSynapse", "weight": 1.0},
                        {"name": "no_vpm_proj", "source": "proj_Thalamocortical_VPM_Source",
                         "target": "hex_O1", "weight": 0.0},
                        {"name": "no_pom_proj", "source": "proj_Thalamocortical_POM_Source",
                         "target": "hex_O1", "weight": 0.0}
                    ]
                    # Write simulation config
                    sim_config = {"run": {"dt": 0.025, "tstop": t_stop, "random_seed": self.seed},
                                  "network": self.circuit_config,
                                  "node_sets_file": jsonf_name,
                                  "node_set": "postcell",
                                  "output": {"output_dir": os.path.join(workdir, "out")},
                                  "inputs": inputs,
                                  "conditions": {
                                        "extracellular_calcium": 2.0,
                                        "v_init": -80.0,
                                        "spike_location": "AIS",
                                        "mechanisms": {"GluSynapse": glusynapse_conditions},
                                        "modifications": [no_sk_e2_modification]
                                    },
                                  "reports": {
                                        "soma": {"cells": "postcell", "type": "compartment",
                                                 "variable_name": "v", "unit": "mV", "dt": 0.1,
                                                 "start_time": 0.0, "end_time": t_stop},
                                        "rho": {"cells": "postcell", "type": "synapse",
                                                 "variable_name": "GluSynapse.rho_GB","sections": "all", "unit": "nd", "dt": 0.1,
                                                 "start_time": 0.0, "end_time": t_stop}
                                    },
                                  "target_simulator": "CORENEURON",
                                  "connection_overrides": connection_overrides}
                    with open(os.path.join(workdir, "simulation_config.json"), "w", encoding="utf-8") as f:
                        json.dump(sim_config, f, indent=4)

                    prefire_inputs["prespikes"] = {"input_type": "spikes", "module": "synapse_replay",
                                                   "node_set": "postcell", "delay": 0.,
                                                   "duration": prefire_t_stop, "source": self.node_pop,
                                                   "spike_file": h5f_prefire_name}
                    prefire_sim_config = {"run": {"dt": 0.025, "tstop": prefire_t_stop, "random_seed": self.seed},
                                  "network": self.circuit_config,
                                  "node_sets_file": jsonf_name,
                                  "node_set": "postcell",
                                  "output": {"output_dir": os.path.join(workdir, "out")},
                                  "inputs": prefire_inputs,
                                  "conditions": {
                                        "extracellular_calcium": 2.0,
                                        "v_init": -80.0,
                                        "spike_location": "AIS",
                                        "mechanisms": {"GluSynapse": glusynapse_conditions},
                                        "modifications": [no_sk_e2_modification]
                                    },
                                  "reports": {
                                        "soma": {"cells": "postcell", "type": "compartment",
                                                 "variable_name": "v", "unit": "mV", "dt": 0.1,
                                                 "start_time": 0.0, "end_time": prefire_t_stop},
                                        "rho": {"cells": "postcell", "type": "synapse",
                                                 "variable_name": "GluSynapse.rho_GB","sections": "all", "unit": "nd", "dt": 0.1,
                                                 "start_time": 0.0, "end_time": t_stop}
                                    },
                                  "target_simulator": "CORENEURON",
                                  "connection_overrides": connection_overrides}
                    with open(os.path.join(workdir, "prefire_simulation_config.json"), "w", encoding="utf-8") as f:
                        json.dump(prefire_sim_config, f, indent=4)

                
                    # Write launch scripts
                    f_name = os.path.join(workdir, "simulation.batch")
                    self.write_batch_sript(f_name, templ, cpu_time)
                    self.write_neurodamus_sbatch(workdir, "00:30:00")
                    all_sims.append((pre_gid, post_gid, freq, dt, f_name))
        sim_idx = pd.DataFrame(all_sims, columns=["pregid", "postgid", "frequency", "dt", "path"])
        sim_idx.to_csv(os.path.join(basedir, "index_%s.csv" % self.label), index=False)

    def read_opt_params(self):
        """Loads latest `bluepyopt` checkpoint file, plots results and return optimal parameter set"""
        basedir = os.path.split(os.path.split(self.out_dir)[0])[0]
        with open(os.path.join(basedir, "checkpoint.pkl"), "rb") as f:
            tmp = pickle.load(f)
        gen = tmp["logbook"][-1]["gen"]
        plot_evolution(tmp["logbook"], os.path.join(FIGS_DIR, "fitting.png"))
        errors = [np.linalg.norm(np.array(ind.fitness.values)) for ind in tmp["halloffame"]]
        fit_params = {param_name: param_value for param_name, param_value
                      in zip(tmp["param_names"], tmp["halloffame"].items[np.argmin(errors)])}
        cachekey = hashlib.md5(str(list(fit_params.values())).encode()).hexdigest()
        with open(os.path.join(basedir, ".cache", "%s.pkl" % cachekey), "rb") as f:
            tmp = pickle.load(f)
        print("In silico EPSP ratios:", tmp["outcome"])
        plot_epsp_ratios(tmp["resdb"], os.path.join(FIGS_DIR, "EPSP_ratios_gen%i.png" % gen))
        return fit_params


class SimWriter(Config):
    """Class to setup single cell simulations for finding C_pre and C_post for all synapses"""
    def write_batch_sript(self, f_name, templ, gid, cpu_time, qos):
        """Writes single cell batch script"""
        with open(f_name, "w+", encoding="latin1") as f:
            f.write(templ.format(name="plast_%i" % gid, cpu_time=cpu_time, qos=qos, log=gid,
                                 env=self.env, run=self.run, args="%s %s" % (self._config_path, gid)))

    def write_sim_files(self):
        """Writes simple `simulation_config.json` used by `bluecellulab` and batch scripts for single cell sims"""
        # create and write simple simulation config
        pathlib.Path(self.sims_dir).mkdir(exist_ok=True)
        sim_config = {"run": {"dt": 0.025, "tstop": 3000.0, "random_seed": self.seed},
                      "network": self.circuit_config,
                      "node_sets_file": self.node_set,
                      "node_set": self.target,
                      "output": {"output_dir": "out"},
                      "connection_overrides": [{"name": "plasticity", "source": self.target, "target": self.target,
                                               "modoverride": "GluSynapse", "weight": 1.0}]}
        with open(self.sim_config, "w", encoding="utf-8") as f:
            json.dump(sim_config, f, indent=4)

        # create folders for batch scripts and output csv files
        sbatch_dir = os.path.join(self.sims_dir, "sbatch")
        pathlib.Path(sbatch_dir).mkdir(exist_ok=True)
        pathlib.Path(os.path.join(self.sims_dir, "out")).mkdir(exist_ok=True)
        # get all EXC gids and write sbatch scripts for all of them
        c = Circuit(self.circuit_config)
        df = c.nodes[self.node_pop].get(self.target, "synapse_class")
        gids = df.loc[df == "EXC"].index.to_numpy()  # just to make sure
        with open(os.path.join("templates", "single_cell.batch.tmpl"), "r") as f:
            templ = f.read()
        f_names = []
        for gid in tqdm(gids, desc="Writing batch scripts for every EXC gid", miniters=len(gids)/100):
            f_name = os.path.join(sbatch_dir, "sim_%i.batch" % gid)
            f_names.append(f_name)
            n_afferents = len(np.intersect1d(c.edges[self.edge_pop].afferent_nodes(gid), gids))
            cpu_time, qos = get_cpu_time(n_afferents)
            self.write_batch_sript(f_name, templ, gid, cpu_time, qos)
        # write master launch scripts in batches of 5k
        idx = np.arange(0, len(f_names), 5000)
        idx = np.append(idx, len(f_names))
        for i, (start, end) in enumerate(zip(idx[:-1], idx[1:])):
            with open(os.path.join(sbatch_dir, "launch_batch%i.sh" % i), "w") as f:
                for f_name in f_names[start:end]:
                    f.write("sbatch %s\n" % f_name)

    def relaunch_failed_jobs(self, error, verbose=False):
        """Checks output files and if they aren't presents checks logs for specific `error`
        and creates master launch script to relaunch all failed jobs"""
        c = Circuit(self.circuit_config)
        df = c.nodes[self.node_pop].get(self.target, "synapse_class")
        gids = df.loc[df == "EXC"].index.to_numpy()  # just to make sure
        f_names = []
        for gid in tqdm(gids, desc="Checking log files", miniters=len(gids)/100):
            if not os.path.isfile(os.path.join(self.sims_dir, "out", "%i.csv" % gid)):
                f_name = os.path.join(self.sims_dir, "sbatch", "sim_%i.log" % gid)
                if os.path.isfile(f_name):
                    if verbose:
                        print(f_name)
                    with open(f_name, "r") as f:
                        if error in f.readlines()[-1]:
                            f_names.append(os.path.join(self.sims_dir, "sbatch", "sim_%i.batch" % gid))
        if len(f_names):
            with open(os.path.join(self.sims_dir, "sbatch", "relaunch_failed.sh"), "w") as f:
                for f_name in f_names:
                    f.write("sbatch %s\n" % f_name)
            if verbose:
                print("Generated relaunch_failed.sh master launch script with %i jobs" % len(f_names))

    def check_failed_thresholds(self):
        """Check log files and returns statistics about failed threshold calibrations (for L6 PCs)"""
        c = Circuit(self.circuit_config)
        df = c.nodes[self.node_pop].get(self.target, ["synapse_class", "layer", "mtype"])
        df = df.loc[(df["synapse_class"] == "EXC") & (df["layer"] == "6")]
        gids, mtypes = df.index.to_numpy(), df["mtype"].to_numpy()
        not_defined_ths = {}
        for gid, mtype in tqdm(zip(gids, mtypes), total=len(gids),
                               desc="Checking log files", miniters=len(gids) / 100):
            f_name = os.path.join(self.sims_dir, "sbatch", "sim_%i.log" % gid)
            with open(f_name, "r") as f:
                if "setting negative thresholds" in f.readlines()[-3]:
                    if mtype in not_defined_ths:
                        not_defined_ths[mtype] += 1
                    else:
                        not_defined_ths[mtype] = 1
        unique_mtypes, counts = np.unique(mtypes, return_counts=True)
        for mtype, count in not_defined_ths.items():
            n = counts[unique_mtypes == mtype][0]
            print("For %s: %i gids (%.2f%% of total) couldn't be calibrated" % (mtype, count, (count/n) * 100))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write simulation files for plastyfire optimization")
    parser.add_argument("--npairs", type=int, default=None,
                        help="Override the number of pairs to find (default: read from config)")
    args = parser.parse_args()

    writer = OptSimWriter("../configs/L5TTPC_L5TTPC_STDP.yaml")
    if args.npairs is not None:
        writer.config["npairs"] = args.npairs
    pairs = writer.find_pairs()
    writer.write_sim_files(pairs)
    # writer = OptSimWriter("../configs/L23PC_L5TTPC_STDP.yaml")
    # pairs = writer.find_pairs()
    # writer.write_sim_files(pairs)

    # writer = OptSimWriter("../configs/L5TTPC_L5TTPC.yaml")
    # fit_params = writer.read_opt_params()
    # TODO: rewrite a couple of files with opt. params and run them with `pairrunner.py`

    # writer = SimWriter("../configs/Zenodo_O1.yaml")
    # writer.write_sim_files()
    # writer.relaunch_failed_jobs("slurmstepd:", True)
    # writer.check_failed_thresholds()

