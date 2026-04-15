"""
Single cell simulations in bluecellulab (based on Giuseppe Chindemi's `BGLibPy` code)
last modified: András Ecker, 06.2024
"""

import importlib
import os
import re
import logging
import multiprocessing
import numpy as np
import pandas as pd
from functools import lru_cache
from libsonata import SpikeReader
from bluepysnap import Simulation
import bluecellulab
from conntility.io.synapse_report import get_presyn_mapping

bluecellulab.set_verbose(2)
bluecellulab.neuron.h.cvode.atolscale("v", .1)
bluecellulab.neuron.load_mechanisms("/project/ctb-emuller/dhuruva/DEES_cell_packages/")
SPIKE_THRESHOLD = -30  # mV
EXTRA_RECIPE_PATH = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/recipe_andras.csv"
with_cache = lru_cache(128)  # set cache for spiking thresholds
logger = logging.getLogger(__name__)  # configure logger
EPG_MODULES = {
    "epg": "plastyfire.epg",
    "epg_dhuruva": "plastyfire.epg_dhuruva",
    "epg_dhuruva_full": "plastyfire.epg_dhuruva_full",
    "epg_dhuruva_custom_ratio": "plastyfire.epg_dhuruva_custom_ratio",
}
# DEBUG flag will be set by modelfitter.py when --debug is used
DEBUG = False  # Default to False, will be overridden by modelfitter
SYNPROPS = ["Cpre", "Cpost", "loc", "Use0_TM", "Dep_TM", "Fac_TM", "Nrrp_TM", "gmax0_AMPA", "gmax_NMDA",
            "volume_CR", "synapseID", "theta_d_GB", "theta_p_GB"]
MOD_PROPS = ["gamma_d_GB", "gamma_p_GB", "g_IP3R_CICR", "tau_IP3R_CICR"]  # 'tau_exp_GB'
#SYNREC = ["rho_GB", "Use_GB", "gmax_AMPA", "cai_CR", "vsyn", "ica_NMDA", "ica_VDCC", "effcai_GB", "shaft_cai"]
SYNREC = ["rho_GB","cai_CR", "vsyn", "effcai_GB", "shaft_cai","ica_NMDA", "ica_VDCC", "g_NMDA", "cai_NMDA_CR", "cai_VDCC_CR"]
# because of the constant ping-pong between GluSynapse.mod, SONATA parameter names
# and synapse helpers both in `neurodamus` and in bluecellulab.synapses.synapse_types/GluSynapse() mimicking it


def _get_params_generator(epg_variant):
    """Resolve the configured parameter generator implementation."""
    module_name = EPG_MODULES.get(epg_variant)
    if module_name is None:
        raise ValueError(
            f"Unknown epg_variant '{epg_variant}'. Expected one of: {', '.join(sorted(EPG_MODULES))}"
        )
    module = importlib.import_module(module_name)
    return module.ParamsGenerator

# some variable names have to be patched (to match the current state of GluSynapse.mod)
PARAM_MAP = {"Use_d_TM": "Use_d", "Use_p_TM": "Use_p", "Use0_TM": "Use",
             "Dep_TM": "Dep", "Fac_TM": "Fac", "Nrrp_TM": "Nrrp"}


def _get_spikes(t, v, dt_int=0.025):
    """Interprets fix `dt` time (and voltage) from CVode results, and finds spike times (skips 200 ms)"""
    tdense = np.linspace(min(t), max(t), int((max(t)-min(t))/dt_int))
    vdense = np.interp(tdense, t, v)
    return np.array([tdense[i+1] for i in range(int(200 / dt_int), len(vdense) - 1)
                     if vdense[i] < SPIKE_THRESHOLD and vdense[i+1] >= SPIKE_THRESHOLD])


def _runsinglecell_proc(sim_config, post_gid, stimulus, results, node_pop, fixhp):
    """Multiprocessing subprocess for `runsinglecell()`"""
    sim = bluecellulab.CircuitSimulation(sim_config)
    sim.instantiate_gids([(node_pop, post_gid)])
    cell = sim.cells[(node_pop, post_gid)]
    if fixhp:  # hyperpolarization workaround
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
    # add stimuli and run sim
    tstim = bluecellulab.neuron.h.TStim(0.5, sec=cell.soma)
    stim_duration = (stimulus["nspikes"] - 1) * 1000. / stimulus["freq"] + stimulus["width"]
    tstim.train(stimulus["offset"], stim_duration, stimulus["amp"], stimulus["freq"], stimulus["width"])
    cell.persistent.append(tstim)
    sim.run(stimulus["offset"] + stim_duration + 200., cvode=True)
    # get soma voltage and simulation time vector and extract spike times
    t = np.array(sim.get_time())
    v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
    spikes = _get_spikes(t, v)
    # store results
    results["t"] = t
    results["v"] = v
    results["t_spikes"] = spikes
    results["t_stimuli"] = np.array(tstim.tvec)[:-1:4]
    results.update(stimulus)


def runsinglecell(sim_config, post_gid, stimulus, node_pop, fixhp):
    """Runs single cell simulation with given stimulus"""
    manager = multiprocessing.Manager()
    logger.debug("Submitting simulation: post_gid={}, f={}, a={}".format(post_gid, stimulus["freq"], stimulus["amp"]))
    results = manager.dict()
    p = multiprocessing.Process(target=_runsinglecell_proc, args=(sim_config, post_gid, stimulus, results, node_pop, fixhp))
    p.start()
    p.join()
    return dict(results)


@with_cache
def spike_threshold_finder(sim_config, post_gid, nspikes, freq, width, offset, min_amp, max_amp, nlevels,
                           node_pop="S1nonbarrel_neurons", fixhp=False):
    """
    Finds the min amplitude of stimulus current (within the range [`min_amp`, `max_amp`])
    that makes the `post_gid` fire `nspikes` APs (given `freq`, `width`, `offset`)
    using a binary search (parametrized by `nlevels`).
    """
    candidate_amp = np.linspace(min_amp, max_amp, nlevels)
    # find suitable amplitude (binary search, leftmost element)
    L = 0
    R = nlevels
    simres = {}
    while L < R:
        m = int(np.floor((L+R)/2.))
        stim = {"nspikes": nspikes,
                "freq": freq,
                "width": width,
                "offset": offset,
                "amp": candidate_amp[m]}
        simres[m] = runsinglecell(sim_config, post_gid, stim, node_pop, fixhp)
        logger.debug("Number of spikes = %d" % len(simres[m]["t_spikes"]))
        if len(simres[m]["t_spikes"]) < nspikes:
            L = m + 1
        else:
            R = m
    # L is the index of the best amplitude, but we don't know if the match was exact
    t_stim = 1000. / freq * np.array(range(nspikes)) + offset
    t_lim = t_stim + width + 5
    if L == nlevels:
        logger.debug("Max stimulation intensity too weak")
    elif len(simres[L]["t_spikes"]) != nspikes:
        logger.debug("Search grid too coarse")
    elif np.any(simres[L]["t_spikes"] < t_stim) or np.any(simres[L]["t_spikes"] > t_lim):
        logger.debug("Spikes out of order")
    else:
        logger.debug("Correct spike count for stimulation amplitude = %.3f nA" % simres[L]["amp"])
        return simres[L]
    return None


# RANGE params in GluSynapse.mod that cannot be set as HOC globals
# These must be set per-synapse in _set_local_params instead
_RANGE_PARAMS = {"enable_CICR_GluSynapse"}

def _set_global_params(allparams):
    """Sets global parameters of the simulation"""
    logger.debug("Setting global parameters")
    for param_name, param_val in allparams.items():
        if param_name in _RANGE_PARAMS:
            continue  # RANGE params are set per-synapse in _set_local_params
        if re.match(".*_GluSynapse$", param_name):
            try:
                setattr(bluecellulab.neuron.h, param_name, param_val)
            except Exception:
                # HOC global doesn't exist in the loaded .mod (e.g. Vmax_CICR_GluSynapse
                # when simulations were built without the CICR extension) — skip gracefully.
                logger.warning("Skipping unknown HOC global: %s (not defined in loaded .mod)", param_name)
            # logger.debug("\t%s = %f", param_name, getattr(bluecellulab.neuron.h, param_name))


def _set_local_params(synapse, fit_params, extra_params, c_pre=0., c_post=0.):
    """Sets synaptic parameters in bluecellulab"""
    for key, val in extra_params.items():  # update basic synapse parameters
        if key in PARAM_MAP:
            setattr(synapse.hsynapse, PARAM_MAP[key], val)
        else:
            if key  == "loc":
                continue
            setattr(synapse.hsynapse, key, val)
    if fit_params is not None:  # update thresholds
        if all(key in fit_params for key in ["a00", "a01"]) and extra_params["loc"] == "basal":
            # set basal depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a00"] * c_pre + fit_params["a01"] * c_post
        if all(key in fit_params for key in ["a10", "a11"]) and extra_params["loc"] == "basal":
            # set basal potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a10"] * c_pre + fit_params["a11"] * c_post
        if all(key in fit_params for key in ["a20", "a21"]) and extra_params["loc"] == "apical":
            # set apical depression threshold
            synapse.hsynapse.theta_d_GB = fit_params["a20"] * c_pre + fit_params["a21"] * c_post
        if all(key in fit_params for key in ["a30", "a31"]) and extra_params["loc"] == "apical":
            # set apical potentiation threshold
            synapse.hsynapse.theta_p_GB = fit_params["a30"] * c_pre + fit_params["a31"] * c_post
        # Set RANGE params that can't be set as globals
        if "enable_CICR_GluSynapse" in fit_params:
            synapse.hsynapse.enable_CICR = fit_params["enable_CICR_GluSynapse"]
    


def _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop):
    """Compared to `bluepysnap` which has one global synapse ID for all synapses in the circuit,
    `bluecellulab` (just as `neurodamus`) re-indexes synapses for each postsynaptic cell starting at 0
    (which ID is used for seeding the synapses). This helper gets the mapping between the two indexing versions
    It's unfortunate that this is here (and gets called so many times)... but I couldn't find a better way"""
    return get_presyn_mapping(Simulation(sim_config).circuit, edge_pop,
                              pd.MultiIndex.from_tuples([(post_gid, syn_id) for syn_id in syn_idx]))


def _c_pre_finder_process(sim_config, fit_params, syn_extra_params, pre_gid, post_gid,
                          node_pop, edge_pop, fixhp):
    """
    Multiprocessing subprocess for `c_pre_finder()`
    Delivers spike from `pre_gid` and measures the Ca++ transient at the synapses on `post_gid`
    """
    logger.debug("c_pre finder process")
    sim = bluecellulab.CircuitSimulation(sim_config)
    sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, add_minis=False,
                          pre_spike_trains={(node_pop, pre_gid): [1000.]},
                          intersect_pre_gids=[(node_pop, pre_gid)])
    cell = sim.cells[(node_pop, post_gid)]
    if fixhp:  # hyperpolarization workaround
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
    if fit_params is not None:  # setup global parameters
        _set_global_params(fit_params)
    # initialize [Ca^{2+}]_i and vsyn recorders
    recorder, recorder_cai, recorder_shaft_cai, recorder_vsyn, recorder_ica = {}, {}, {}, {}, {}
    recorder_ica_NMDA, recorder_ica_VDCC, syn_idx = {}, {}, []
    recorder_cai_NMDA_CR, recorder_cai_VDCC_CR = {}, {}
    recorder_t = bluecellulab.neuron.h.Vector()
    recorder_t.record(bluecellulab.neuron.h._ref_t, 1.0)

    for syn_id, synapse in cell.synapses.items():
        syn_idx.append(syn_id[1])
        recorder[syn_id[1]] = bluecellulab.neuron.h.Vector()
        recorder[syn_id[1]].record(synapse.hsynapse._ref_effcai_GB, 1.0)
        recorder_cai[syn_id[1]] = bluecellulab.neuron.h.Vector()
        recorder_cai[syn_id[1]].record(synapse.hsynapse._ref_cai_CR, 1.0)

        # recorder_cai_NMDA_CR[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_cai_NMDA_CR[syn_id[1]].record(synapse.hsynapse._ref_cai_NMDA_CR, 1.0)
        # recorder_cai_VDCC_CR[syn_id[1]] = bluecellulab.neuron.h.Vector()    
        # recorder_cai_VDCC_CR[syn_id[1]].record(synapse.hsynapse._ref_cai_VDCC_CR, 1.0)

        # recorder_shaft_cai[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_shaft_cai"):
        #     recorder_shaft_cai[syn_id[1]].record(synapse.hsynapse._ref_shaft_cai)
        # # vsyn: synaptic voltage state variable in GluSynapse (CVode-adaptive time axis)
        # recorder_vsyn[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_vsyn[syn_id[1]].record(synapse.hsynapse._ref_vsyn)
        # # ica: total calcium current density (mA/cm2) at the dendritic segment
        # # This is what CaDynamics.mod reads to compute shaft calcium
        # seg = synapse.hsynapse.get_segment()
        # recorder_ica[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_ica[syn_id[1]].record(seg._ref_ica)
        # # ica_NMDA and ica_VDCC: synapse-level calcium current components
        # recorder_ica_NMDA[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_ica_NMDA"):
        #     recorder_ica_NMDA[syn_id[1]].record(synapse.hsynapse._ref_ica_NMDA)
        # recorder_ica_VDCC[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_ica_VDCC"):
        #     recorder_ica_VDCC[syn_id[1]].record(synapse.hsynapse._ref_ica_VDCC)
    # setup synapses and run simulation
    df = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
    for syn_id, synapse in cell.synapses.items():
        # Reduced logging: removed verbose synapse debug
        # logger.debug("Configuring synapse %d", syn_id[1])
        if syn_extra_params is not None:  # configure local parameters
            _set_local_params(synapse, fit_params,
                              syn_extra_params[df.loc[df["local_syn_idx"] == syn_id[1]].index[0]])
        synapse.hsynapse.rho0_GB = 1  # override rho (not sure if it's needed)
        synapse.hsynapse.Use_p = 1  # override potentiated U_SE (to guarantee release)
        synapse.hsynapse.Use = synapse.hsynapse.Use_p  # override U_SE
        synapse.hsynapse.gmax0_AMPA = synapse.hsynapse.gmax_p_AMPA
        #synapse.hsynapse.gmax_NMDA = synapse.hsynapse.gmax_p_AMPA * 0.55
        synapse.hsynapse.theta_d_GB = -1  # disable LTD
        synapse.hsynapse.theta_p_GB = -1  # disable LTP
    sim.run(1500, cvode=True)
    logger.debug("Simulation completed")
    # get soma voltage and extract spike times
    t = np.array(sim.get_time())
    v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
    # compute c_pre and store results
    c_pre = {df.loc[df["local_syn_idx"] == syn_id].index[0]: recorder[syn_id].max() for syn_id in syn_idx}
    results = {"c_pre": c_pre,
               "c_trace": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
                           np.array(recorder[syn_id].to_python()) for syn_id in syn_idx},
               "cai_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
                           np.array(recorder_cai[syn_id].to_python()) for syn_id in syn_idx},
            #    "cai_NMDA_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
            #                     np.array(recorder_cai_NMDA_CR[syn_id].to_python()) for syn_id in syn_idx},
            #    "cai_VDCC_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
            #                     np.array(recorder_cai_VDCC_CR[syn_id].to_python()) for syn_id in syn_idx},
               "t": t,
               "v": v,
               "t_recorded": np.array(recorder_t.to_python())}
    return results


def c_pre_finder(sim_config, fit_params, syn_extra_params, pre_gid, post_gid,
                 node_pop="S1nonbarrel_neurons", edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical",
                 fixhp=False):
    """Replays spike from `pre_gid` and measures Ca++ transient in synapses on `post_gid`"""
    pool = multiprocessing.Pool(processes=1)
    results = pool.apply(_c_pre_finder_process, [sim_config, fit_params, syn_extra_params, pre_gid, post_gid,
                                                 node_pop, edge_pop, fixhp])
    pool.terminate()
    logger.debug("C_pre: %s", str(results["c_pre"]))
    return results["c_pre"]


def _c_post_finder_process(sim_config, fit_params, syn_extra_params, pre_gid, post_gid, stimulus,
                           node_pop, edge_pop, fixhp):
    """
    Multiprocessing subprocess for `c_post_finder()`
    Injects (precalculated) stimulus to the `post_gid` to make it fire 1 AP and measures the Ca++ transient
    (from the backpropagiting AP) at the synapses made by `pre_gid`
    """
    logger.debug("c_post finder process")
    sim = bluecellulab.CircuitSimulation(sim_config)
    sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, add_minis=False,
                          intersect_pre_gids=[(node_pop, pre_gid)])
    cell = sim.cells[(node_pop, post_gid)]
    if fixhp:  # hyperpolarization workaround
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
    # add stimuli
    tstim = bluecellulab.neuron.h.TStim(0.5, sec=cell.soma)
    stim_duration = (stimulus["nspikes"] - 1) * 1000./stimulus["freq"] + stimulus["width"]
    tstim.train(stimulus["offset"], stim_duration, stimulus["amp"], stimulus["freq"], stimulus["width"])
    cell.persistent.append(tstim)
    if fit_params is not None:  # setup global parameters
        _set_global_params(fit_params)
    # initialize [Ca^{2+}]_i and vsyn recorders
    recorder, recorder_cai, recorder_shaft_cai, recorder_vsyn, recorder_ica = {}, {}, {}, {}, {}
    recorder_ica_NMDA, recorder_ica_VDCC, syn_idx = {}, {}, []
    recorder_cai_NMDA_CR, recorder_cai_VDCC_CR = {}, {}
    recorder_t = bluecellulab.neuron.h.Vector()
    recorder_t.record(bluecellulab.neuron.h._ref_t, 1.0)

    for syn_id, synapse in cell.synapses.items():
        syn_idx.append(syn_id[1])
        recorder[syn_id[1]] = bluecellulab.neuron.h.Vector()
        recorder[syn_id[1]].record(synapse.hsynapse._ref_effcai_GB, 1.0)
        recorder_cai[syn_id[1]] = bluecellulab.neuron.h.Vector()
        recorder_cai[syn_id[1]].record(synapse.hsynapse._ref_cai_CR, 1.0)

        # recorder_cai_NMDA_CR[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_cai_NMDA_CR[syn_id[1]].record(synapse.hsynapse._ref_cai_NMDA_CR, 1.0)
        # recorder_cai_VDCC_CR[syn_id[1]] = bluecellulab.neuron.h.Vector()    
        # recorder_cai_VDCC_CR[syn_id[1]].record(synapse.hsynapse._ref_cai_VDCC_CR, 1.0)
        # recorder_shaft_cai[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_shaft_cai"):
        #     recorder_shaft_cai[syn_id[1]].record(synapse.hsynapse._ref_shaft_cai)
        # # vsyn: local dendritic membrane voltage at the synapse (CVode-adaptive time axis)
        # recorder_vsyn[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_vsyn[syn_id[1]].record(synapse.hsynapse._ref_vsyn)
        # ica: total calcium current density (mA/cm2) at the dendritic segment
        # This is what CaDynamics.mod reads to compute shaft calcium
        # seg = synapse.hsynapse.get_segment()
        # recorder_ica[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # recorder_ica[syn_id[1]].record(seg._ref_ica)
        # # ica_NMDA and ica_VDCC: synapse-level calcium current components
        # recorder_ica_NMDA[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_ica_NMDA"):
        #     recorder_ica_NMDA[syn_id[1]].record(synapse.hsynapse._ref_ica_NMDA)
        # recorder_ica_VDCC[syn_id[1]] = bluecellulab.neuron.h.Vector()
        # if hasattr(synapse.hsynapse, "_ref_ica_VDCC"):
        #     recorder_ica_VDCC[syn_id[1]].record(synapse.hsynapse._ref_ica_VDCC)
    # setup synapses and run simulation
    df = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
    for syn_id, synapse in cell.synapses.items():
        # Reduced logging: removed verbose synapse debug
        # logger.debug("Configuring synapse %d", syn_id[1])
        if syn_extra_params is not None:  # configure local parameters
            _set_local_params(synapse, fit_params,
                              syn_extra_params[df.loc[df["local_syn_idx"] == syn_id[1]].index[0]])
        synapse.hsynapse.theta_d_GB = -1  # disable LTD
        synapse.hsynapse.theta_p_GB = -1  # disable LTP
    sim.run(1500, cvode=True)
    logger.debug("Simulation completed")
    # get soma voltage and simulation time vector and extract spike times
    t = np.array(sim.get_time())
    v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
    spikes = _get_spikes(t, v)
    # compute c_post and store results
    c_post = {df.loc[df["local_syn_idx"] == syn_id].index[0]: recorder[syn_id].max() for syn_id in syn_idx}
    results = {"c_post": c_post,
               "c_trace": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
                           np.array(recorder[syn_id].to_python()) for syn_id in syn_idx},
               "cai_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
                           np.array(recorder_cai[syn_id].to_python()) for syn_id in syn_idx},
            #    "cai_NMDA_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
            #                     np.array(recorder_cai_NMDA_CR[syn_id].to_python()) for syn_id in syn_idx},
            #    "cai_VDCC_CR": {df.loc[df["local_syn_idx"] == syn_id].index[0]:
            #                     np.array(recorder_cai_VDCC_CR[syn_id].to_python()) for syn_id in syn_idx},
               "t": t,
               "v": v,
               "t_spikes": spikes,
               "t_stimuli": np.array(tstim.tvec)[:-1:4],
               "t_recorded": np.array(recorder_t.to_python())}
    return results


def c_post_finder(sim_config, fit_params, syn_extra_params, pre_gid, post_gid, stimulus,
                  node_pop="S1nonbarrel_neurons", edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical",
                  fixhp=False):
    """
    Finds c_post - the calcium transient from a postsynaptic spike in `post_gid` at synapses made by `pre_gid`.
    Injects current that makes the cell fire a single AP and then measures the Ca++ transient of the backpropagating AP.
    (To do so the necessary current that makes `post_gid` fire a single AP must be calculated beforehand
    - see `spike_threshold_finder()` above - and passed as a `stimulus` dictionary.)
    """
    # find c_post
    logger.debug("Stimulating cell with {} nA pulse ({} ms)".format(stimulus["amp"], stimulus["width"]))
    pool = multiprocessing.Pool(processes=1)
    results = pool.apply(_c_post_finder_process, [sim_config, fit_params, syn_extra_params,
                                                  pre_gid, post_gid, stimulus, node_pop, edge_pop, fixhp])
    pool.terminate()
    # validate number of spikes
    logger.debug("Spike timing: {}".format(results["t_spikes"]))
    logger.debug("C_post: %s", str(results["c_post"]))
    if len(results["t_spikes"]) < 1:
        # special case, small integration differences with threshold detection sim
        logger.debug("Cell not spiking as expected during c_post, "
                     "attempting to bump stimulus amplitude before failing...")
        # find c_post
        amp = stimulus["amp"] + 0.05
        logger.debug("Stimulating cell with %f nA pulse", amp)
        stimulus = {"nspikes": 1, "freq": 0.1, "width": stimulus["width"], "offset": 1000., "amp": amp}
        pool = multiprocessing.Pool(processes=1)
        results = pool.apply(_c_post_finder_process, [sim_config, fit_params, syn_extra_params,
                                                      pre_gid, post_gid, stimulus, node_pop, edge_pop, fixhp])
        pool.terminate()
        logger.debug("C_post: %s", str(results["c_post"]))
        return results["c_post"] if len(results["t_spikes"]) == 1 else None
    return results["c_post"]

def _runconnectedpair_prefire_process(conn, workdir, fit_params, syn_extra_params, pre_gid, post_gid, t_end, c_pre, c_post,
                              syn_rec_lst, fastforward, node_pop, edge_pop, fixhp):
    """
    Multiprocessing subprocess for `runconnectedpair()`.
    Sends the results dict through `conn` (a one-way multiprocessing.Pipe connection) on success,
    or closes `conn` and re-raises on failure (causing a non-zero exit code the parent can detect).
    """
    import time
    import logging

    # Set up proper logging for child process
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(f'{workdir}/child_process.log')
        ]
    )
    logger = logging.getLogger(__name__)
    sim_start_time = time.time()

    try:
        sim_config = os.path.join(workdir, "prefire_simulation_config.json")
        logger.info(f"Loading simulation from {sim_config}")
        sim = bluecellulab.CircuitSimulation(sim_config)

        # Get presynaptic spike times from file (heavily relies on naming conventions
        # and the fact that spikes are delivered from a single presynaptic cell)
        pre_spikes = SpikeReader(os.path.join(workdir, "prefire_prespikes.h5"))[node_pop].get_dict()["timestamps"]

        # Instantiate gid with stimuli from sim. config and presynaptic spike times read from file
        sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, add_minis=False, add_pulse_stimuli=True,
                             intersect_pre_gids=[(node_pop, pre_gid)],
                             pre_spike_trains={(node_pop, pre_gid): pre_spikes})
        cell = sim.cells[(node_pop, post_gid)]
        if fixhp:  # hyperpolarization workaround
            for sec in cell.somatic + cell.axonal:
                sec.uninsert("SK_E2")
        if fit_params is not None:  # setup global parameters
            _set_global_params(fit_params)
        syn_rec_lst = SYNREC if syn_rec_lst is None else syn_rec_lst
        # Check what is actually available on the first synapse to record
        actual_syn_rec = []
        if len(cell.synapses) > 0:
            first_syn = list(cell.synapses.values())[0]
            for key in syn_rec_lst:
                if hasattr(first_syn.hsynapse, "_ref_%s" % key):
                    actual_syn_rec.append(key)
            # Add CICR variables if they exist (Li-Rinzel IP3R model, minimal set)
            for k in ["ca_cicr_CICR", "Popen_CICR", "IP3_CICR", "J_cicr_CICR"]:
                if hasattr(first_syn.hsynapse, "_ref_%s" % k) and k not in actual_syn_rec:
                    actual_syn_rec.append(k)

        syn_rec, syn_idx = {key: [] for key in actual_syn_rec}, []
        for syn_id, synapse in cell.synapses.items():
            syn_idx.append(syn_id[1])
            if len(actual_syn_rec) != 0:
                for key, lst in syn_rec.items():  # set up recordings
                    recorder = bluecellulab.neuron.h.Vector()
                    recorder.record(getattr(synapse.hsynapse, "_ref_%s" % key))
                    lst.append(recorder)
        df = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
        syn_props = {key: [] for key in SYNPROPS}
        for syn_id, synapse in cell.synapses.items():
            if syn_extra_params is not None:  # configure local parameters
                global_syn_id = df.loc[df["local_syn_idx"] == syn_id[1]].index[0]
                _set_local_params(synapse, fit_params, syn_extra_params[global_syn_id],
                                  c_pre[global_syn_id], c_post[global_syn_id])
            for key, lst in syn_props.items():  # store synapse properties
                if key == "Cpre":
                    lst.append(c_pre[global_syn_id])
                elif key == "Cpost":
                    lst.append(c_post[global_syn_id])
                elif key == "loc":
                    lst.append(syn_extra_params[global_syn_id]["loc"])
                else:
                    if key in PARAM_MAP:
                        lst.append(getattr(synapse.hsynapse, PARAM_MAP[key]))
                    else:
                        lst.append(getattr(synapse.hsynapse, key))
        # Run
        bluecellulab.neuron.h.cvode_active(1)
        sim.run(t_end, cvode=True)
        logger.debug("Simulation completed")
        # Collect all properties
        syn_props.update({key: getattr(bluecellulab.neuron.h, "%s_GluSynapse" % key) for key in MOD_PROPS if hasattr(bluecellulab.neuron.h, "%s_GluSynapse" % key)})
        # Collect Results
        t = np.array(sim.get_time())
        v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
        results = {"c_pre": c_pre, "c_post": c_post, "c_trace": None,
                   "t": t, "syn_idx": syn_idx, "syn_props": syn_props,
                   "v": v, "prespikes": pre_spikes, "postspikes": _get_spikes(t, v)}
        for key in syn_rec_lst:
            if key in syn_rec:
                results[key] = {df.loc[df["local_syn_idx"] == syn_id].index[0]:
                                np.array(syn_rec[key][i].to_python()) for i, syn_id in enumerate(syn_idx)}
        sim_duration = time.time() - sim_start_time
        logger.info(f"Simulation walltime: {sim_duration:.2f} seconds")
        # Send results back to parent through the pipe
        conn.send(results)
    except Exception as e:
        logger.error(f"Child process failed: {e}")
        raise
    finally:
        conn.close()


def _runconnectedpair_process(results, workdir, fit_params, syn_extra_params, pre_gid, post_gid, t_end, c_pre, c_post,
                              syn_rec_lst, fastforward, node_pop, edge_pop, fixhp, log_queue=None):
    """
    Multiprocessing subprocess for `runconnectedpair()`
    Injects periodic current pulses (read from simulation_config.json) that makes the postsynaptic cell fire APs,
    while delivering spikes from the (non-simulated) presynaptic cell (read from file
    written by `plastyfire/simwriter` as well) after setting up all custom synapse parameters
    """
    import time
    sim_start_time = time.time()
    # Configure logging for child process
    import logging
    child_logger = logging.getLogger(__name__)
    if not child_logger.handlers:
        # Copy logging configuration from parent
        parent_logger = logging.getLogger()
        for handler in parent_logger.handlers:
            child_logger.addHandler(handler)
        child_logger.setLevel(parent_logger.level)
    sim_config = os.path.join(workdir, "simulation_config.json")
    sim = bluecellulab.CircuitSimulation(sim_config)
    if log_queue:
        log_queue.put(f"Loaded simulation from {sim_config}")
    logger.debug("Loaded simulation")
    # Get presynaptic spike times from file (heavily relies on naming conventions
    # and the fact that spikes are delivered from a single presynaptic cell)
    pre_spikes = SpikeReader(os.path.join(workdir, "prespikes.h5"))[node_pop].get_dict()["timestamps"]
    # Instantiate gid with stimuli from sim. config and presynaptic spike times read from file
    sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, add_minis=False, add_pulse_stimuli=True,
                         intersect_pre_gids=[(node_pop, pre_gid)],
                         pre_spike_trains={(node_pop, pre_gid): pre_spikes})
    cell = sim.cells[(node_pop, post_gid)]
    if fixhp:  # hyperpolarization workaround
        for sec in cell.somatic + cell.axonal:
            sec.uninsert("SK_E2")
    if fit_params is not None:  # setup global parameters
        _set_global_params(fit_params)
    syn_rec_lst = SYNREC if syn_rec_lst is None else syn_rec_lst
    # Check what is actually available on the first synapse to record
    actual_syn_rec = []
    if len(cell.synapses) > 0:
        first_syn = list(cell.synapses.values())[0]
        for key in syn_rec_lst:
            if hasattr(first_syn.hsynapse, "_ref_%s" % key):
                actual_syn_rec.append(key)
        # Add V7 specific variables if they exist
        for k in ["IP3_CICR", "P_CICR", "ca_ryr_CICR", "ca_ip3r_CICR", "S_CICR", "T_CICR", "J_serca_CICR", "J_leak_CICR", "ip3_gate_CICR", "ca_gate_CICR", "J_ryr_CICR", "J_ip3r_CICR", "Nves_CICR"]:
            if hasattr(first_syn.hsynapse, "_ref_%s" % k) and k not in actual_syn_rec:
                actual_syn_rec.append(k)

    syn_rec, syn_idx = {key: [] for key in actual_syn_rec}, []
    for syn_id, synapse in cell.synapses.items():
        syn_idx.append(syn_id[1])
        if len(actual_syn_rec) != 0:
            for key, lst in syn_rec.items():  # set up recordings
                recorder = bluecellulab.neuron.h.Vector()
                recorder.record(getattr(synapse.hsynapse, "_ref_%s" % key))
                lst.append(recorder)
    df = _map_syn_idx(sim_config, post_gid, syn_idx, edge_pop)
    syn_props = {key: [] for key in SYNPROPS}
    for syn_id, synapse in cell.synapses.items():
        # Reduced logging: removed verbose synapse debug
        # logger.debug("Configuring synapse %d", syn_id[1])
        if syn_extra_params is not None:  # configure local parameters
            global_syn_id = df.loc[df["local_syn_idx"] == syn_id[1]].index[0]
            _set_local_params(synapse, fit_params, syn_extra_params[global_syn_id],
                              c_pre[global_syn_id], c_post[global_syn_id])
        for key, lst in syn_props.items():  # store synapse properties
            if key == "Cpre":
                lst.append(c_pre[global_syn_id])
            elif key == "Cpost":
                lst.append(c_post[global_syn_id])
            elif key == "loc":
                lst.append(syn_extra_params[global_syn_id]["loc"])
            else:
                if key in PARAM_MAP:
                    lst.append(getattr(synapse.hsynapse, PARAM_MAP[key]))
                else:
                    lst.append(getattr(synapse.hsynapse, key))
        # for attr in dir(synapse.hsynapse):  # show all params
        #     if re.match('__.*', attr) is None:
        #         logger.debug("%s = %s", attr, str(getattr(synapse.hsynapse, attr)))
    # Run
    t_end = 3 if DEBUG else t_end
    if (fastforward is not None):
        # Run until fastforward point
        if log_queue:
            log_queue.put(f"Fastforward enabled, simulating {fastforward/1000.:.1f} seconds...")
        # Reduced logging: removed verbose simulation debug
    # logger.debug("Fastforward enabled, simulating %.1f seconds...", fastforward / 1000.)
        sim.run(fastforward, cvode=True)
        # Fastforward synapses
        if log_queue:
            log_queue.put("Updating synapses...")
        # Reduced logging: removed verbose synapse debug
        # logger.debug("Updating synapses...")
        for syn_id, synapse in cell.synapses.items():
            # Reduced logging: removed verbose synapse debug
        # logger.debug("Configuring synapse %d", syn_id[1])
            if synapse.hsynapse.rho_GB >= 0.5:
                synapse.hsynapse.rho_GB = 1.
                synapse.hsynapse.Use = synapse.hsynapse.Use_p
                synapse.hsynapse.Use_GB = synapse.hsynapse.Use_p
                synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_p_AMPA
                synapse.hsynapse.gmax0_AMPA = synapse.hsynapse.gmax_p_AMPA
            else:
                synapse.hsynapse.rho_GB = 0.
                synapse.hsynapse.Use = synapse.hsynapse.Use_d
                synapse.hsynapse.Use_GB = synapse.hsynapse.Use_d
                synapse.hsynapse.gmax_AMPA = synapse.hsynapse.gmax_d_AMPA
                synapse.hsynapse.gmax0_AMPA = synapse.hsynapse.gmax_d_AMPA
        # Complete run
        if log_queue:
            log_queue.put(f"Simulating remaining {(t_end - fastforward)/1000.:.1f} seconds...")
        # Reduced logging: removed verbose simulation debug
        # logger.debug("Simulating remaining %.1f seconds...", (t_end - fastforward) / 1000.)
        bluecellulab.neuron.h.cvode_active(1)
        bluecellulab.neuron.h.continuerun(t_end)
    else:
        if log_queue:
            log_queue.put(f"Simulating {t_end/1000.:.1f} seconds...")
        # Reduced logging: removed verbose simulation debug
        # logger.debug("Simulating %.1f seconds...", t_end / 1000.)
        sim.run(t_end, cvode=True)
    if log_queue:
        log_queue.put("Simulation completed")
    logger.debug("Simulation completed")
    # Collect all properties
    syn_props.update({key: getattr(bluecellulab.neuron.h, "%s_GluSynapse" % key) for key in MOD_PROPS if hasattr(bluecellulab.neuron.h, "%s_GluSynapse" % key)})
    # Collect Results
    t = np.array(sim.get_time())
    v = np.array(sim.get_voltage_trace((node_pop, post_gid)))
    results["t"] = t
    results["v"] = v
    results["prespikes"] = pre_spikes
    results["postspikes"] = _get_spikes(t, v)
    results["synprop"] = syn_props
    if len(syn_rec_lst) != 0:
        for key, lst in syn_rec.items():
            results[key] = np.transpose([np.array(rec.to_python()) for rec in lst])
    
    # Log total simulation walltime
    sim_duration = time.time() - sim_start_time
    if log_queue:
        log_queue.put(f"Simulation walltime: {sim_duration:.2f} seconds")
    # Reduced logging: only log walltime if needed
    # logger.debug("Simulation walltime: %.2f seconds", sim_duration)


def runconnectedpair_induction(workdir, fit_params=None, syn_rec_lst=None, fastforward=None,
                               node_pop="S1nonbarrel_neurons",
                               edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical",
                               fixhp=True, recipe_path=None, epg_variant="epg_dhuruva"):
    """
    Run the induction simulation using `_runconnectedpair_process`.

    Reads from ``simulation_config.json`` and ``prespikes.h5`` (as opposed to the
    prefire variants).  Supports the ``fastforward`` parameter to fast-forward
    the simulation to a given point (ms) before the recording window.
    """
    sim_config = os.path.join(workdir, "simulation_config.json")
    sim = Simulation(sim_config)
    pre_gid = sim.node_sets.content["precell"]["node_id"][0]
    post_gid = sim.node_sets.content["postcell"]["node_id"][0]
    t_end = sim.config["run"]["tstop"]
    width = sim.config["inputs"]["pulse0"]["width"]
    amp = sim.config["inputs"]["pulse0"]["amp_start"]
    stimulus = {"nspikes": 1, "freq": 0.1, "width": width, "offset": 1000, "amp": amp}

    recipe_file = recipe_path if recipe_path else EXTRA_RECIPE_PATH
    params_generator_cls = _get_params_generator(epg_variant)
    pgen = params_generator_cls(sim.circuit, node_pop, edge_pop, recipe_file)
    syn_extra_params = pgen.generate_params(pre_gid, post_gid)

    baseline_params = fit_params.copy() if fit_params is not None else {}
    c_pre = c_pre_finder(sim_config, baseline_params, syn_extra_params, pre_gid, post_gid,
                         node_pop=node_pop, edge_pop=edge_pop, fixhp=fixhp)
    c_post = c_post_finder(sim_config, baseline_params, syn_extra_params, pre_gid, post_gid, stimulus,
                           node_pop=node_pop, edge_pop=edge_pop, fixhp=fixhp)

    with multiprocessing.Manager() as manager:
        results = manager.dict()
        proc = multiprocessing.Process(
            target=_runconnectedpair_process,
            args=(results, workdir, fit_params, syn_extra_params, pre_gid, post_gid, t_end,
                  c_pre, c_post, syn_rec_lst, fastforward, node_pop, edge_pop, fixhp))
        proc.start()
        proc.join()
        if proc.exitcode != 0:
            raise RuntimeError(f"Child process failed with exit code {proc.exitcode}")
        return dict(results)  # copy out before manager shuts down


def runconnectedpair(workdir, fit_params=None, syn_rec_lst=None, fastforward=None,
                     node_pop="S1nonbarrel_neurons", edge_pop="S1nonbarrel_neurons__S1nonbarrel_neurons__chemical",
                     fixhp=True, recipe_path=None, epg_variant="epg_dhuruva"):
    """
    Pairs spikes (at given frequency and dt written in `plastyfire/simwriter`) and triggers EPSP before and after that
    (to test the effect of the STDP protocol on the EPSP amplitude)
    """
    # Read necessary params from the sim. config
    sim_config = os.path.join(workdir, "simulation_config.json")
    sim = Simulation(sim_config)
    pre_gid, post_gid = sim.node_sets.content["precell"]["node_id"][0], sim.node_sets.content["postcell"]["node_id"][0]
    t_end = sim.config["run"]["tstop"]
    width, amp = sim.config["inputs"]["pulse0"]["width"], sim.config["inputs"]["pulse0"]["amp_start"]
    stimulus = {"nspikes": 1, "freq": 0.1, "width": width, "offset": 1000, "amp": amp}
    # Get reference Cpre and Cpost values (used to derive depression and potentiation thresholds)
    
    # Use custom recipe path if provided, otherwise use default
    recipe_file = recipe_path if recipe_path else EXTRA_RECIPE_PATH
    params_generator_cls = _get_params_generator(epg_variant)
    pgen = params_generator_cls(sim.circuit, node_pop, edge_pop, recipe_file)
    syn_extra_params = pgen.generate_params(pre_gid, post_gid)

    baseline_params = fit_params.copy() if fit_params is not None else {}

    c_pre = c_pre_finder(sim_config, baseline_params, syn_extra_params, pre_gid, post_gid,
                         node_pop=node_pop, edge_pop=edge_pop, fixhp=fixhp)
    c_post = c_post_finder(sim_config, baseline_params, syn_extra_params, pre_gid, post_gid, stimulus,
                           node_pop=node_pop, edge_pop=edge_pop, fixhp=fixhp)

    # --- Diagnostic: print calibration values before induction starts ---
    print("=" * 60)
    print(f"[runconnectedpair] Calibration values for {workdir}")
    if fit_params:
        print(f"  tau_effca_GB_GluSynapse = {fit_params.get('tau_effca_GB_GluSynapse', 'NOT SET')}")
    for syn_id in sorted(c_pre.keys()):
        cp  = c_pre[syn_id]
        cpo = c_post.get(syn_id, float("nan")) if c_post else float("nan")
        # Compute thresholds using the a-values in fit_params (if present)
        td = tp = None
        if fit_params:
            loc = syn_extra_params[syn_id].get("loc", "basal") if syn_extra_params else "basal"
            if loc == "basal" and all(k in fit_params for k in ("a00", "a01", "a10", "a11")):
                td = fit_params["a00"] * cp + fit_params["a01"] * cpo
                tp = fit_params["a10"] * cp + fit_params["a11"] * cpo
            elif loc == "apical" and all(k in fit_params for k in ("a20", "a21", "a30", "a31")):
                td = fit_params["a20"] * cp + fit_params["a21"] * cpo
                tp = fit_params["a30"] * cp + fit_params["a31"] * cpo
        td_str = f"{td:.6f}" if td is not None else "N/A (no a-values)"
        tp_str = f"{tp:.6f}" if tp is not None else "N/A (no a-values)"
        print(f"  syn {syn_id}: cpre={cp:.6f}  cpost={cpo:.6f}  theta_d={td_str}  theta_p={tp_str}")
    print("=" * 60)
    # --- end diagnostic ---

    # Run main simulation
    # Reduced logging: removed verbose simulation info
    # logger.info("Simulating %s...", workdir)
    # Use a one-way Pipe instead of Manager: no extra server process, results freed when child exits
    parent_conn, child_conn = multiprocessing.Pipe(duplex=False)
    child_proc = multiprocessing.Process(target=_runconnectedpair_prefire_process,
                                         args=(child_conn, workdir, fit_params, syn_extra_params, pre_gid, post_gid, t_end,
                                               c_pre, c_post, syn_rec_lst, fastforward, node_pop, edge_pop, fixhp))
    child_proc.start()
    child_conn.close()  # close child end in parent — required so parent recv() sees EOF on child crash
    try:
        results = parent_conn.recv()  # blocks until child sends or closes conn
    except EOFError:
        results = None  # child crashed before sending — exit code check below will raise
    finally:
        parent_conn.close()
    child_proc.join()
    if child_proc.exitcode != 0:
        raise RuntimeError(f"Child process failed with exit code {child_proc.exitcode}")
    return results
