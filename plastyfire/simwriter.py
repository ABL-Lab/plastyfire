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
from plastyfire.simulator import spike_threshold_finder, runsinglecell, doublet_tolerant, tolerant_amp_grid, any_firing
from plastyfire.pipette import afferent_table, rank_pre_gids, measure_epsp, calibrate_group

MECHANISMS_PATH = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"  # same as simulator_edges

MIN2MS = 60 * 1000.
OPT_CPU_TIME = 2.  # heuristics: it takes ~2x compute time (w/ CVode w/ reporting w/o fastforward) as biological time
CPU_TIME = 1.  # heuristics: c_pre and c_post for a single connection takes <1 minute to simulate/calculate
FIGS_DIR = "/project/rrg-emuller/dhuruva/figures/plastyfire"

# Single-AP threshold search (for c_post), mirroring thresholdfinder.py::run()
# and simulator_edges.C_POST_*. Distinct from the induction-train search, which
# takes its amplitude grid from the stimulus config.
# Widths stop at 5 ms on purpose: spike_threshold_finder's binary search assumes
# the leftmost amplitude firing >=nspikes also fires ==nspikes, which for
# nspikes=1 holds only while the pulse is short enough that rheobase gives a
# single AP rather than a train. Do not extend this list upward.
SINGLE_AP_PULSE_WIDTHS_MS = [1.5, 3, 5]
SINGLE_AP_FREQ_HZ         = 0.1
SINGLE_AP_MIN_AMP_NA      = 0.05
SINGLE_AP_MAX_AMP_NA      = 5.
SINGLE_AP_AMP_LEVELS      = 100   # -> 0.05 nA grid over [0.05, 5.] nA


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
    (and save params. of current injection that makes it fire)

    Two searches are run per cell and both are stored in the same pkl:

    - one per entry in `stim_config["freq"]`, keyed by that float frequency,
      for the `nspikes`-spike induction train. This is the one that gates pair
      selection: no amplitude found -> the cell is rejected (return False).
    - one under the "single_ap" key, for the single AP that c_post needs
      (`_find_cpre_cpost` in simulator_edges.py). Cached here so it is paid
      once per cell at calibration instead of once per pair on every rebuild
      of the c_pre/c_post cache — tau_effca / gamma_d / gamma_p are GluSynapse
      params and do not change whether the soma fires.

    BOTH searches gate pair selection. A cell that cannot be made to fire
    exactly one AP has no valid c_post, so `_find_cpre_cpost` raises on it and
    the pair is lost anyway — but only after the basis and the whole c_pre/c_post
    sweep have been paid for. Rejecting it here costs one extra search per cell
    and keeps the pair list to cells the full pipeline can actually use.
    """
    pklf_name = os.path.join(save_dir, "%i.pkl" % gid)
    if os.path.isfile(pklf_name):  # check if these sims were already run...
        with open(pklf_name, "rb") as f:
            cached = pickle.load(f)
        # A pkl written before the single-AP gate existed carries only the
        # induction-train entries. Its presence proves the cell fires `nspikes`
        # APs at the induction frequency -- NOT that c_post is defined. gids
        # 185042 and 197123 passed on exactly that and then failed later in
        # `_find_cpre_cpost` ("could not fire ... exactly one AP at any width
        # up to 5.0 nA"), costing 2 of the 100 pairs. Re-run the searches
        # rather than trusting the file's existence.
        if "single_ap" in cached and set(stim_config["freq"]).issubset(cached):
            return True
        warnings.warn("Stale stimulus cache for gid %i (keys %s) - re-running searches"
                      % (gid, sorted(cached, key=str)))
    results = {}
    for freq in stim_config["freq"]:
        simres = spike_threshold_finder(sim_config, gid, stim_config["nspikes"],
                                        freq, stim_config["width"], stim_config["offset"], stim_config["amp_min"],
                                        stim_config["amp_max"], stim_config["amp_lev"], fixhp=True)
        if simres is None and doublet_tolerant():  # opt-in fallback: 1 AP or a <= 6 ms doublet per pulse
            simres = spike_threshold_finder(sim_config, gid, stim_config["nspikes"],
                                            freq, stim_config["width"], stim_config["offset"],
                                            *tolerant_amp_grid(stim_config["amp_min"], stim_config["amp_max"],
                                                               stim_config["amp_lev"]), fixhp=True, doublet_ok=True)
        if simres is None:
            return False
        else:
            results[freq] = simres
    # Single-AP threshold for c_post. Widths / amplitude range / level count
    # mirror thresholdfinder.py::run() and simulator_edges.C_POST_* — keep them
    # in sync. The 1-AP amplitude is lower than the induction one, so the
    # induction amp_min/amp_lev grid is too coarse and starts too high to use.
    for pulse_width in SINGLE_AP_PULSE_WIDTHS_MS:
        simres = spike_threshold_finder(sim_config, gid, 1, SINGLE_AP_FREQ_HZ,
                                        pulse_width, stim_config["offset"], SINGLE_AP_MIN_AMP_NA,
                                        SINGLE_AP_MAX_AMP_NA, SINGLE_AP_AMP_LEVELS, fixhp=True)
        if simres is not None:
            results["single_ap"] = simres
            break
    if "single_ap" not in results and doublet_tolerant():  # opt-in fallback, after every strict width failed
        for pulse_width in SINGLE_AP_PULSE_WIDTHS_MS:
            simres = spike_threshold_finder(sim_config, gid, 1, SINGLE_AP_FREQ_HZ,
                                            pulse_width, stim_config["offset"], SINGLE_AP_MIN_AMP_NA,
                                            SINGLE_AP_MAX_AMP_NA, SINGLE_AP_AMP_LEVELS, fixhp=True, doublet_ok=True)
            if simres is not None:
                results["single_ap"] = simres
                break
    if "single_ap" not in results:
        # No single AP at any width -> c_post is undefined for this cell. Reject
        # it now rather than letting _find_cpre_cpost raise on it later.
        warnings.warn("No single-AP stimulus found for gid %i at any width in %s ms "
                      "up to %.1f nA - rejecting cell (c_post would be invalid)"
                      % (gid, SINGLE_AP_PULSE_WIDTHS_MS, SINGLE_AP_MAX_AMP_NA))
        return False
    # Store results for manual validation and simulation setup
    with open(pklf_name, "wb") as f:
        pickle.dump(results, f, -1)
    return True


def protocol_stim_key(proto, stim_config):
    """Key of the post stimulus of `proto` in the single_cells pkl. Protocols of the classic freq x dt grid
    keep the float frequency key written by `check_electrical_constraint` (searched in `find_pairs()` on
    `pre_simwriter_config`), so existing caches stay valid. Explicit `protocols:` always get their own
    string key, so they are calibrated on `calibration_circuit_config` even when the same pkl already has
    a grid entry for that frequency. None if there is no post stimulus."""
    if proto["n_post"] == 0:
        return None
    if proto.get("post_type", "pulse") == "step":
        return step_stim_key(proto, stim_config)
    if proto["grid"]:
        return float(proto["freq"])
    key = "%iap_%gHz_%gms" % (proto["n_post"], proto["freq"], proto["width"])
    if proto["amp"] is not None:
        return key + "_%gnA" % proto["amp"]
    amp_range = [proto.get(k, stim_config[k]) for k in ("amp_min", "amp_max", "amp_lev")]
    if amp_range != [stim_config[k] for k in ("amp_min", "amp_max", "amp_lev")]:
        key += "_search%g-%gnA_%i" % tuple(amp_range)
    return key + "_burstwin" if proto.get("burst_window_search") else key


def step_stim_key(proto, stim_config):
    """Key of a depolarising step (`post_type` / `pre_type` "step") in the single_cells pkl. The same key is
    used for the pre and the post cell, so a gid that is pre in one pair and post in another is searched once."""
    key = "step_%gms" % proto["step_duration"]
    if proto["amp"] is not None:
        key += "_%gnA" % proto["amp"]
    else:
        key += "_search%g-%gnA_%i" % tuple(proto.get(k, stim_config[k]) for k in ("amp_min", "amp_max", "amp_lev"))
    key += "_%i-%iap" % (proto["step_nspikes"], proto.get("step_nspikes_max", proto["step_nspikes"]))
    if proto.get("step_below_nA"):  # subthreshold step (Sjostrom 2004 dLTD): threshold amplitude less step_below_nA
        key += "_below%gnA" % proto["step_below_nA"]
    return key


def protocol_pre_stim_key(proto, stim_config):
    """Key of the presynaptic step of `proto` (`pre_type: step`) in the pre gid's single_cells pkl, else None"""
    if proto["n_pre"] == 0 or proto.get("pre_type", "train") != "step":
        return None
    return step_stim_key(proto, stim_config)


def step_rate_threshold_finder(sim_config, gid, duration, offset, min_amp, max_amp, nlevels, nspikes_min,
                               nspikes_max, node_pop="S1nonbarrel_neurons", fixhp=True, below_nA=0.):
    """Lowest amplitude on np.linspace(min_amp, max_amp, nlevels) at which a single `duration` ms step fires
    >= `nspikes_min` APs within [onset, onset + duration + 5 ms] (binary search; the step spike count grows with
    the amplitude), accepted if that count is <= `nspikes_max`. With nlevels 1 this only checks a fixed
    amplitude. Returns the `runsinglecell` results (+ "t_step_spikes", "n_step_spikes") or None."""
    candidate_amp = np.linspace(min_amp, max_amp, nlevels)
    simres = {}

    def n_step(m):
        if m not in simres:
            stim = {"nspikes": 1, "freq": SINGLE_AP_FREQ_HZ, "width": duration, "offset": offset,
                    "amp": candidate_amp[m]}
            simres[m] = runsinglecell(sim_config, gid, stim, node_pop, fixhp)
            spikes = simres[m]["t_spikes"]
            simres[m]["t_step_spikes"] = spikes[(spikes >= offset) & (spikes <= offset + duration + 5.)]
            simres[m]["n_step_spikes"] = len(simres[m]["t_step_spikes"])
        return simres[m]["n_step_spikes"]

    L, R = 0, nlevels
    while L < R:
        m = int(np.floor((L + R) / 2.))
        if n_step(m) < nspikes_min:
            L = m + 1
        else:
            R = m
    if L < nlevels and nspikes_min <= n_step(L) <= nspikes_max:
        if below_nA:
            # subthreshold step (Sjostrom 2004: "minimal suprathreshold current less 20 pA"): the grid level
            # `below_nA` under the threshold level; it must fire no AP in the step (else None = skipped)
            m = int(np.argmin(np.abs(candidate_amp - (candidate_amp[L] - below_nA))))
            if m < L and n_step(m) == 0:
                return simres[m]
            return None
        return simres[L]
    return None


def burst_window_threshold_finder(sim_config, post_gid, nspikes, freq, width, offset, min_amp, max_amp, nlevels,
                                  node_pop="S1nonbarrel_neurons", fixhp=True, max_scan=10):
    """Like `spike_threshold_finder`, but only APs inside the burst window [first pulse onset, last pulse onset +
    `width` + 5 ms] are counted, and later spikes are allowed. For high-frequency bursts (Letzkus 2006: 3 APs at
    200 Hz) the total spike count is not monotonic in amplitude: low amplitudes evoke late (Ca2+ spike driven)
    spikes, so the plain search lands on a low amplitude and its exact-count check rejects the cell.
    Binary search for the lowest amplitude with >= `nspikes` APs in the window, then scan up to `max_scan` levels
    upwards for one with exactly `nspikes` in-window APs, each after its own pulse onset. Returns the
    `runsinglecell` results (all spikes in "t_spikes", in-window ones first as they are sorted) or None."""
    candidate_amp = np.linspace(min_amp, max_amp, nlevels)
    t_stim = offset + 1000. / freq * np.arange(nspikes)
    t_end = t_stim[-1] + width + 5.
    simres = {}

    def n_in_window(m):
        if m not in simres:
            stim = {"nspikes": nspikes, "freq": freq, "width": width, "offset": offset, "amp": candidate_amp[m]}
            simres[m] = runsinglecell(sim_config, post_gid, stim, node_pop, fixhp)
        spikes = simres[m]["t_spikes"]
        return np.sum((spikes >= t_stim[0]) & (spikes <= t_end))

    L, R = 0, nlevels
    while L < R:
        m = int(np.floor((L + R) / 2.))
        if n_in_window(m) < nspikes:
            L = m + 1
        else:
            R = m
    for m in range(L, min(L + max_scan, nlevels)):
        if n_in_window(m) == nspikes:
            spikes = simres[m]["t_spikes"]
            if np.all(spikes[:nspikes] >= t_stim):
                return simres[m]
    return None


def calibrate_protocol_stimuli(sim_config, gid, protocols, stim_config, save_dir, pre_protocols=()):
    """Makes sure the single_cells pkl of `gid` has an entry for every post stimulus used in `protocols`
    (and, as presynaptic cell, for every pre step of `pre_protocols`, see `protocol_pre_stim_key`).
    Missing ones are searched (`amp` None) or checked (fixed `amp`: the cell has to fire exactly `n_post`
    APs at it). A failed search is stored as None, so it is not retried and the protocol is skipped for
    this cell in `write_sim_files`. Returns the stimulus keys that failed."""
    pklf_name = os.path.join(save_dir, "%i.pkl" % gid)
    cached = {}
    if os.path.isfile(pklf_name):
        with open(pklf_name, "rb") as f:
            cached = pickle.load(f)
    updated = False
    jobs = [(p, protocol_stim_key(p, stim_config)) for p in protocols] + \
           [(p, protocol_pre_stim_key(p, stim_config)) for p in pre_protocols]
    tolerant = doublet_tolerant()
    for proto, key in jobs:
        # a failed (None) search is retried only by the opt-in doublet-tolerant fallback
        retry = tolerant and key in cached and cached[key] is None
        # any_firing(): a stimulus that an earlier tolerant run scanned up above threshold is searched again
        retry = retry or (any_firing() and isinstance(cached.get(key), dict)
                          and cached[key].get("threshold_amp", cached[key]["amp"]) != cached[key]["amp"])
        if key is None or (key in cached and not retry):
            continue
        freq = proto["freq"] if proto["freq"] > 0 else SINGLE_AP_FREQ_HZ
        if proto["amp"] is None:
            amp_min, amp_max, amp_lev = [proto.get(k, stim_config[k]) for k in ("amp_min", "amp_max", "amp_lev")]
        else:
            amp_min, amp_max, amp_lev = proto["amp"], proto["amp"], 1
        if isinstance(key, str) and key.startswith("step_"):  # grid keys are floats
            cached[key] = step_rate_threshold_finder(sim_config, gid, proto["step_duration"], stim_config["offset"],
                                                     amp_min, amp_max, amp_lev, proto["step_nspikes"],
                                                     proto.get("step_nspikes_max", proto["step_nspikes"]),
                                                     below_nA=proto.get("step_below_nA", 0.))
        else:
            finder = burst_window_threshold_finder if proto.get("burst_window_search") else spike_threshold_finder
            res = None if retry else finder(sim_config, gid, proto["n_post"], freq, proto["width"],
                                            stim_config["offset"], amp_min, amp_max, amp_lev, fixhp=True)
            if res is None and tolerant and finder is spike_threshold_finder and proto["n_post"]:
                grid = (amp_min, amp_max, amp_lev) if amp_lev == 1 else tolerant_amp_grid(amp_min, amp_max, amp_lev)
                res = spike_threshold_finder(sim_config, gid, proto["n_post"], freq, proto["width"],
                                             stim_config["offset"], *grid, fixhp=True, doublet_ok=True)
            cached[key] = res
        updated = True
    if updated:
        with open(pklf_name, "wb") as f:
            pickle.dump(cached, f, -1)
    return [key for key in {key for _, key in jobs} if key is not None and cached[key] is None]


def step_pairing_times(proto, simres, base, pre_simres):
    """`pairing_times` for step protocols (`post_type` / `pre_type` "step"). The post step starts at `base`;
    a pre step starts `dt` ms before the post step (dt_ref "stim") and the pre spikes are the APs the pre cell
    fired during its own calibration step (`pre_simres`), shifted to that onset. A pre spike train
    (`pre_type` "train") next to a post step starts `dt` before the reference as in `pairing_times`."""
    post_onsets = np.array([base] if proto["n_post"] else [], dtype=np.float64)
    if proto["n_pre"] and proto.get("pre_type", "train") == "step":
        assert proto["dt_ref"] == "stim", "pre steps are timed to the post step onset (dt_ref: stim) (%s)" % proto["id"]
        pre = base - proto["dt"] + (pre_simres["t_step_spikes"] - pre_simres["t_stimuli"][0])
        return post_onsets, np.asarray(pre, dtype=np.float64)
    if proto["n_post"] and proto["dt_ref"].startswith("ap"):
        ref = base + (simres["t_step_spikes"] - simres["t_stimuli"][0])
        ref = ref[::-1] if proto["dt_ref"].endswith("_last") else ref
    else:
        ref = post_onsets
    t0 = ref[0] if len(ref) else base
    pre_isi = 1000.0 / proto["pre_freq"] if proto["pre_freq"] > 0 else 0.
    return post_onsets, np.array([t0 - proto["dt"] + k * pre_isi for k in range(proto["n_pre"])])


def pairing_times(proto, simres, base, pre_simres=None):
    """Post current pulse onsets and pre spike times of one pairing repetition starting at `base` (ms).
    AP times are pulse onsets + the AP delays measured in the calibration run (`simres`).
    Step protocols go to `step_pairing_times` (`pre_simres`: the pre cell's step calibration)."""
    if "step" in (proto.get("post_type"), proto.get("pre_type")):
        return step_pairing_times(proto, simres, base, pre_simres)
    isi = 1000.0 / proto["freq"] if proto["freq"] > 0 else 0.
    post_onsets = np.array([base + i * isi for i in range(proto["n_post"])])
    if proto["n_post"]:
        # slice both: burst-window calibrations may carry later (Ca2+ burst) spikes after the n_post APs
        spike_delay = simres["t_spikes"][:proto["n_post"]] - simres["t_stimuli"][:proto["n_post"]]
        post_aps = post_onsets + spike_delay
    else:
        spike_delay, post_aps = None, np.array([])
    ref = post_aps if proto["dt_ref"].startswith("ap") and proto["n_post"] else post_onsets
    if proto["dt_ref"].endswith("_last"):
        ref = ref[::-1]
    if proto["n_pre"] == proto["n_post"] and proto["pre_freq"] == proto["freq"] and proto["n_pre"] \
            and proto["dt_ref"] in ("ap", "stim"):
        # pre spike i paired with post AP i (Markram 1997); same operation order as the original simwriter
        pre = post_onsets - proto["dt"] + (spike_delay[:proto["n_post"]] if proto["dt_ref"] == "ap" else 0.)
    else:
        t0 = ref[0] if len(ref) else base
        pre_isi = 1000.0 / proto["pre_freq"] if proto["pre_freq"] > 0 else 0.
        pre = np.array([t0 - proto["dt"] + k * pre_isi for k in range(proto["n_pre"])])
    return post_onsets, pre


def save_spikes(h5f_name, prefix, spike_times, spiking_gids):
    """Save spikes to SONATA format"""
    assert (spiking_gids.shape == spike_times.shape)
    order = np.argsort(spike_times, kind="stable")  # pipette groups: one train per gid, merged in time
    spike_times, spiking_gids = spike_times[order], spiking_gids[order]
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
    
    def _pipette_post_gid(self, args):
        """Worker function for `find_pipette_groups()`: electrical constraint, pipette site and compound EPSP
        calibration for one post_gid. Everything is logged to `single_cells/pipette_<gid>.json` (reused if present)"""
        post_gid, i, sim_config_path, save_dir = args
        logf = os.path.join(save_dir, "pipette_%i.json" % post_gid)
        if os.path.isfile(logf):
            with open(logf) as f:
                log = json.load(f)
            return (tuple(log["pre_gids"]), post_gid) if log["ok"] else None
        if not check_electrical_constraint(sim_config_path, post_gid, self.config["stimulus"], save_dir):
            return None
        pip = self.pipette
        df = afferent_table(Circuit(self.circuit_config), self.node_pop, self.edge_pop, post_gid)
        site, ranked = rank_pre_gids(df, np.random.default_rng(self.seed + i), pip.get("path_dist", [50., 150.]))
        log = {"post_gid": int(post_gid), "ok": False, "site": site, "pre_gids": []}
        if site is not None:
            # EPSP sims: GluSynapse on every afferent (paircells = whole population), same conditions as prefire
            pip_dir = os.path.join(save_dir, "pipette")
            os.makedirs(pip_dir, exist_ok=True)
            ns_name = os.path.join(pip_dir, "%i_node_sets.json" % post_gid)
            with open(ns_name, "w", encoding="utf-8") as f:
                json.dump({"postcell": {"node_id": [int(post_gid)], "population": self.node_pop},
                           "paircells": {"population": self.node_pop}}, f, indent=4)
            spike_times = [1000. * (k + 1) for k in range(pip.get("n_epsp_spikes", 5))]
            cfg_name = os.path.join(pip_dir, "%i_epsp_config.json" % post_gid)
            with open(cfg_name, "w", encoding="utf-8") as f:
                json.dump(self._sim_config(pip_dir, ns_name, spike_times[-1] + 200., {}, 0., 0.), f, indent=4)
            ranked_gids = ranked.index.to_numpy(dtype=int)
            pre_gids, res, trials = calibrate_group(
                ranked_gids, lambda gids: measure_epsp(cfg_name, post_gid, gids, self.node_pop, spike_times,
                                                       MECHANISMS_PATH),
                pip.get("epsp_target", 2.), pip.get("epsp_range", [1., 3.]), pip.get("max_pre", 200))
            pre_gids = np.array(pre_gids, dtype=int)
            grp = df.loc[df["@source_node"].isin(pre_gids)]
            lo, hi = pip.get("path_dist", [50., 150.])
            log.update({"pre_gids": pre_gids.tolist(), "epsp": res["epsp"], "epsps": res["epsps"],
                        "spiked": res["spiked"], "n_syn": res["n_syn"], "n_edges": len(grp),
                        "n_basal_in_range": int(((grp["afferent_section_type"] == 2) &
                                                 grp["path_dist"].between(lo, hi)).sum()),
                        "trials": {str(k): v["epsp"] for k, v in trials.items()},
                        "pre_mtype_cells": grp.groupby("mtype", observed=True)["@source_node"].nunique().to_dict(),
                        "pre_mtype_syns": grp.groupby("mtype", observed=True).size().to_dict(),
                        "rho0_missing": grp.loc[grp["rho0_GB"].isna()].groupby("mtype", observed=True)
                                           .size().to_dict()})
            rng_mv = pip.get("epsp_range", [1., 3.])
            log["ok"] = bool(rng_mv[0] <= res["epsp"] <= rng_mv[1] and not res["spiked"])
        with open(logf, "w", encoding="utf-8") as f:
            json.dump(log, f, indent=4, default=float)
        return (tuple(log["pre_gids"]), post_gid) if log["ok"] else None

    def find_pipette_groups(self, post_gids=None):
        """Pipette mode of `find_pairs()`: (tuple of pre gids, post_gid) for `npairs` post cells of `post_mtype`"""
        save_dir = os.path.join(self.out_dir, "single_cells")
        os.makedirs(save_dir, exist_ok=True)
        sim_config = {"run": {"dt": 0.025, "tstop": self.T, "random_seed": self.seed},
                      "network": self.pre_simwriter_circuit_config,
                      "node_sets_file": self.node_set,
                      "node_set": self.target,
                      "output": {"output_dir": os.path.join(self.sims_dir, "out")}}
        sim_config_path = os.path.join(save_dir, "simulation_config.json")
        with open(sim_config_path, "w", encoding="utf-8") as f:
            json.dump(sim_config, f, indent=4)
        all_post = self.pilot_post_gids(None)
        if post_gids is None:  # `post_gids` (e.g. a pilot) keep their index in the shuffled order -> same rng seeds
            post_gids = all_post
        pos = {int(gid): i for i, gid in enumerate(all_post)}
        groups = []
        args_list = [(int(gid), pos[int(gid)], sim_config_path, save_dir) for gid in post_gids]
        with ProcessPoolExecutor(max_workers=min(30, self.npairs)) as pool:
            futures = [pool.submit(self._pipette_post_gid, args) for args in args_list]
            try:
                for future in tqdm(futures, desc="Pipette groups"):
                    result = future.result()
                    if result is not None:
                        groups.append(result)
                        if len(groups) >= self.npairs:
                            break
            finally:
                for future in futures:
                    future.cancel()
        if len(groups) < self.npairs:
            warnings.warn("Not enough pipette groups found")
        return groups

    def pilot_post_gids(self, n):
        """First `n` (all if None) post gids in `find_pipette_groups()` order"""
        nrn = Circuit(self.pre_simwriter_circuit_config).nodes[self.node_pop].get(self.target, ["mtype"])
        post_gids = nrn.index[nrn["mtype"].isin(self.post_mtype)].to_numpy()
        np.random.seed(self.seed)
        np.random.shuffle(post_gids)
        return post_gids[:n]

    def find_pairs(self):
        """Finds connected pairs of gids based on constraints specified in the config"""
        if self.pipette:
            return self.find_pipette_groups()
        save_dir = os.path.join(self.out_dir, "single_cells")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # Write simple simulation config (for checking electrical constraints)
        sim_config = {"run": {"dt": 0.025, "tstop": self.T, "random_seed": self.seed},
                      "network": self.pre_simwriter_circuit_config,
                      "node_sets_file": self.node_set,
                      "node_set": self.target,
                      "output": {"output_dir": os.path.join(self.sims_dir, "out")}}
        sim_config_path = os.path.join(save_dir, "simulation_config.json")
        with open(sim_config_path, "w", encoding="utf-8") as f:
            json.dump(sim_config, f, indent=4)

        c = Circuit(self.pre_simwriter_circuit_config)
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
        account = self.config.get("simulator", {}).get("account", "rrg-emuller")
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

    def load_pairs(self):
        """(pregid, postgid) pairs of the `pairs_from` index csv (in order of first appearance)"""
        df = pd.read_csv(self.pairs_from)
        return list(dict.fromkeys(zip(df["pregid"].astype(int), df["postgid"].astype(int))))

    def calibrate_stimuli(self, pairs, max_workers=30):
        """Searches (in parallel over postsynaptic cells) the current amplitudes for every protocol
        post stimulus missing from the single_cells pkls. Returns {post_gid: [failed stimulus keys]}"""
        save_dir = os.path.join(self.out_dir, "single_cells")
        os.makedirs(save_dir, exist_ok=True)
        sim_config = {"run": {"dt": 0.025, "tstop": self.T, "random_seed": self.seed},
                      "network": self.calibration_circuit_config,
                      "node_sets_file": self.node_set,
                      "node_set": self.target,
                      "output": {"output_dir": os.path.join(self.sims_dir, "out")}}
        sim_config_path = os.path.join(save_dir, "calibration_simulation_config.json")
        with open(sim_config_path, "w", encoding="utf-8") as f:
            json.dump(sim_config, f, indent=4)
        post_gids = sorted({post_gid for _, post_gid in pairs})
        # Presynaptic steps (`pre_type: step`) are calibrated on the pre cells, in the same task as the post
        # stimuli when a gid is both (one writer per pkl). Without such protocols this is the post-only loop.
        pre_protocols = [p for p in self.protocols if protocol_pre_stim_key(p, self.config["stimulus"])]
        pre_gids = sorted({int(g) for pre_gid, _ in pairs for g in np.atleast_1d(pre_gid)}) if pre_protocols else []
        failed = {}
        gids = sorted(set(post_gids) | set(pre_gids))
        with ProcessPoolExecutor(max_workers=min(max_workers, len(gids))) as pool:
            futures = {pool.submit(calibrate_protocol_stimuli, sim_config_path, int(gid),
                                   self.protocols if gid in post_gids else [], self.config["stimulus"], save_dir,
                                   pre_protocols if gid in pre_gids else ()): gid for gid in gids}
            for future in tqdm(futures, desc="Calibrating protocol stimuli"):
                keys = future.result()
                if keys:
                    failed[futures[future]] = keys
        for gid, keys in failed.items():
            warnings.warn("gid %i: no valid stimulus for %s - these protocols are skipped for it" % (gid, keys))
        return failed

    def _sim_config(self, workdir, jsonf_name, tstop, inputs, soma_end, rho_end):
        """SONATA simulation config shared by the full and the prefire simulations"""
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
            "type": "configure_all_sections",
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
        return {"run": {"dt": 0.025, "tstop": tstop, "random_seed": self.seed},
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
                             "start_time": 0.0, "end_time": soma_end},
                    "rho": {"cells": "postcell", "type": "synapse",
                            "variable_name": "GluSynapse.rho_GB", "sections": "all", "unit": "nd", "dt": 0.1,
                            "start_time": 0.0, "end_time": rho_end}
                },
                "target_simulator": "CORENEURON",
                "connection_overrides": connection_overrides}

    def write_sim_files(self, pairs):
        """Writes pair and protocol specific `simulation_config.json` (full: C01 test pulses + pairing + C02)
        and `prefire_simulation_config.json` (pairing only) used by `bluecellulab` / neurodamus,
        and batch scripts to launch single cell sims. Protocols come from `self.protocols` (the classic
        freq x dt grid unless the yaml has a `protocols:` list); workdirs are `<pre>-<post>/<protocol id>`"""
        # Copy config file to sims dir (just to make sure that one can more or less know what was run...)
        basedir = os.path.split(os.path.split(self.out_dir)[0])[0]
        if not os.path.exists(basedir):
            os.makedirs(basedir)
        shutil.copyfile(self._config_path, os.path.join(basedir, "%s.yaml" % self.index_label))
        stim_config = self.config["stimulus"]
        failed = self.calibrate_stimuli(pairs, max_workers=getattr(self, "calibration_workers", 30))
        # Presynaptic test spike times (for C01 and C02) are protocol independent
        n_spikes_before = int(self.C01_duration * MIN2MS / self.C01_T)
        n_spikes_after = int(self.C02_duration * MIN2MS / self.C02_T)
        before_pre_spikes = [self.offset + i * self.C01_T for i in range(n_spikes_before)]
        before_duration = n_spikes_before * self.C01_T
        with open(os.path.join("templates", "simulation.batch.tmpl"), "r") as f:
            templ = f.read()

        all_sims = []
        for proto in self.protocols:
            nreps, T = proto["nreps"], proto["T"]
            after_pre_spikes = [self.offset + n_spikes_before * self.C01_T + nreps * T +
                                i * self.C02_T for i in range(n_spikes_after)]
            pairing_duration = nreps * T
            t_stop = before_duration + pairing_duration + n_spikes_after * self.C02_T
            prefire_t_stop = self.offset + pairing_duration + 1000.0
            # CPU time heuristics: estimated sim time + 2 hour buffer
            h, m = np.divmod(OPT_CPU_TIME * t_stop / 1000. + 7200, 3600)
            m, s = np.divmod(m, 60)
            cpu_time = "%.2i:%.2i:%.2i" % (h, m, s)
            stim_key = protocol_stim_key(proto, stim_config)
            pre_key = protocol_pre_stim_key(proto, stim_config)  # None unless `pre_type: step`
            # post step (`post_type: step`): one "step<i>" pulse input of width step_duration per repetition;
            # the "step" name tells simulator_edges' induction guardrail to count APs per step, not per pulse
            post_step = proto.get("post_type", "pulse") == "step"
            input_name, width = ("step%i", proto["step_duration"]) if post_step else ("pulse%i", proto["width"])
            if post_step and proto.get("step_below_nA"):
                input_name = "substep%i"  # subthreshold step: simulator_edges' guardrail requires NO AP in it
            for pre_gid, post_gid in pairs:
                simres, amplitude, pre_simres = None, 0., None
                if stim_key is not None:
                    if post_gid in failed and stim_key in failed[post_gid]:
                        continue
                    with open(os.path.join(self.out_dir, "single_cells", "%i.pkl" % post_gid), "rb") as f:
                        simres = pickle.load(f)[stim_key]
                    amplitude = simres["amp"]
                if pre_key is not None:
                    assert not isinstance(pre_gid, tuple), "pre steps need a single presynaptic cell"
                    if pre_gid in failed and pre_key in failed[pre_gid]:
                        continue
                    with open(os.path.join(self.out_dir, "single_cells", "%i.pkl" % pre_gid), "rb") as f:
                        pre_simres = pickle.load(f)[pre_key]
                # One pairing repetition (pulses are periodic with period T, pre spikes are repeated explicitly)
                post_onsets, pre_rep = pairing_times(proto, simres, self.offset + before_duration, pre_simres)
                prefire_post_onsets, prefire_pre_rep = pairing_times(proto, simres, self.offset, pre_simres)
                events = np.concatenate([post_onsets - self.offset - before_duration,
                                         post_onsets - self.offset - before_duration + width,
                                         pre_rep - self.offset - before_duration])
                if len(events) and (events.max() - min(events.min(), 0.) >= T or events.min() < -self.offset):
                    raise ValueError("Protocol %s does not fit in its period T=%.0f ms (events %.1f..%.1f ms)"
                                     % (proto["id"], T, events.min(), events.max()))
                # pre_gid is a tuple for pipette groups (see `find_pipette_groups()`): all of them fire together
                pre_gids = [int(g) for g in np.atleast_1d(pre_gid)]
                pair = "pip0-%i" % post_gid if isinstance(pre_gid, tuple) else "%i-%i" % (pre_gid, post_gid)
                workdir = os.path.join(self.out_dir, pair, proto["id"])
                if not os.path.exists(workdir):
                    os.makedirs(workdir)
                # Write node set with idx of pre- and postsynaptic neurons
                jsonf_name = os.path.join(workdir, "node_sets.json")
                node_sets = {"precell": {"node_id": pre_gids, "population": self.node_pop},
                             "postcell": {"node_id": [int(post_gid)], "population": self.node_pop},
                             "paircells": {"node_id": pre_gids + [int(post_gid)], "population": self.node_pop}}
                with open(jsonf_name, "w", encoding="utf-8") as f:
                    json.dump(node_sets, f, indent=4)
                inputs = {input_name % i: {"input_type": "current_clamp", "module": "pulse",
                                          "node_set": "postcell",
                                          "delay": post_spike, "duration": pairing_duration, "amp_start": amplitude,
                                          "width": width, "frequency": 1000. / T}
                          for i, post_spike in enumerate(post_onsets)}
                prefire_inputs = {input_name % i: {"input_type": "current_clamp", "module": "pulse",
                                                  "node_set": "postcell",
                                                  "delay": post_spike, "duration": pairing_duration,
                                                  "amp_start": amplitude,
                                                  "width": width, "frequency": 1000. / T}
                                  for i, post_spike in enumerate(prefire_post_onsets)}
                # Generate (full) presynaptic spike train used as spike replay stimulus
                pre_spikes = [t + j * T for j in range(nreps) for t in pre_rep]
                pre_spikes = np.array(before_pre_spikes + pre_spikes + after_pre_spikes)
                h5f_name = os.path.join(workdir, "prespikes.h5")
                if not self.prefire_only:
                    save_spikes(h5f_name, self.node_pop, np.tile(pre_spikes, len(pre_gids)),
                                np.repeat(pre_gids, len(pre_spikes)))
                prefire_pre_spikes = np.array([t + j * T for j in range(nreps) for t in prefire_pre_rep])
                h5f_prefire_name = os.path.join(workdir, "prefire_prespikes.h5")
                save_spikes(h5f_prefire_name, self.node_pop, np.tile(prefire_pre_spikes, len(pre_gids)),
                            np.repeat(pre_gids, len(prefire_pre_spikes)))
                # Add spike replay to inputs for neurodamus (needed for correct synapse creation).
                # node_set is the TARGET (postsynaptic) cell; source is the population for replay matching.
                inputs["prespikes"] = {"input_type": "spikes", "module": "synapse_replay", "node_set": "postcell",
                                       "delay": 0., "duration": t_stop, "source": self.node_pop,
                                       "spike_file": h5f_name}
                if not self.prefire_only:
                    sim_config = self._sim_config(workdir, jsonf_name, t_stop, inputs, t_stop, t_stop)
                    with open(os.path.join(workdir, "simulation_config.json"), "w", encoding="utf-8") as f:
                        json.dump(sim_config, f, indent=4)
                prefire_inputs["prespikes"] = {"input_type": "spikes", "module": "synapse_replay",
                                               "node_set": "postcell", "delay": 0.,
                                               "duration": prefire_t_stop, "source": self.node_pop,
                                               "spike_file": h5f_prefire_name}
                prefire_sim_config = self._sim_config(workdir, jsonf_name, prefire_t_stop, prefire_inputs,
                                                      prefire_t_stop, t_stop)
                with open(os.path.join(workdir, "prefire_simulation_config.json"), "w", encoding="utf-8") as f:
                    json.dump(prefire_sim_config, f, indent=4)
                # Write launch scripts
                if self.prefire_only:  # index points at the prefire config instead of the batch script
                    f_name = os.path.join(workdir, "prefire_simulation_config.json")
                else:
                    f_name = os.path.join(workdir, "simulation.batch")
                    self.write_batch_sript(f_name, templ, cpu_time)
                self.write_neurodamus_sbatch(workdir, "00:30:00")
                all_sims.append((pair, pre_gids[0] if len(pre_gids) == 1 else -1, post_gid, proto["freq"], proto["dt"], f_name, proto["id"],
                                 proto["n_pre"], proto["n_post"], nreps, T))
        sim_idx = pd.DataFrame(all_sims, columns=["pair", "pregid", "postgid", "frequency", "dt", "path", "protocol_id",
                                                  "n_pre", "n_post", "nreps", "T"])
        sim_idx.to_csv(os.path.join(basedir, "index_%s.csv" % self.index_label), index=False)

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
    parser.add_argument("--config", default="../configs/Sabrina_L5TTPC_L5TTPC_STDP.yaml",
                        help="yaml config (default: %(default)s)")
    parser.add_argument("--npairs", type=int, default=None,
                        help="Override the number of pairs to find (default: read from config)")
    parser.add_argument("--pairs-from", default=None,
                        help="Index csv whose pairs to reuse (overrides the yaml's `pairs_from`), e.g. to rewrite "
                             "workdirs after a protocol change without searching pairs again")
    parser.add_argument("--pilot", type=int, default=0,
                        help="Pipette mode only: run the first N post cells (of the full run's order) only")
    parser.add_argument("--workers", type=int, default=30,
                        help="Processes for the per-cell stimulus calibration in write_sim_files (default: %(default)s)")
    args = parser.parse_args()

    writer = OptSimWriter(args.config)
    writer.calibration_workers = args.workers
    if args.pairs_from:
        writer.config["pairs_from"] = args.pairs_from
    if args.npairs is not None:
        writer.config["npairs"] = args.npairs
    # `pairs_from` in the yaml reuses the pairs of an existing index csv instead of searching new ones
    if args.pilot:  # pipette mode: the first `pilot` post cells of the full run's order, same out_dir
        pairs = writer.find_pipette_groups(post_gids=writer.pilot_post_gids(args.pilot))
    else:
        pairs = writer.load_pairs() if writer.pairs_from else writer.find_pairs()
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

