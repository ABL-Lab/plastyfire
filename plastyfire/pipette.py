"""
Extracellular (pipette) stimulation groups, for protocols recorded with a stimulation
electrode instead of a paired recording (Nevian & Sakmann 2006: pipette near the basal
dendrites of an L2/3 PC, 50-150 um from the soma, compound EPSP of 1-3 mV, bicuculline).

A pipette group is a set of excitatory presynaptic cells that all fire the same spike
train. The pipette site is a basal EXC synapse of the post cell at 50-150 um path
distance; presynaptic cells are recruited in order of their nearest basal synapse to the
site (the stimulus excites their axons, so ALL their synapses on the post cell are
active, not just the basal ones near the site). Cells are added until the somatic
compound EPSP reaches the target amplitude. Inhibitory inputs are never recruited
(bicuculline in the experiments).
"""

import os
import json
import tempfile
import multiprocessing
import numpy as np
import pandas as pd

SONATA_BASAL = 2  # afferent_section_type + 1 == NeuriteType (see thresholdfinder.BRANCH_TYPE_OFFSET)
SONATA_APICAL = 3


def syn_path_distances(morph, section_ids, section_pos):
    """Path distance (um) from the soma of synapses given by SONATA (afferent_section_id,
    afferent_section_pos). SONATA section 0 is the soma, section k is morphio section k-1."""
    start = {}
    for sec in morph.iter():  # depth first: parents come before children
        start[sec.id] = 0. if sec.is_root else start[sec.parent.id] + _sec_len(sec.parent)
    out = np.zeros(len(section_ids))
    for i, (sid, pos) in enumerate(zip(section_ids, section_pos)):
        if sid > 0:
            sec = morph.section(int(sid) - 1)
            out[i] = start[sec.id] + pos * _sec_len(sec)
    return out


def _sec_len(sec):
    return float(np.linalg.norm(np.diff(sec.points, axis=0), axis=1).sum())


def afferent_table(circuit, node_pop, edge_pop, post_gid):
    """All chemical afferents of `post_gid` with presynaptic mtype / synapse class and path distance"""
    edges = circuit.edges[edge_pop]
    nodes = circuit.nodes[node_pop]
    props = ["@source_node", "afferent_section_id", "afferent_section_pos", "afferent_section_type",
             "afferent_center_x", "afferent_center_y", "afferent_center_z", "rho0_GB"]
    df = edges.afferent_edges(int(post_gid), props)
    pre = nodes.get(df["@source_node"].unique(), ["mtype", "synapse_class"])
    df = df.join(pre, on="@source_node")
    morph = nodes.morph.get(int(post_gid), transform=False, extension="asc")
    df["path_dist"] = syn_path_distances(morph, df["afferent_section_id"].to_numpy(),
                                         df["afferent_section_pos"].to_numpy())
    return df


def rank_pre_gids(df, rng, path_range):
    """Picks the pipette site (a basal EXC synapse at `path_range` um path distance) and ranks
    EXC presynaptic cells by the distance of their nearest basal synapse to it"""
    exc = df.loc[df["synapse_class"] == "EXC"]
    basal = exc.loc[exc["afferent_section_type"] == SONATA_BASAL]
    cand = basal.loc[basal["path_dist"].between(*path_range)]
    if not len(cand):
        return None, None
    site_edge = cand.index[rng.integers(len(cand))]
    xyz = ["afferent_center_x", "afferent_center_y", "afferent_center_z"]
    site = basal.loc[site_edge, xyz].to_numpy(dtype=float)
    d = np.linalg.norm(basal[xyz].to_numpy(dtype=float) - site, axis=1)
    nearest = pd.Series(d, index=basal.index).groupby(basal["@source_node"]).min().sort_values()
    return {"edge_id": int(site_edge), "xyz": site.tolist(),
            "path_dist": float(basal.loc[site_edge, "path_dist"])}, nearest


def _epsp_process(conn, sim_config, post_gid, pre_gids, node_pop, spike_times, mech_path, fixhp):
    """Subprocess: somatic compound EPSP (mean peak over `spike_times`) of `pre_gids` firing together"""
    try:
        import bluecellulab
        bluecellulab.neuron.load_mechanisms(mech_path)
        with open(sim_config) as f:
            cfg = json.load(f)
        glusyn = dict(cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", {}))
        cfg.pop("reports", None)
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(cfg, tmp)
        tmp.close()
        sim = bluecellulab.CircuitSimulation(tmp.name, base_seed=cfg["run"]["random_seed"])
        sim.instantiate_gids([(node_pop, post_gid)], add_synapses=True, add_minis=False,
                             add_pulse_stimuli=False,
                             intersect_pre_gids=[(node_pop, int(g)) for g in pre_gids],
                             pre_spike_trains={(node_pop, int(g)): list(spike_times) for g in pre_gids})
        cell = sim.cells[(node_pop, post_gid)]
        if fixhp:
            for sec in cell.somatic + cell.axonal:
                sec.uninsert("SK_E2")
        for p in ["cao_CR", "init_depleted", "minis_single_vesicle"]:  # the ones that change the EPSP
            if p in glusyn:
                setattr(bluecellulab.neuron.h, "%s_GluSynapse" % p, float(glusyn[p]))
        # Stochastic release seeded as in simulator_edges (the Neurodamus Random123 pattern): unseeded,
        # every spike gives the identical EPSP
        synapse_seed = int(cfg.get("run", {}).get("synapse_seed") or 0)
        for syn_id, synapse in cell.synapses.items():
            synapse.hsynapse.setRNG(post_gid, 100000 + int(syn_id[1]), 200 + synapse_seed)
        n_syn = len(cell.synapses)
        sim.run(spike_times[-1] + 200., cvode=False)  # fixed step: CVode leaves no samples in the baseline window
        t, v = np.array(sim.get_time()), np.array(sim.get_voltage_trace((node_pop, post_gid)))
        os.unlink(tmp.name)
        peaks = []
        for ts in spike_times:
            base = v[(t >= ts - 20.) & (t < ts)].mean()
            peaks.append(v[(t >= ts) & (t < ts + 60.)].max() - base)
        peak_avg, rise = _avg_waveform(t, v, spike_times)
        conn.send({"epsp": float(np.mean(peaks)), "epsps": peaks, "n_syn": n_syn,
                   "spiked": bool(v.max() > 0.), "peak_avg": peak_avg, "rise_10_90": rise})
    except Exception as e:
        import traceback
        conn.send({"error": "%r\n%s" % (e, traceback.format_exc())})


def _avg_waveform(t, v, spike_times, dt=0.025):
    """Peak (mV) and 10-90% rise time (ms, as Letzkus 2006 Methods) of the spike-averaged,
    baseline-subtracted EPSP (failures included, as in an averaged recording)"""
    grid = np.arange(-20., 60., dt)
    w = np.mean([np.interp(ts + grid, t, v) - v[(t >= ts - 20.) & (t < ts)].mean() for ts in spike_times], axis=0)
    w, g = w[grid >= 0], grid[grid >= 0]
    ip = int(np.argmax(w))
    peak = float(w[ip])
    if peak <= 0:
        return peak, float("nan")

    def cross(frac):  # first upward crossing of frac * peak before the peak, linearly interpolated
        k = int(np.argmax(w[:ip + 1] >= frac * peak))
        if k == 0:
            return g[0]
        return g[k - 1] + (frac * peak - w[k - 1]) / (w[k] - w[k - 1]) * dt

    return peak, float(cross(.9) - cross(.1))


def measure_epsp(sim_config, post_gid, pre_gids, node_pop, spike_times, mech_path, fixhp=True):
    parent, child = multiprocessing.Pipe(duplex=False)
    proc = multiprocessing.Process(target=_epsp_process, args=(child, sim_config, int(post_gid),
                                                                list(pre_gids), node_pop, spike_times,
                                                                mech_path, fixhp))
    proc.start()
    child.close()
    res = parent.recv()
    proc.join()
    if "error" in res:
        raise RuntimeError(res["error"])
    return res


def calibrate_group(ranked, epsp_fn, target, epsp_range, max_pre, max_skip_tests=20):
    """Pre gids (nearest first) whose compound EPSP is closest to `target` (mV). First the prefix of
    `ranked`: doubling until the EPSP passes the target, then bisection. If the best prefix is still out
    of `epsp_range` because one strong cell jumps over it, continue from the last prefix below range in
    `ranked` order, skipping every cell that would push the EPSP above the range (at most
    `max_skip_tests` extra sims). Returns (pre gids, result, trials); trials are keyed by group size,
    with "<size>+s<k>" for the skip search."""
    trials = {}

    def f(n):
        if n not in trials:
            trials[n] = epsp_fn(ranked[:n])
        return trials[n]["epsp"]

    n_max = min(max_pre, len(ranked))
    hi = 1
    while hi < n_max and f(hi) < target:
        hi = min(2 * hi, n_max)
    if f(hi) < target:  # even all recruitable cells stay below target
        return list(ranked[:hi]), trials[hi], trials
    lo = hi // 2 if hi > 1 else 0
    while hi - lo > 1:  # invariant: f(hi) >= target, f(lo) < target (or lo == 0)
        mid = (lo + hi) // 2
        if f(mid) >= target:
            hi = mid
        else:
            lo = mid
    n = hi if lo == 0 or abs(f(hi) - target) <= abs(f(lo) - target) else lo
    def in_range(m):
        return epsp_range[0] <= trials[m]["epsp"] <= epsp_range[1] and not trials[m]["spiked"]

    if not in_range(n) and in_range(hi):
        n = hi
    if in_range(n) or trials[hi]["epsp"] <= epsp_range[1]:
        return list(ranked[:n]), trials[n], trials
    # f(lo) < target and f(hi) > range max: cell `hi` alone jumps over the range -> skip search
    group = list(ranked[:lo])
    best = (list(group), trials[lo]) if lo else (None, None)
    for k, g in enumerate(ranked[hi:n_max][:max_skip_tests]):
        res = epsp_fn(group + [g])
        trials["%i+s%i" % (len(group) + 1, k)] = res
        if res["epsp"] > epsp_range[1] or res["spiked"]:
            continue  # this cell overshoots on its own: skip it
        group.append(g)
        best = (list(group), res)
        if res["epsp"] >= epsp_range[0]:
            break
    if best[0] is None:  # nothing below the range: report the closest prefix (rejected by the caller)
        return list(ranked[:n]), trials[n], trials
    return best[0], best[1], trials
