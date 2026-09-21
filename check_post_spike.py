#!/usr/bin/env python3
"""
Does the post cell actually fire under the c_post stimulus?

c_post is calcium from a BACKPROPAGATING AP. No somatic spike means no bAP, no VDCC
influx, and effcai_GB stays at its resting floor (~3e-5 instead of ~0.07). That is a
silent failure: _find_cpre_cpost reports the floor as a c_post value, so a whole cache
can be built from cells that never spiked.

_find_cpre_cpost scrapes the stimulus amplitude from the workdir config, which was
calibrated against whichever emodels were current then. Modified ion channels change
rheobase, so the same current may no longer reach threshold.

This replicates that exact stimulus (IClamp, delay=1000, dur=3.0, amp=amp_start),
reports spikes and peak Vm, and — with --sweep — finds the amplitude that does fire
exactly one AP.

Usage:
    python check_post_spike.py --pairs 180164-197248 180236-203656 \
        --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json
    python check_post_spike.py --pairs 180164-197248 --sweep \
        --circuit-config data/dhuruva_modified_ion_channels_circuit_config.json
"""

import argparse
import json
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BASE = os.path.join(
    HERE, "refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations")
NODE_POP = "S1nonbarrel_neurons"
PROTOS = ["10Hz_-10ms", "10Hz_10ms", "10Hz_5ms", "10Hz_-30ms",
          "10Hz_30ms", "10Hz_-50ms", "10Hz_50ms"]


def find_workdir(pair):
    for proto in PROTOS:
        wd = os.path.join(BASE, pair, proto)
        if os.path.isfile(os.path.join(wd, "prefire_simulation_config.json")):
            return wd
    return None


def run_one(workdir, circuit_config, amp=None, fixhp=True):
    """Fire the c_post stimulus once; return (n_spikes, peak_mV, amp_used)."""
    import bluecellulab
    from bluepysnap import Simulation
    from plastyfire.simulator import _get_spikes
    import plastyfire.simulator_edges as se

    bluecellulab.neuron.load_mechanisms(se.MECHANISMS_PATH)

    with open(os.path.join(workdir, "prefire_simulation_config.json")) as f:
        cfg = json.load(f)
    ns = os.path.join(workdir, "node_sets.json")
    if os.path.exists(ns):
        cfg["node_sets_file"] = os.path.abspath(ns)
    if circuit_config:
        cfg["network"] = os.path.abspath(circuit_config)
    cfg.pop("reports", None)
    cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    json.dump(cfg, tmp)
    tmp.close()
    try:
        snap = Simulation(tmp.name)
        post_gid = int(snap.node_sets.content["postcell"]["node_id"][0])

        if amp is None:  # same scrape _find_cpre_cpost does
            amp = 1.0
            for stim in cfg.get("inputs", cfg.get("stimuli", {})).values():
                if stim.get("input_type") == "current_clamp" and "amp_start" in stim:
                    amp = float(stim["amp_start"])
                    break

        sim = bluecellulab.CircuitSimulation(
            tmp.name, base_seed=cfg["run"]["random_seed"])
        sim.instantiate_gids([(NODE_POP, post_gid)], add_synapses=False,
                             add_minis=False, add_pulse_stimuli=False)
        cell = sim.cells[(NODE_POP, post_gid)]
        if fixhp:
            for sec in cell.somatic + cell.axonal:
                sec.uninsert("SK_E2")

        stim = bluecellulab.neuron.h.IClamp(0.5, sec=cell.soma)
        stim.delay, stim.dur, stim.amp = 1000.0, 3.0, amp
        cell.persistent.append(stim)

        sim.run(1500.0, cvode=True)
        t = np.array(sim.get_time())
        v = np.array(sim.get_voltage_trace((NODE_POP, post_gid)))
        return len(_get_spikes(t, v)), float(v.max()), amp, post_gid
    finally:
        os.unlink(tmp.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", nargs="+", required=True)
    ap.add_argument("--circuit-config", default=None)
    ap.add_argument("--sweep", action="store_true",
                    help="scan amplitudes to find one that fires exactly 1 AP")
    ap.add_argument("--sweep-max", type=float, default=6.0)
    args = ap.parse_args()

    print(f"config: {args.circuit_config or '(workdir default)'}")
    print("=" * 72)
    for pair in args.pairs:
        wd = find_workdir(pair)
        if wd is None:
            print(f"{pair}: no workdir")
            continue
        n, vmax, amp, gid = run_one(wd, args.circuit_config)
        verdict = "SPIKES" if n else "NO SPIKE — c_post invalid"
        print(f"\n{pair}  (post {gid})")
        print(f"  config amp = {amp:.4f} nA -> {n} spike(s), peak Vm {vmax:.1f} mV"
              f"   [{verdict}]")

        if args.sweep and n == 0:
            print("  sweeping for threshold …")
            a = amp
            while a <= args.sweep_max:
                a = round(a + 0.25, 4)
                sn, sv, _, _ = run_one(wd, args.circuit_config, amp=a)
                print(f"    {a:>5.2f} nA -> {sn} spike(s), peak {sv:.1f} mV", flush=True)
                if sn >= 1:
                    print(f"  => fires at {a:.2f} nA "
                          f"({a / amp:.2f}x the config value)")
                    break
            else:
                print(f"  => no spike up to {args.sweep_max} nA")
    return 0


if __name__ == "__main__":
    sys.exit(main())
