"""Does the post cell fire one AP per somatic pulse in the induction protocols? (all pairs, first repetition)

For each pair and each protocol whose config has current-clamp pulses: the post cell alone (no synapses),
the protocol's own pulse inputs (bluecellulab builds them from the config), run the first repetition window,
count somatic spikes vs pulses delivered in that window.

    python glusynapse_v2/count_post_spikes.py --pairs 186261-208490 --out glusynapse_v2/results/post_spikes_x.csv
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

SIMS = ("/project/rrg-emuller/dhuruva/plastyfire/refitting_results/fitting/n120/seed20262009/"
        "Sabrina_L5TTPC_L5TTPC_STDP/simulations")
MECH = "/project/rrg-emuller/dhuruva/DEES_cell_packages/"
WIN = 1200.0


def one(cfgp, post):
    import bluecellulab
    from bluecellulab import CircuitSimulation
    cfg = json.load(open(cfgp))
    pulses = {k: v for k, v in cfg.get("inputs", {}).items() if v.get("input_type") == "current_clamp"}
    if not pulses:
        return None
    t0 = min(v["delay"] for v in pulses.values())
    times = []
    for v in pulses.values():
        f = v.get("frequency", 0) or 0
        n = int(v["duration"] * f / 1000.0) + 1 if f else 1
        times += [v["delay"] + k * 1000.0 / f if f else v["delay"] for k in range(n)]
    times = sorted(x for x in times if t0 <= x < t0 + WIN)
    cfg["inputs"] = pulses; cfg.pop("reports", None); cfg["run"]["tstop"] = t0 + WIN
    tmp = os.path.join(os.environ.get("SLURM_TMPDIR", "/tmp"), f"cps_{os.getpid()}.json")
    json.dump(cfg, open(tmp, "w"))
    sim = CircuitSimulation(tmp)
    pop = "S1nonbarrel_neurons"
    sim.instantiate_gids([(pop, post)], add_synapses=False, add_minis=False, add_pulse_stimuli=True)
    sim.run(t0 + WIN, cvode=False, dt=0.025)
    t = np.array(sim.get_time()); v = np.array(sim.get_voltage_trace((pop, post)))
    up = t[1:][(v[1:] > -10) & (v[:-1] <= -10)]
    up = up[up >= t0 - 1]
    w = [p["width"] for p in pulses.values()][0]; a = [p["amp_start"] for p in pulses.values()][0]
    # per pulse: a spike within [pulse start, next pulse start)
    hit = []
    for i, x in enumerate(times):
        x1 = times[i + 1] if i + 1 < len(times) else x + 60
        hit.append(int(np.any((up >= x) & (up < x1))))
    return dict(n_pulses=len(times), n_spikes=len(up), pulses_with_spike=int(sum(hit)),
                missed=[round(times[i] - t0, 2) for i, hh in enumerate(hit) if not hh],
                spikes=[round(s - t0, 2) for s in up], width=w, amp=a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="comma list or 'all'")
    ap.add_argument("--protos", default="", help="substring filter, comma list (default: all with pulses)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import bluecellulab
    bluecellulab.neuron.load_mechanisms(MECH)
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    pairs = sorted(os.listdir(SIMS)) if args.pairs == "all" else args.pairs.split(",")
    rows = []
    for pair in pairs:
        post = int(pair.split("-")[1])
        for cfgp in sorted(glob.glob(os.path.join(SIMS, pair, "*", "prefire_simulation_config.json"))):
            proto = os.path.basename(os.path.dirname(cfgp))
            if args.protos and not any(s in proto for s in args.protos.split(",")):
                continue
            with ctx.Pool(1, maxtasksperchild=1) as pool:
                try:
                    r = pool.apply(one, (cfgp, post))
                except Exception as e:                     # keep going; record the failure
                    r = dict(error=repr(e)[:200])
            if r is None:
                continue
            rows.append(dict(pair=pair, proto=proto, **r))
            print(pair, proto, {k: r.get(k) for k in ("n_pulses", "pulses_with_spike", "n_spikes", "missed", "error")},
                  flush=True)
            pd.DataFrame(rows).to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
