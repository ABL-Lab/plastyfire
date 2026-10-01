"""K for t_drive 4 = 0.2% of the median per-bAP peak of -ica_VDCC, measured on the Sjostrom 2007 post-only prefire
traces (post APs without pre spikes: the pure bAP influx) for the ljp-0 and ljp-25 hashes, 24 subset pairs.
Run in Slurm:  python glusynapse_v2/k_ca_measure.py"""
import os, pickle, numpy as np
ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
SIMS = f"{ROOT}/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations"
PAIRS = open(f"{ROOT}/glusynapse_v2/subset24_pairs.txt").read().strip().split(",")
for h in ("delta-prefire-vseg", "delta-ljp25-prefire-vseg"):
    syn_med, pooled = [], []
    for p in PAIRS:
        d = f"{SIMS}/{p}/sjostrom07_step200ms_post_only"
        tf, rf = f"{d}/bluecellulab_results_{h}/simulation_traces.pkl", f"{d}/simulation_edges_{h}.pkl"
        if not (os.path.isfile(tf) and os.path.isfile(rf)):
            continue
        tr = pickle.load(open(tf, "rb")); post = np.asarray(pickle.load(open(rf, "rb"))["postspikes"], float)
        t0 = float(tr["t"][0]); dt = 0.025
        for sid, ica in tr["ica_VDCC"].items():
            q = -np.asarray(ica)
            pk = [q[max(int(round((tp - t0) / dt)) - 40, 0):int(round((tp - t0) / dt)) + 400].max() for tp in post]
            syn_med.append(np.median(pk)); pooled += pk
        del tr
    syn_med = np.array(syn_med); pooled = np.array(pooled)
    print(f"{h}: {len(syn_med)} synapses, {len(pooled)} bAPs; median of per-synapse median peak {np.median(syn_med):.4e} nA "
          f"-> K(0.2%) {0.002 * np.median(syn_med):.4e}; pooled median {np.median(pooled):.4e} -> K {0.002 * np.median(pooled):.4e}", flush=True)
