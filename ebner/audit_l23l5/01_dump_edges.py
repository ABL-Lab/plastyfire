"""Dump all edge properties (shared edges h5, read-only) for the 794 L23->L5 synapses and the Sabrina L5->L5 synapses
(cache keys of sabrina_n120_delta.pkl) to edges_l23l5.csv / edges_l5l5.csv, with pre/post mtype and path distance."""
import os, pickle, numpy as np, pandas as pd, h5py
R = "/project/rrg-emuller/dhuruva/plastyfire"; O = os.path.join(R, "ebner/audit_l23l5")
E = os.path.join(R, "data/dhuruva_modified_edges.h5"); P = "edges/S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/0"
for tag, pk in (("l23l5", "ebner_l23l5_delta"), ("l5l5", "sabrina_n120_delta")):
    d = pickle.load(open(f"{R}/cpre_cpost_cache/{pk}.pkl", "rb"))
    rows = []
    for pair, v in d.items():
        for s in v["c_pre"]: rows.append((pair[0], pair[1], s, v["c_pre"][s], v["c_post"][s]))
    df = pd.DataFrame(rows, columns=["pre", "post", "syn", "c_pre", "c_post"]).sort_values("syn").reset_index(drop=True)
    idx = df.syn.values
    with h5py.File(E, "r") as f:
        g = f[P]
        for k in ["Use_d_TM", "Use_p_TM", "gmax_d_AMPA", "gmax_p_AMPA", "rho0_GB", "conductance", "conductance_scale_factor",
                  "u_syn", "n_rrp_vesicles", "volume_CR", "decay_time", "depression_time", "facilitation_time", "u_hill_coefficient",
                  "theta_d", "theta_p", "afferent_section_type", "spine_length", "syn_type_id", "afferent_center_y", "delay"]:
            df[k] = g[k][idx]
    df.to_csv(f"{O}/edges_{tag}.csv", index=False); print(tag, df.shape)
