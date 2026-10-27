"""(a) rho0 -> final-state transition probabilities per set and c_post-floor group (full-rho, delta-cooker, S&H 50 Hz +10).
(b) per-pair rho0 split; (c) mean-field Tsodyks-Markram release of a 5-spike 50 Hz train for the edge parameters of the
two pathways (release per spike relative to c_pre's Use=1 reference), i.e. how much pre-driven Ca per c_pre the release dynamics give."""
import numpy as np, pandas as pd
R = "/home/dhuruva/ws/plastyfire"; O = R + "/ebner/audit_l23l5"
m = pd.read_csv(O + "/syn_table_merged.csv")
m["r0"] = m.rho0 > .5; m["pot"] = m.rho_f > .5
print("(a) transitions  P(0->1), P(1->0), frac potentiated before/after")
for (s, fl), g in m.groupby(["set", "floor"]):
    print(f"  {s:6s} floor={fl!s:5s} n={len(g):4d}  P(0->1)={g[~g.r0].pot.mean():.2f}  P(1->0)={(~g[g.r0].pot).mean():.2f}  pot before {g.r0.mean():.2f} after {g.pot.mean():.2f}")
print("(b) per-pair rho0 fraction (L23L5):", m[m.set == "L23L5"].groupby("pair").r0.mean().describe()[["mean", "50%", "min", "max"]].round(2).to_dict())
d = m[m.set == "L23L5"]
print("    mean conductance rho0=1 vs 0 (nS):", d.groupby("r0").conductance.mean().round(3).to_dict())
# (c) mean-field TM (GluSynapse: u facilitates, R depletes), per edge-parameter sets
def tm(U, D, F, n=5, dt=20.):
    u, Rv, out = 0., 1., []
    for k in range(n):
        if k == 0: u = U
        else:
            u = u * np.exp(-dt / F) ; Rv = 1 - (1 - Rv) * np.exp(-dt / D) if False else Rv
        out.append(None)
    return out
def release(U, D, F, n=5, dt=20.):
    u_s, Rv, rel = U, 1.0, []
    for k in range(n):
        if k > 0:
            u_s = U + (u_s - U) * np.exp(-dt / F)      # facilitation decays to U
            Rv = 1 - (1 - Rv) * np.exp(-dt / D)
        ueff = u_s; r = ueff * Rv; rel.append(r)
        Rv -= r; u_s = u_s + U * (1 - u_s)
    return np.array(rel)
TAU = 278.3
for name, s in (("L23L5", "l23l5"), ("L5L5", "l5l5")):
    e = pd.read_csv(f"{O}/edges_{s}.csv")
    U, D, F = e.u_syn.median(), e.depression_time.median(), e.facilitation_time.median()
    rel = release(U, D, F); t_end = 4 * 20.
    acc = (rel * np.exp(-(t_end - 20. * np.arange(5)) / TAU)).sum()
    print(f"(c) {name}: median U={U:.2f} D={D:.0f} F={F:.1f} -> release/spike {np.round(rel,2)} sum {rel.sum():.2f}; decay-weighted (tau {TAU}) {acc:.2f} x c_pre-equivalents")
    # per-synapse
    rr = np.array([release(u, dd, max(ff, 1e-3)) for u, dd, ff in zip(e.u_syn, e.depression_time, e.facilitation_time)])
    acc = (rr * np.exp(-(80. - 20. * np.arange(5)) / TAU)[None, :]).sum(1)
    print(f"    per-synapse decay-weighted sum: mean {acc.mean():.2f} median {np.median(acc):.2f}; first spike {rr[:,0].mean():.2f}")
