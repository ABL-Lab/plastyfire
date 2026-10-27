"""Rough counterfactual for the rho0 assignment. Linear EPSP proxy E = sum_i gmax_state_i * Use_state_i (no nrrp, no dendritic
integration). Compares, for L23->L5 (rho0 from the circuit, independent of conductance) the ratio with (a) the actual rho0 and
(b) rho0 re-assigned by the L5->L5 rule (per pair, rho0=1 for the synapses with conductance >= the pair median), holding the
final state rho_f (from the full-rho run) fixed. Synapse gmax_d/gmax_p and Use_d/Use_p are re-derived from the original
conductance c and u with k_gsyn=2, k_u=0.2 (create_dhuruva_edges.py:85-100), so they follow the assigned rho0. NOT a re-simulation."""
import numpy as np, pandas as pd
R = "/home/dhuruva/ws/plastyfire"; O = R + "/ebner/audit_l23l5"
m = pd.read_csv(O + "/syn_table_merged.csv"); d = m[m.set == "L23L5"].copy()
K_G, K_U = 2.0, 0.2
def params(c, u, r0):
    Ud = np.where(r0, u ** (1 / K_U), u); Up = np.where(r0, u, u ** K_U)
    gd = np.where(r0, c / K_G, c); gp = np.where(r0, c, K_G * c)
    return gd, gp, Ud, Up
out = []
for pair, g in d.groupby("pair"):
    c, u = g.conductance.values, g.u_syn.values; rf = g.rho_f.values > .5; r0a = g.rho0.values > .5
    r0b = c >= np.percentile(c, 50)
    for tag, r0 in (("actual", r0a), ("median_split", r0b)):
        gd, gp, Ud, Up = params(c, u, r0)
        E = lambda s: np.sum(np.where(s, gp * Up, gd * Ud))
        out.append(dict(pair=pair, tag=tag, ratio=E(rf) / E(r0), frac_r0=r0.mean()))
o = pd.DataFrame(out)
print(o.groupby("tag").agg(ratio=("ratio", "mean"), sem=("ratio", "sem"), frac_r0=("frac_r0", "mean")).round(3))
# same proxy on the L5L5 Ebner set with its own (median split) rho0 from the edges
e = pd.read_csv(O + "/edges_l5l5.csv").drop_duplicates("syn")
l5 = m[m.set == "L5L5"].merge(e[["syn", "Use_d_TM", "Use_p_TM", "gmax_d_AMPA", "gmax_p_AMPA"]], on="syn", suffixes=("", "_e"))
r = []
for pair, g in l5.groupby("pair"):
    E = lambda s: np.sum(np.where(s, g.gmax_p_AMPA * g.Use_p_TM, g.gmax_d_AMPA * g.Use_d_TM))
    r.append(E(g.rho_f.values > .5) / E(g.rho0.values > .5))
print("L5L5 linear proxy ratio (actual edges rho0):", round(np.mean(r), 3), "+-", round(pd.Series(r).sem(), 3), "n", len(r))
