"""delta-bd: og-delta with distance-dependent basal Ka and LVA (runtime version of the new hoc's basal block).

Basal segment at path distance d (um, from soma(0.5), as the hoc's distribute_distance):
    gNaTgbar_NaTg          = na                                   (uniform)
    gbar_Ka_kampa          = ka0 * (1 + d / lam)                  (lam = inf -> uniform)
    gCa_LVAstbar_Ca_LVAst  = OG_LVA * (1 + m * sigmoid((d - dh) / 15))
HVA, BK and every non-basal section keep the og-delta values of the circuit's own hoc.
"""
import numpy as np

# og-delta basal block (SSCx-AAD-delta-emodels/cADpyr_L5TPC.hoc; og_delta_values.json)
OG = dict(na=0.003, ka0=0.002, lam=np.inf, m=0.0, dh=150.0)
OG_LVA = 0.0023897659621395684
ANTIC = dict(OG, na=0.0143, ka0=0.025)                       # antic-delta: basal Na, Ka only
BOUNDS = dict(na=(0.0, 0.03), ka0=(0.0, 0.05), lam=(20.0, 1e4), m=(0.0, 8.0), dh=(60.0, 250.0))
SIG_W = 15.0


def densities(P, d):
    """Basal densities (S/cm2) at path distance d (array ok)."""
    P = {**OG, **P}
    d = np.asarray(d, float)
    ka = P["ka0"] * (1.0 + d / P["lam"]) if np.isfinite(P["lam"]) else np.full_like(d, P["ka0"])
    lva = OG_LVA * (1.0 + P["m"] / (1.0 + np.exp(-(d - P["dh"]) / SIG_W)))
    return dict(gNaTgbar_NaTg=np.full_like(d, P["na"]), gbar_Ka_kampa=ka, gCa_LVAstbar_Ca_LVAst=lva)


def basal_segments(cell, h):
    """[(seg, distance)] for every basal segment of a bluecellulab cell."""
    out = []
    for sec in cell.basal:
        for seg in sec:
            out.append((seg, h.distance(cell.soma(0.5), seg)))
    return out


def apply(segs, P, ka_off=False, lva_off=False):
    """Write delta-bd densities into the basal segments; ka_off / lva_off = in-silico 4-AP / Ni2+ (basal only)."""
    d = np.array([x for _, x in segs])
    g = densities(P, d)
    if ka_off:
        g["gbar_Ka_kampa"] = np.zeros_like(d)
    if lva_off:
        g["gCa_LVAstbar_Ca_LVAst"] = np.zeros_like(d)
    for i, (seg, _) in enumerate(segs):
        for k, v in g.items():
            setattr(seg, k, float(v[i]))


def hoc_lines(P):
    """distribute_distance lines for the new hoc's biophys() (after the basal forsec block)."""
    P = {**OG, **P}
    ka = (f"{P['ka0']!r}*(1+(%.17g)/{P['lam']!r})" if np.isfinite(P["lam"]) else f"{P['ka0']!r}+0*(%.17g)")
    lva = f"{OG_LVA!r}*(1+{P['m']!r}/(1+exp(-((%.17g)-{P['dh']!r})/{SIG_W!r})))"
    return [f'  distribute_distance(CellRef.basal, "gNaTgbar_NaTg", "{P["na"]!r}+0*(%.17g)")',
            f'  distribute_distance(CellRef.basal, "gbar_Ka_kampa", "{ka}")',
            f'  distribute_distance(CellRef.basal, "gCa_LVAstbar_Ca_LVAst", "{lva}")']
