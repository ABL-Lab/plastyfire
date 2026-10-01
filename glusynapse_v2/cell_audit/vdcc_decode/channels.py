"""Gating curves of the Ca channels as compiled in DEES_cell_packages/x86_64 (celsius 34), vs their sources.
Login node, numpy only.  python cell_audit/vdcc_decode/channels.py  -> figs/channels.png, channels.json"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
CEL = 34.0
V = np.linspace(-100, 40, 1401)


def hva(v, hshift=12.3, q10=True):
    """Ca_HVA2 (other_mods/modified_mechanisms, compiled): Reuveni 1993 rates, h at v+12.3, qt = 2.3^((34-36)/10)."""
    v = np.where(v == -27, v + 1e-4, v)
    a = 0.055 * (-27 - v) / (np.exp((-27 - v) / 3.8) - 1); b = 0.94 * np.exp((-75 - v) / 17)
    qt = 2.3 ** ((CEL - 36) / 10) if q10 else 1.0
    vh = v + hshift
    ha = 0.000457 * np.exp((-13 - vh) / 50); hb = 0.0065 / (np.exp((-vh - 15) / 28) + 1)
    return a / (a + b), 1 / (a + b) / qt, ha / (ha + hb), 1 / (ha + hb) / qt


def lva(v, shift=10.0):
    """Ca_LVAst (Hay 2011 file, compiled): curves at v+10 ('-10 mV LJP'), qt = 2.3^((34-21)/10)."""
    qt = 2.3 ** ((34 - 21) / 10); vv = v + shift
    return (1 / (1 + np.exp((vv + 30) / -6)), (5 + 20 / (1 + np.exp((vv + 25) / 5))) / qt,
            1 / (1 + np.exp((vv + 80) / 6.4)), (20 + 50 / (1 + np.exp((vv + 40) / 7))) / qt)


def rtype(v, ljp=0.0):
    """GluSynapse spine VDCC (Chindemi 2022 from Magee & Johnston 1995): Boltzmann, fixed taus 1 / 27 ms, no q10."""
    one = np.ones_like(v)
    return (1 / (1 + np.exp(((-5.9 - ljp) - v) / 9.5)), one, 1 / (1 + np.exp(((-39 - ljp) - v) / -9.2)), 27 * one)


def half(v, y):                       # V where y crosses 0.5 (y monotone)
    i = np.where(np.diff(np.sign(y - 0.5)))[0]
    return float(v[i[0]]) if len(i) else np.nan


def main():
    chans = {"Ca_HVA2 (as used)": hva(V), "Ca_HVA2 source (Reuveni 1993 = Hay 2011 Ca_HVA)": hva(V, 0.0, False),
             "Ca_LVAst (as used)": lva(V), "Ca_LVAst source (Avery & Johnston 1996, before -10 mV)": lva(V, 0.0),
             "spine VDCC (as used, ljp 0)": rtype(V), "spine VDCC ljp 5 (s5)": rtype(V, 5.0),
             "spine VDCC, 110 mM Ba -> 2 mM Ca (-25 mV)": rtype(V, 25.0)}
    out = {}
    rest = -78.0
    for k, (mi, mt, hi, ht) in chans.items():
        at = lambda y, x: float(np.interp(x, V, y))
        out[k] = dict(V50_m=half(V, mi), V50_m2=half(V, mi ** 2), V50_h=half(V, hi), h_rest=at(hi, rest),
                      m2_at={x: at(mi ** 2, x) for x in (-70, -60, -50, -44, -27, -10, 0)},
                      mtau_at_m30=at(mt, -30.0), htau_at_m30=at(ht, -30.0))
    json.dump(out, open(os.path.join(HERE, "channels.json"), "w"), indent=1)
    for k, d in out.items():
        print(f"{k:55s} V50 m {d['V50_m']:6.1f} m2 {d['V50_m2']:6.1f} h {d['V50_h']:6.1f}  h(-78) {d['h_rest']:.2f} "
              f"mtau(-30) {d['mtau_at_m30']:.2f} htau(-30) {d['htau_at_m30']:.0f}  m2@-60/-44/-27/0 "
              + " ".join(f"{d['m2_at'][x]:.2g}" for x in (-60, -44, -27, 0)))
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    col = {"HVA": "C0", "LVA": "C1", "spine": "C3"}
    for k, (mi, mt, hi, ht) in chans.items():
        c = col["HVA" if "HVA" in k else "LVA" if "LVA" in k else "spine"]
        ls = "-" if "as used" in k else (":" if "source" in k else "--")
        lw = 1.0 if "ljp 5" in k else 1.8
        ax[0, 0].plot(V, mi ** 2, c=c, ls=ls, lw=lw, label=k); ax[0, 1].plot(V, hi, c=c, ls=ls, lw=lw)
        ax[1, 0].plot(V, mt, c=c, ls=ls, lw=lw); ax[1, 1].semilogy(V, mi ** 2 * np.interp(rest, V, hi), c=c, ls=ls, lw=lw)
    for a in ax.flat:
        a.set_xlabel("V (mV)"); a.axvline(rest, c="k", lw=0.5)
        for x, lab in ((-44, ">150 um bAP pk"), (-27, "60-150"), (0.1, "<60")):
            a.axvline(x, c="grey", lw=0.5, ls="--")
    ax[0, 0].set_ylabel("m_inf^2"); ax[0, 1].set_ylabel("h_inf"); ax[1, 0].set_ylabel("tau_m (ms, 34 C)")
    ax[1, 1].set_ylabel("m_inf^2 * h_inf(rest -78)  (open fraction for a bAP from rest)"); ax[1, 1].set_ylim(1e-5, 1.5)
    ax[1, 0].set_ylim(0, 8)
    ax[0, 0].legend(fontsize=7, loc="upper left")
    ax[0, 0].set_title("grey dashed: median bAP peak V at >150 / 60-150 / <60 um (NMDA_VDCC_CALIBRATION)", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "channels.png"), dpi=110)


if __name__ == "__main__":
    main()
