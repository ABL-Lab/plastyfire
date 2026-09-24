"""Measure the EFFECTIVE tau_effca: which value reproduces NEURON's recorded rho?

Speed note: effcai MUST be built at the trace's native 0.025 ms (lfilter, C
speed), but the rho ODE does not — validate.py showed 2 ms and 0.25 ms give
identical rho. So decimate effcai after filtering, then step rho vectorised over
all synapses at once. That turns ~280M python iterations into ~21k.

    python fit_tau.py [n_files] [--workers N]
"""
import argparse, glob, pickle
import numpy as np
from scipy.signal import lfilter
from concurrent.futures import ProcessPoolExecutor
from constants import (MIN_CA_CR, RHO_STAR_GB, TAU_IND_GB,
                       GAMMA_D_GB, GAMMA_P_GB, TIED_IC4)

ROOT = "/lustre06/project/6077694/dhuruva/plastyfire"
SIMS = f"{ROOT}/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
RHO_DT_MS = 2.0
TAUS = [100, 150, 200, 225, 250, 278.3177658387, 300, 350, 400, 500]


def load(nf):
    cache = pickle.load(open(f"{ROOT}/cpre_cpost_cache/ion_channels_tau278.pkl", "rb"))
    cp, cq = {}, {}
    for v in cache.values():
        for s, x in v["c_pre"].items():  cp[int(s)] = x
        for s, x in v["c_post"].items(): cq[int(s)] = x
    U, td, tp, r0, rf = [], [], [], [], []
    dt = None
    for f in sorted(glob.glob(f"{SIMS}/*/10Hz_10ms/*/simulation_traces.pkl"))[:nf]:
        tr = pickle.load(open(f, "rb"))
        t = tr["t"]; dt = float(np.median(np.diff(t[:1000])))
        for sid, rr in tr["rho_GB"].items():
            sid = int(sid)
            if sid not in cp: continue
            rr = np.asarray(rr, dtype=np.float64)
            U.append(np.asarray(tr["cai_CR"][sid], dtype=np.float64) - MIN_CA_CR)
            td.append(TIED_IC4["a00"]*cp[sid] + TIED_IC4["a01"]*cq[sid])
            tp.append(TIED_IC4["a10"]*cp[sid] + TIED_IC4["a11"]*cq[sid])
            r0.append(rr[0]); rf.append(rr[-1])
    return np.array(U), dt, np.array(td), np.array(tp), np.array(r0), np.array(rf)


def score(args):
    tau, U, dt, td, tp, r0, obs = args
    a = np.exp(-dt/tau); b = tau*(1-a)
    step = max(1, int(round(RHO_DT_MS/dt)))
    eff = lfilter([0., b], [1., -a], U, axis=1)[:, ::step]     # (n_syn, T)
    h = (dt*step)/(1e3*TAU_IND_GB)
    rho = r0.astype(np.float64).copy()
    for k in range(eff.shape[1]):
        e = eff[:, k]
        rho += h*(-rho*(1-rho)*(RHO_STAR_GB-rho)
                  + (e > tp)*GAMMA_P_GB*(1-rho) - (e > td)*(~(e > tp))*GAMMA_D_GB*rho)
        np.clip(rho, 0., 1., out=rho)
    return tau, np.abs(rho-obs).mean(), (rho-obs).mean(), 100*((rho>=.5)==(obs>=.5)).mean()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("nfiles", nargs="?", type=int, default=6)
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()

    U, dt, td, tp, r0, obs = load(a.nfiles)
    print(f"{len(U)} synapses from {a.nfiles} files, dt={dt:.4f} ms, "
          f"rho step {RHO_DT_MS} ms, {a.workers} workers\n")
    jobs = [(t, U, dt, td, tp, r0, obs) for t in TAUS]
    with ProcessPoolExecutor(min(a.workers, len(TAUS))) as ex:
        res = sorted(ex.map(score, jobs))

    print(f"{'tau':>10} {'mean|err|':>10} {'bias':>9} {'binary%':>9}")
    for tau, e, bias, bina in res:
        flag = "  <- current" if abs(tau-278.3177658387) < 1e-6 else ""
        print(f"{tau:>10.1f} {e:>10.4f} {bias:>+9.4f} {bina:>8.1f}%{flag}")
    best = min(res, key=lambda r: r[1])
    print(f"\nbest tau: {best[0]:.1f} ms   mean|err| {best[1]:.4f}   bias {best[2]:+.4f}")
