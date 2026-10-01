"""GPU objective v7x: gpu_v7_rho (not modified, stays reproducible) + the live arrival-time convention (V1 of
bcl_validation/diag_dt0_v7.py, BCL_W2_N3.md dt 0 diagnosis). SET t_exact 1 (needs veto_T > 0); without it every call
is the gpu_v7_rho path.

The v7 kernel reads the trigger pool at the first grid sample t_k at or after the arrival t_a, so it includes the
whole input of bin j = [t_k-1, t_k], also its part after t_a, which live GluSynapseV7 weights 0 (bglu = 1 from t_a on).
Its veto window starts at t_k and misses that part. v7x reads both at the exact arrival by linear sub-bin
interpolation, with f = (t_k - t_a) / h_j per arrival entry (earliest arrival of the sample; precomputed on the host,
aligned with astep like gpu_v7_rho.veto_q):
    trigger  W(t_a) = W(t_k) - f x (bV_j wb_j s_j)        (v5_mode 1: V(t_a) = V(t_k) - f x bV_j s_j)
    veto     Q(t_a) = Q(t_k) + f x bV_j s_j               (window [t_a, t_a + veto_T], as live)
Everything else is the v7 kernel. Left-point sampling inside bin j stays (the residual against live at dt 0, where
the bAP reaches the synapse 0.6-4 ms before the pre arrival and W rises steeply inside the bin).
Run: python gpu_v7x_rho.py <fit_v6 args> (as gpu_v7_rho.py; options in --set json).
"""
import math, os, sys
import numpy as np
from numba import cuda

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import gpu_v7_rho
from gpu_v7_rho import GPUModelV7, RHO_STAR_GB, KT
import model_v2 as MV


def _make_kernel_v7x(weighted=False, e2=False):
    """gpu_v7_rho._make_kernel_v7 plus the sub-bin arrival corrections (fa per arrival entry)."""
    RS = float(RHO_STAR_GB)
    TRIG_W = bool(weighted) or bool(e2)
    E2 = bool(e2)

    @cuda.jit
    def kern(E, S, Hg, eoff, slen, hoff, aptr, astep, acnt, rho0, order, td, tp, prm, out_rho, out_d, eptr, estep,
             uV, uE, vq, fa):
        p = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
        q = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
        if p >= td.shape[0] or q >= order.shape[0]:
            return
        i = order[q]
        GD = prm[p, 0]; GP = prm[p, 1]; thV0 = prm[p, 2]; iscV = prm[p, 3]; tauE1 = prm[p, 4]
        dmin = prm[p, 5]; Ae = prm[p, 6]; dp0 = prm[p, 7]
        thE0 = prm[p, 8]; tauD = prm[p, 9]
        thV = thV0 * uV[i]
        if thE0 <= 0.0:
            thE = thV
        else:
            thE = thE0 * uE[i]
        a_td = td[p, i]; a_tp = tp[p, i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        ip = 0.0; up = 0.0                                  # W / V input of the previous step
        e0 = eoff[i]; n = slen[i]; h0 = hoff[i]
        ka = aptr[i]; kae = aptr[i + 1]
        na = astep[ka] if ka < kae else -1
        ke = eptr[i]; kee = eptr[i + 1]
        ne = estep[ke] if ke < kee else -1
        gE = 1.0
        for k in range(n):
            e = E[e0 + k] * 1.0
            s = S[e0 + k] * 1.0
            h = Hg[h0 + k]
            hm = h * KT
            c = 0.0
            qa = 0.0
            fx = 0.0
            if k == na:
                c = acnt[ka]; qa = vq[ka]; fx = fa[ka]; ka += 1
                na = astep[ka] if ka < kae else -1
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            if pot > 0.0 and thV0 > 0.0 and not (V > thV):
                pot = 0.0
            r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - dep * (1 - pot) * GD * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if TRIG_W:
                if c > 0.0:
                    Bg = Bg + c
                trig = (W - fx * ip) > thE                   # W at the exact arrival
            else:
                trig = (V - fx * up) > thE
            if c > 0.0 and Ae > 0.0 and trig and not ((qa + fx * up) > thE):   # veto window from t_a
                d0 = dmin + (d0 - dmin) * math.pow(1.0 - Ae, c)
            if hm != chV:
                chV = hm
                aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / iscV
            if TRIG_W:
                wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
                if E2:
                    if k == ne:
                        gE = 0.0 if V > thV else 1.0
                        ke += 1
                        ne = estep[ke] if ke < kee else -1
                    ip = bV * (wb * gE * s)
                else:
                    ip = bV * (wb * s)
                W = aV * W + ip
                if hm != chD:
                    chD = hm
                    aD = math.exp(-hm / tauD)
                Bg = Bg * aD
            up = bV * s
            V = aV * V + up
        out_rho[p, i] = r
        out_d[p, i, 0] = d0; out_d[p, i, 1] = dp0; out_d[p, i, 2] = d0

    return kern


def arrival_frac(B):
    """Per arrival entry (gpu_v7_rho.veto_q order: records, synapses, ascending arrival samples): f = (t_k - t_a)/h_j,
    t_a = the earliest arrival (prespike + delay) mapped to sample k = arr, h_j = t_k - t_k-1; 0 at k = 0; in [0, 1]."""
    out = []
    for r in B.recs:
        t = np.asarray(r["t"], np.float64); n, T = r["effcai"].shape
        arr = np.asarray(r["arr"]).reshape(n, -1)
        delay = MV.edge_params(r["syn"])["delay"] if "delay" not in r else r["delay"]
        ta = np.asarray(r["prespikes"], np.float64).ravel()[None, :] + np.asarray(delay, np.float64)[:, None]
        for ii in range(n):
            k = arr[ii]; tai = ta[ii]
            ks = np.unique(k)                                    # = flatnonzero(cnt[ii]) of veto_q
            first = np.full(T, np.inf); np.minimum.at(first, k, tai)
            tk = t[ks]; tj = t[np.maximum(ks - 1, 0)]
            hj = tk - tj
            f = np.where(hj > 0, (tk - first[ks]) / np.where(hj > 0, hj, 1.0), 0.0)
            out.append(np.clip(f, 0.0, 1.0))
    f = np.concatenate(out) if out else np.zeros(0)
    return f if len(f) else np.zeros(1)


class GPUModelV7X(GPUModelV7):
    def __init__(self, B, targets, filters, **kw):
        super().__init__(B, targets, filters, **kw)
        self.t_exact = int(self.P0.get("t_exact", 0))
        if self.t_exact:
            assert self.veto_T > 0, "t_exact needs veto_T > 0 (the v7 veto path)"
            fa = arrival_frac(B)
            na = int(np.asarray(self.d["aptr"].copy_to_host())[-1])
            assert len(fa) == max(na, 1), (len(fa), na)
            self.fa = cuda.to_device(np.ascontiguousarray(fa, np.float64))
            self._kern7x = _make_kernel_v7x(weighted=self.v5_mode == 2, e2=self.v5_mode == 3)
            print(f"v7x t_exact: {na} arrivals, f q10/50/90 {np.quantile(fa, [.1, .5, .9]).round(3)}, "
                  f"f > 0 in {np.mean(fa > 0):.3f}", flush=True)

    def rho_dpre_syn(self, td, tp, Ps):
        if not (self.v5 and self.veto_T > 0 and self.t_exact):
            return super().rho_dpre_syn(td, tp, Ps)
        P, N = td.shape
        out_rho = cuda.device_array((P, N), np.float64); out_d = cuda.device_array((P, N, 3), np.float64)
        bx, by = self.block
        g = self.d
        self._kern7x[(-(-P // bx), -(-N // by)), (bx, by)](
            g["E"], g["S"], g["Hg"], g["eoff"], g["slen"], g["hoff"], g["aptr"], g["astep"], g["acnt"], g["rho0"],
            g["order"], cuda.to_device(np.ascontiguousarray(td, np.float64)),
            cuda.to_device(np.ascontiguousarray(tp, np.float64)), cuda.to_device(self._params5(Ps)), out_rho, out_d,
            *self.e5, *self.u6, self.vq, self.fa)
        return out_rho.copy_to_host(), out_d.copy_to_host()


if __name__ == "__main__":
    import json
    _s = json.loads(sys.argv[sys.argv.index("--set") + 1]) if "--set" in sys.argv else {}
    if int(_s.get("ecb_ref", 0)) == 2:                  # as gpu_v7_rho.__main__
        os.environ.setdefault("CEXP_DIR", "/scratch/dhuruva/split1/cexp")
        _e = os.environ.get("SCALE_E") or '{"vdcc_q_post": 100000.0}'
        assert json.loads(_e) == {"vdcc_q_post": 1e5}, f"ecb_ref 2 needs SCALE_E vdcc_q_post x 1e5, got {_e}"
        os.environ["SCALE_E"] = _e
    elif int(_s.get("ecb_ref", 0)) == 1:
        assert not os.environ.get("SCALE_E"), "ecb_ref 1 is not combined with SCALE_E"
    import fit_v6
    fit_v6.gpu_v6_rho.GPUModelV6 = GPUModelV7X
    fit_v6.main()
