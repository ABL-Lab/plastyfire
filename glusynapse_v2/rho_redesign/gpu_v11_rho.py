"""GPU objective v11 = gpu_v10_rho (copy; v10 not modified) + eCB veto window, peak veto and the MVD depression arm (SPEC_V11.md).
With every v11 option at its default each call is the v10 kernel operation for operation (v10 doc below).

v11 options (SET; all synapse-local, uniform across pathways):
  veto_t0 (ms, 0), veto_Tv (ms, = veto_T), veto_peak (0): eCB veto window (t_a + veto_t0, t_a + veto_Tv] on the own spine VDCC
           influx (t_a the exact arrival with t_exact 1, else the arrival step). Any non-default value moves the veto into the
           kernel (summed per candidate, overlap-weighted steps). veto_peak 0: Q = sum bV_j s_j over the window > theta_eCB,i
           (as v7). veto_peak 1: Q = max s_j over the window > veto_peak_k x I1_i, I1_i = own single-bAP PEAK VDCC current
           (step mean, median over isolated bAPs of the records, window [x - 1, x + 8] ms, v7 ecb_ref 1 isolation; tier 1 no own
           arrival in [x - 150, x + 8], tier 2 any; NaN -> pooled median). veto_peak_k default 0.5 (placeholder; free-able).
           Spine eCB only (ecb_src 0).
  (1 - b)-unweighted eCB trigger: already v5_mode 1 (TRIG_W off: trigger on V, the unweighted pool; nothing else changes).
  mvd_mode 1: MVD ("mGluR-gated VDCC depression"). M = [lo < V <= hi P1_i] [b > b_c], V the unweighted own spine VDCC
           pool (tau_E1), b = min(Bg, 1) the own glutamate-bound state (jump 1 per own arrival, so b_peak = 1; decay
           tau_d_NMDA). b_c = exp(-mvd_W / tau_d_NMDA), mvd_W 50 ms (mGluR1 deactivation, Marcaggi 2009: 0.49 at 70 ms);
           mvd_b >= 0 overrides b_c. Depression indicator D = max(dep, M); mvd_pot 0: rho' = ... - D (1 - pot) gamma_d rho
           (MVD silent while licensed pot, like dep); mvd_pot 1: rho' = ... - max(dep (1 - pot), M) gamma_d rho (gamma_d once).
           lo: mvd_ref 0 = theta_MVD_lo x P1_i (ecb_ref 2 single-bAP pool; 0.5 placeholder); 1 = theta_MVD_lo x P1_i
           sum_{k<n} e^(-k isi / tau_E1) (mvd_ref_n 10, mvd_ref_isi 20 ms); 2 = theta_MVD_abs, absolute, in c_VDCC pool units
           (1e5 x VDCC charge, the P1_i units; placeholder 100, box 1-1e4 log; set from C_MVD_CALIB.md). hi = theta_MVD_hi x
           P1_i (1e30 = no upper edge, ANCHORS_V11). All free-able (--free-filters). Needs ecb_ref 2.
  Counter 6 (v10_counts 1): MVD-depressing steps. Saved json gets a "v11" block.
  ecb_lp_mode (0; SPEC_ECB_LP.md): eCB trigger read from L, a first-order low-pass (ecb_lp_tau, default 10 ms, Ebner 2019
           pre-LTD) of the same own spine VDCC input as the trigger pool (s, or (1 - b) s with v5_mode 2), in pool units
           (L' = -L / ecb_lp_tau + x / i_scale, so one isolated bAP lifts L by ~P1_i). At each own arrival (exact with t_exact
           1): trigger = L(t_a) > theta_eCB x P1_i (thE, no new threshold). 1: no veto (Ebner); 2: the configured veto (v7
           [t_a, t_a + veto_T] integral or the v11 window) is kept. Spine eCB only (ecb_src 0). 0 = v11 unchanged.

v10 options (SET):
  lic_src  LTP licence; a potentiating crossing without the licence counts as depression (v8 convention).
           absent: v8 gate_src semantics (0 = Vg, the v7x path; 1/2/3 as v8).
           0  none (pot never gated);  1/2/3  as v8 gate_src 1/2/3 (gate_win, gate_theta, kappa);
           4  peak own spine Ca above rest (cai_CR - min_ca_CR, from effcai as v8 src 1, uM) over [t - lic_win, t]
              >= theta_L,i = a20 Cpre_i + a21 Cpost_i, Cpre/Cpost = cexp cpre_cai / cpost_cai - 0.07 uM (single-EPSP /
              single-bAP peak free spine Ca above rest; CEXP_DIR);
           5  peak own shaft dCa above rest (v8 src 2 signal, uM) over [t - lic_win, t] >= kappa_post bapsh_i + kappa_pre
              presh_i. bapsh_i = median peak shaft dCa over [x - 1, x + 30] ms of the synapse's bAPs with no other post
              spike in [x - 100, x + 20] ms; tier 1 also has no own arrival in [x - 150, x + 8] ms (v8 bap_sh), tier 2
              (used only for synapses without a tier-1 bAP) allows own arrivals, since the shaft carries no synaptic Ca
              (GATE_IMPL (2)) and the single-EPSP shaft dCa presh_i is ~0. presh_i: the same peak for isolated own arrivals
              (no post spike in [t_a - 300, t_a + 30], no other own arrival in [t_a - 100, t_a + 30]).
           6  peak c* (effcai) over [t - lic_win, t] >= theta_L,i = a20 c_pre,i + a21 c_post,i (Chindemi units).
           lic_win (ms, default 100; > 0 for 4 / 5; 6 also takes 0 = the current c* only). Free: a20, a21 (4, 6), kappa_pre, kappa_post (5).
           lic_lp 1 (lic_src 4 / 5): the licence quantity is a first-order low-pass of the same signal instead of the
           window peak, dL/dt = (x - L) / tau_L (tau_L default 100 ms, free; Nevian 2006), licence = [L >= theta_L,i].
  ecb_src  eCB-LTD trigger source (step, veto window veto_T and d_min as v8).
           0  v8: W (own spine VDCC pool, (1 - b) weighted), theta = k_E P1_i (theta_eCB x uE), veto on own VDCC influx.
           1  shaft pool W_sh: dW_sh/dt = ((1 - b) x_sh - W_sh) / tau_E1, x_sh = own shaft dCa above rest (uM, v8 src 2
              signal), so W_sh is a 100 ms low-pass in uM. Trigger at own arrival if W_sh > k_E W1_i, W1_i = peak W_sh of a
              bAP at this synapse (bapsh tiers, window [x - 1, min(x + 300, next post spike)) ms, unweighted, <= 5 per tier,
              median per synapse). Veto if the unweighted W_sh input over [t_a, t_a + veto_T], sum (1 - e^(-h/tau_E1)) x_sh,
              exceeds the same threshold (the shaft analogue of the v7 veto Q). k_E = theta_eCB.
           2  as 1 with the absolute threshold ecb_theta_uM (uM; settable or free).
  no_mode  NO-dependent presynaptic LTP (Sjostrom 2007). d = clip(d_eCB + d_NO, d_min, d_NO_max) (d_eCB the v8 eCB state),
           d_NO capped at d_NO_max. Lanes: control d; mglu_block clip(dpre0 + d_NO) (no eCB step, NO kept); no_block = d_eCB
           (L-NAME / cPTIO: NO off). nmdar_block stays d = 0 (gpu_v3), post_nmdar reads the control lane (no_drive 0 lane
           layout; not a fit target).
           0  off.
           1/2: N: dN/dt = -N / tau_NO + (cai_CR - min_ca_CR), own spine Ca (Chindemi c* units: N = c* at tau_NO 278.3 ms),
              theta_NO,i = a_NO_pre c_pre,i + a_NO_post c_post,i. 1 graded: at each own arrival d_NO += c A_NO max(0, N -
              theta_NO,i); 2 all-or-none: at each own arrival with N > theta_NO,i, d_NO += c dNO_step.
              Free: tau_NO (default 100), A_NO, dNO_step, d_NO_max, a_NO_pre, a_NO_post.
           3  licence-gated NO: N += c on each step with an own arrival (c) and licensed potentiation (pot after the
              licence), dN/dt = -N / tau_NO (tau_NO 6.7 ms, Hall & Garthwaite 2006), no NO threshold;
              dd_NO/dt = A_NO N in tau_ind units (A_NO a rate like gamma_p: one licensed arrival adds A_NO tau_NO / 70 s,
              exact per step), d_NO <= d_NO_max. Free: A_NO (box 32-1123, default 190), d_NO_max (box 0.4-1.0, default 0.5).
  tau_E1 free (--free-filters tau_E1, box 25-125 ms, every seed starts at tau_E1_start, default 50; 0 keeps the seed's):
           the v7 veto Q (and the ecb_src 1/2 shaft veto) is summed in the kernel per candidate instead of the host
           precompute at the filters' tau_E1, and ecb_src 1 W1_i is interpolated (log tau) from a host table at 25-125 ms.
  v10_counts 1: kernel counters, printed per model for single-candidate calls (the final table): licensed / blocked
           potentiation steps, eCB steps, vetoed eCB triggers, NO steps.
With t_exact 1 the shaft trigger / veto and N (modes 1/2) read the exact arrival (v7x sub-bin correction, as W / V). NaN
per-synapse calibrations (no event, NaN cexp row) get the pooled median over all models (fit_v6 rule). Every drive is the
synapse's own spine / shaft Ca or its own pre arrivals. Saved json gets a "v10" block. Run (v11): python gpu_v11_rho.py <fit_v6 args>.

GPU objective v8 (GATE_ALT.md, GATE_IMPL.md): v7 (gpu_v7_rho, not modified) with the LTP licence read from a chosen
synapse-local Ca quantity. Option gate_src in --set json:
  0  Vg (own spine VDCC pool, tau_E1): every call is the v7 code path unchanged (W2_N3 repro at MAXITER 0).
  1  spine [Ca] above rest (cai_CR - min_ca_CR, recovered exactly from effcai, batch_v2.cacr_from_effcai), uM.
  2  shaft [Ca] above the synapse's pre-stimulus rest (npz shaft_cai = cai of the segment carrying the synapse), uM.
  3  c* (effcai) itself, uM s (mM ms), licence = [c* > theta_G].
  gate_win > 0 (ms, src 1/2): peak-in-window, licence = [max of the signal over [t - gate_win, t] > theta_G] (GATE_ALT G1);
  gate_win 0 (src 1/2): 100 ms low-pass, dG/dt = (x - G) / tau_E1, licence = [G > theta_G] (GATE_ALT G2).
Threshold theta_G (natural unit, one value for every pathway): gate_theta > 0 fixes it (theta_V then unused); otherwise
theta_G = kappa x theta_V (DE slot theta_V), kappa = Q50 / gate_conv_thV, Q50 = median over synapses (pooled over all
models) of each synapse's median gate quantity at the moments the Vg gate opens (Vg crosses gate_conv_thV, default W2_N3
theta_V 3.575), so the W2_N3 seed starts at the converted threshold of the median synapse. gate_kappa > 0 pins kappa.
The saved json gets a "v8" block (gate_src, gate_win, kappa, Q50, theta_G and unit). For gate_src != 0 the fit_v6
repro check (seed rule spec) fails by design; that REPRO DIFF exit is caught here and the run exits 0.
src 2 calibration print (GATE_ALT checklist): single isolated bAP peak shaft Delta[Ca] per synapse (bap_ref selection,
peak over [x - 1, x + 30] ms) by cexp loc / proximal (<= 100 um), Spearman vs cexp vdcc_q_post; 3AP 200 Hz record peaks.
t_exact 1 (SET): GPUModelV8 inherits gpu_v7x_rho.GPUModelV7X; the eCB trigger and veto use the v7x exact-arrival corrections
(live convention); t_exact 0 = the v7 arithmetic unchanged (G8 results bit-for-bit). v5_mode 1/2 only (no E2 latch). Shaft is stored on the GPU as uint16 nM above rest (clip 0..65535 nM).
"""
import json, math, os, sys
import numpy as np
from numba import cuda, njit

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import gpu_v7_rho, gpu_v7x_rho
from gpu_v7_rho import KT
from gpu_v7x_rho import GPUModelV7X, arrival_frac
from batch_v2 import RHO_STAR_GB, TAU_EFFCA

UNIT = {1: "uM", 2: "uM", 3: "uM s"}
_INST = []
_STATE = dict(kappa=None, fin=False, te1_free=False)
CA_REST_UM = 0.07                      # min_ca_CR (GluSynapse*.mod 70e-6 mM); cexp *_cai are absolute peaks
# v10 parameters: defaults (used only by the option that reads them), DE boxes (fit_v2 FILTER_BOX convention), log-spaced ones
V10_DEFAULTS = dict(lic_win=100.0, a20=0.0, a21=1.0, kappa_pre=0.0, kappa_post=0.3, ecb_theta_uM=0.05, tau_L=100.0,
                    tau_NO=100.0, A_NO=0.1, dNO_step=0.05, d_NO_max=0.5, a_NO_pre=2.0, a_NO_post=4.0,
                    theta_MVD_lo=0.5, theta_MVD_hi=1e30, veto_peak_k=0.5, theta_MVD_abs=100.0)            # v11 (placeholders, SPEC_V11.md)
V10_NO3 = dict(tau_NO=6.7, A_NO=190.0, d_NO_max=0.5)          # no_mode 3 defaults (Hall & Garthwaite 2006 NO lifetime)
V10_BOX = dict(a20=(0.0, 5.0), a21=(0.0, 5.0), kappa_pre=(0.0, 3.0), kappa_post=(0.0, 3.0), ecb_theta_uM=(1e-3, 2.0),
               tau_L=(10.0, 1000.0), tau_NO=(10.0, 1000.0), A_NO=(1e-3, 10.0), dNO_step=(1e-3, 0.5), d_NO_max=(0.0, 2.0),
               a_NO_pre=(0.0, 20.0), a_NO_post=(0.0, 20.0), tau_E1=(25.0, 125.0),
               theta_MVD_lo=(0.2, 1.0), theta_MVD_hi=(0.5, 10.0), veto_peak_k=(0.1, 3.0), theta_MVD_abs=(1.0, 1e4))
V10_BOX_NO3 = dict(A_NO=(32.0, 1123.0), d_NO_max=(0.4, 1.0))
V10_LOG = {"ecb_theta_uM", "tau_L", "tau_NO", "A_NO", "dNO_step", "theta_MVD_abs"}
TE1_GRID = np.array([25.0, 35.0, 50.0, 70.0, 100.0, 125.0])  # tau_E1 free: ecb_src 1 W1_i table (log-interpolated)
NCNT = 6                               # counters: licensed pot steps, blocked pot steps, eCB steps, vetoed triggers, NO steps,
                                       # v11 MVD-depressing steps
V11_CONST = dict(veto_t0=0.0, veto_peak=0, mvd_mode=0, mvd_pot=0, mvd_ref=0, mvd_ref_n=10, mvd_ref_isi=20.0, mvd_b=-1.0, mvd_W=50.0,
                 ecb_lp_mode=0, ecb_lp_tau=10.0)  # eCB LP: tau 10 ms = Ebner 2019 pre-LTD low-pass (SPEC_ECB_LP.md)
V11_KEYS = ("theta_MVD_lo", "theta_MVD_hi", "veto_peak_k", "theta_MVD_abs")   # v11 per-candidate parameters (prm slots 20-23)


def v11_opts(P, veto_T):
    """v11 kernel constants from a SET dict; vw = the veto window differs from the v7 [t_a, t_a + veto_T] integral."""
    o = {k: type(v)(P.get(k, v)) for k, v in V11_CONST.items()}
    o["veto_Tv"] = float(P.get("veto_Tv", veto_T))
    o["vw"] = bool(o["veto_t0"] != 0.0 or o["veto_Tv"] != float(veto_T) or o["veto_peak"])
    if o["veto_t0"] < 0.0 or (o["vw"] and o["veto_Tv"] > 0.0 and o["veto_Tv"] <= o["veto_t0"]):
        raise ValueError(f"v11 veto window ({o['veto_t0']}, {o['veto_Tv']}]")
    if o["mvd_mode"] not in (0, 1) or o["mvd_ref"] not in (0, 1, 2) or o["veto_peak"] not in (0, 1) or o["mvd_pot"] not in (0, 1):
        raise ValueError(f"v11 options {o}")
    if o["ecb_lp_mode"] not in (0, 1, 2) or not o["ecb_lp_tau"] > 0.0:
        raise ValueError(f"v11 ecb_lp_mode {o['ecb_lp_mode']} ecb_lp_tau {o['ecb_lp_tau']}")
    if o["ecb_lp_mode"] == 1:                    # no veto at all: the window options are not read
        o["vw"] = False
    return o


def v10_opts(P):
    """(LIC, ECB, NOM) kernel codes from a SET / filters dict. LIC -1 = Vg licence (v8 gate_src 0)."""
    gs = int(P.get("gate_src", 0))
    lic = int(P["lic_src"]) if P.get("lic_src") is not None else (gs if gs else -1)
    ecb, nom = int(P.get("ecb_src", 0)), int(P.get("no_mode", 0))
    if lic not in (-1, 0, 1, 2, 3, 4, 5, 6) or ecb not in (0, 1, 2) or nom not in (0, 1, 2, 3):
        raise ValueError(f"v10 lic_src {lic} ecb_src {ecb} no_mode {nom}")
    return lic, ecb, nom


def v10_default(k, nom):
    return V10_NO3[k] if nom == 3 and k in V10_NO3 else V10_DEFAULTS[k]


def v10_keys(lic, ecb, nom, llp=0, mvd=0, vpk=0, mref=0):
    """The v10 / v11 parameters read by the active options (others are never read)."""
    k = []
    if mvd:
        k += ["theta_MVD_abs" if mref == 2 else "theta_MVD_lo", "theta_MVD_hi"]
    if vpk:
        k += ["veto_peak_k"]
    if lic in (4, 6):
        k += ["a20", "a21"]
    elif lic == 5:
        k += ["kappa_pre", "kappa_post"]
    if llp and lic in (4, 5):
        k += ["tau_L"]
    if ecb == 2:
        k += ["ecb_theta_uM"]
    if nom == 3:
        k += ["tau_NO", "A_NO", "d_NO_max"]
    elif nom:
        k += ["tau_NO", "A_NO" if nom == 1 else "dNO_step", "d_NO_max", "a_NO_pre", "a_NO_post"]
    return k


def _make_kernel_v10(weighted, lic, win, tex=False, ecb=0, nom=0, cnt=False, llp=False, tef=False, tv=0.0, v11=None):
    """_make_kernel_v8 (gpu_v8_rho) with lic (v10 LIC code), ecb (eCB trigger source), nom (NO arm), llp (low-pass licence),
    tef (tau_E1 free: veto sums in the kernel, tv = veto_T ms); see module doc. lic 1/2/3, ecb 0, nom 0, llp / tef off is the
    v8 kernel operation for operation. v11 (v11_opts dict): veto window / peak veto (vw) and MVD; None or defaults = v10."""
    o = dict(V11_CONST, veto_Tv=tv, vw=False) if v11 is None else v11
    VW = bool(o["vw"]) and int(ecb) == 0
    VT0 = float(o["veto_t0"]); VTV = float(o["veto_Tv"]); VPK = bool(o["veto_peak"]) and VW
    MVD = bool(o["mvd_mode"]); MVP = bool(o["mvd_pot"]); MREF = int(o["mvd_ref"])
    MRN = int(o["mvd_ref_n"]); MRI = float(o["mvd_ref_isi"]); MB = float(o["mvd_b"]); MW = float(o["mvd_W"])
    LPM = int(o["ecb_lp_mode"]) if int(ecb) == 0 else 0   # eCB LP trigger (1 no veto, 2 configured veto kept)
    TLP = float(o["ecb_lp_tau"])
    BGT = bool(weighted) or MVD            # own glutamate-bound state Bg tracked
    RS = float(RHO_STAR_GB)
    TRIG_W = bool(weighted)
    SRC = int(lic)
    WIN = float(win)
    TS = float(TAU_EFFCA)
    TEX = bool(tex)                       # v7x arrival-time convention (gpu_v7x_rho, SET t_exact 1)
    ECB = int(ecb)
    NOM = int(nom)
    CNT = bool(cnt)
    LLP = bool(llp) and (SRC == 4 or SRC == 5)
    TEF = bool(tef)
    TV = float(tv)
    NG = len(TE1_GRID)
    SPN = SRC == 1 or SRC == 4 or NOM == 1 or NOM == 2   # spine Ca above rest needed (licence 1 / 4 or the NO pool)
    SHX = SRC == 2 or SRC == 5 or ECB > 0                # shaft dCa trace needed

    @cuda.jit
    def kern(E, S, Hg, eoff, slen, hoff, aptr, astep, acnt, rho0, order, td, tp, prm, out_rho, out_d, uV, uE, vq,
             X, gsc, gfix, fa, q1, q2, uEs, vqs, cpr, cpo, ocnt, w1g, ltg, spk):
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
        thG = gfix if gfix > 0.0 else thV0 * gsc
        thL = prm[p, 10] * q1[i] + prm[p, 11] * q2[i]            # v10 licence 4 / 5 / 6
        if ECB == 2:
            thEs = prm[p, 12]
        elif TEF and ECB == 1:                                   # W1_i at this candidate's tau_E1 (log interpolation)
            lt = math.log(tauE1)
            if lt <= ltg[0]:
                w1 = w1g[i, 0]
            elif lt >= ltg[NG - 1]:
                w1 = w1g[i, NG - 1]
            else:
                g = 0
                while ltg[g + 1] < lt:
                    g += 1
                u = (lt - ltg[g]) / (ltg[g + 1] - ltg[g])
                w1 = w1g[i, g] + u * (w1g[i, g + 1] - w1g[i, g])
            thEs = thE0 * w1
        else:
            thEs = thE0 * uEs[i]
        tauN = prm[p, 13]; Ano = prm[p, 14]; dstep = prm[p, 15]; dnmax = prm[p, 16]
        thN = prm[p, 17] * cpr[i] + prm[p, 18] * cpo[i]
        tauL = prm[p, 19]
        thMlo = 0.0; thMhi = 0.0; thP = 0.0; bcM = 0.0
        if MVD:                                                  # v11 MVD band on the own unweighted VDCC pool
            Rm = uE[i]
            if MREF == 1:                                        # own control-train pool, linear sum of P1_i at tau_E1
                gs = 0.0
                for kk in range(MRN):
                    gs += math.exp(-kk * MRI / tauE1)
                Rm = Rm * gs
            if MREF == 2:                                        # absolute lower edge, c_VDCC pool units (1e5 x charge)
                thMlo = prm[p, 23]
            else:
                thMlo = prm[p, 20] * Rm
            thMhi = prm[p, 21] * Rm
            bcM = MB if MB >= 0.0 else math.exp(-MW / tauD)      # own-glutamate cutoff: b > b_peak e^(-W / tau_d_NMDA)
        if VPK:                                                  # v11 peak veto threshold, peak units
            thP = prm[p, 22] * spk[i]
        hpv = 0.0
        a_td = td[p, i]; a_tp = tp[p, i]
        r = rho0[i]; d0 = dp0
        V = 0.0; chV = -1.0; aV = 1.0; bV = 0.0
        W = 0.0; Bg = 0.0; chD = -1.0; aD = 1.0
        Gq = 0.0; clk = 0.0; tl = -1e30; chS = -1.0; aS = 1.0; bS = 0.0
        ip = 0.0; up = 0.0                                      # v7x: W / V input of the previous step
        Ws = 0.0; ips = 0.0; ups = 0.0                          # v10 shaft eCB pool and its last inputs
        Nn = 0.0; inN = 0.0; chN = -1.0; aN = 1.0; bN = 0.0; dn = 0.0   # v10 NO pool, its last input, NO state of d
        Gl = 0.0; chL = -1.0; aL = 1.0                          # v10 low-pass licence quantity
        Lp = 0.0; upL = 0.0; chP = -1.0; aP = 1.0; bP = 0.0      # eCB LP pool (ecb_lp_mode) and its last input
        n_lo = 0.0; n_lb = 0.0; n_ec = 0.0; n_ve = 0.0; n_no = 0.0; n_mv = 0.0
        e0 = eoff[i]; n = slen[i]; h0 = hoff[i]
        ka = aptr[i]; kae = aptr[i + 1]
        na = astep[ka] if ka < kae else -1
        for k in range(n):
            e = E[e0 + k] * 1.0
            s = S[e0 + k] * 1.0
            h = Hg[h0 + k]
            hm = h * KT
            c = 0.0
            qa = 0.0
            fx = 0.0
            qs = 0.0
            if k == na:
                c = acnt[ka]; qa = vq[ka]
                if TEX:
                    fx = fa[ka]
                if ECB > 0:
                    qs = vqs[ka]
                ka += 1
                na = astep[ka] if ka < kae else -1
                if VW:                                           # v11 veto window (t_a + VT0, t_a + VTV], this tau_E1
                    qa = 0.0
                    off = fx * hpv if TEX else 0.0               # t_k - t_a
                    if off > 0.0 and k > 0:                      # the part of step k - 1 after t_a
                        lo = VT0
                        hi = off if off < VTV else VTV
                        if hi > lo:
                            if VPK:
                                qa = S[e0 + k - 1] * 1.0
                            else:
                                qa = up * (hi - lo) / hpv
                    cj = off; j = k
                    while j < n and cj < VTV:
                        hj = Hg[h0 + j] * KT
                        lo = cj if cj > VT0 else VT0
                        hi = (cj + hj) if (cj + hj) < VTV else VTV
                        if hj > 0.0 and hi > lo:
                            sj = S[e0 + j] * 1.0
                            if VPK:
                                if sj > qa:
                                    qa = sj
                            else:
                                aj = math.exp(-hj / tauE1)
                                qa += tauE1 * (1.0 - aj) / iscV * sj * ((hi - lo) / hj)
                        cj += hj
                        j += 1
                elif TEF:                                        # veto inputs over [t_k, t_k + veto_T] at this tau_E1
                    qa = 0.0; qs = 0.0
                    if TV > 0.0:
                        cj = 0.0; j = k
                        while j < n and cj <= TV:
                            hj = Hg[h0 + j] * KT
                            aj = math.exp(-hj / tauE1)
                            qa += tauE1 * (1.0 - aj) / iscV * (S[e0 + j] * 1.0)
                            if ECB > 0:
                                qs += (1.0 - aj) * (X[e0 + j] * 1e-3)
                            cj += hj
                            j += 1
            pot = 1.0 if e > a_tp else 0.0
            dep = 1.0 if e > a_td else 0.0
            if SRC == 6 and e >= thL:
                tl = clk
            if SRC == -1:
                lic = V > thV
            elif SRC == 0:
                lic = True
            elif SRC == 3:
                lic = e > thG
            elif LLP:
                lic = Gl >= thL
            elif WIN > 0.0:
                lic = (clk - tl) <= WIN
            elif SRC == 1 or SRC == 2:
                lic = Gq > thG
            else:
                lic = (clk - tl) <= 0.0
            pq = pot
            if SRC == -1:
                if pot > 0.0 and thV0 > 0.0 and not lic:
                    pot = 0.0
            elif SRC >= 4:
                if pot > 0.0 and not lic:
                    pot = 0.0
            elif SRC > 0:
                if pot > 0.0 and thG > 0.0 and not lic:
                    pot = 0.0
            if CNT and pq > 0.0:
                if pot > 0.0:
                    n_lo += 1.0
                else:
                    n_lb += 1.0
            if BGT and c > 0.0:                                      # own glutamate-bound state (v5b b; read by W and MVD)
                Bg = Bg + c
            if MVD:                                                  # v11 MVD: own VDCC pool in band and own glutamate
                mv = 0.0
                if V > thMlo and V <= thMhi and Bg > bcM:
                    mv = 1.0
                if MVP:
                    Dd = dep * (1 - pot)
                    if mv > Dd:
                        Dd = mv
                        if CNT:
                            n_mv += 1.0
                    r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - Dd * GD * r)
                else:
                    Dd = dep
                    if mv > Dd:
                        Dd = mv
                        if CNT and pot == 0.0:
                            n_mv += 1.0
                    r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - Dd * (1 - pot) * GD * r)
            else:
                r = r + h * (-r * (1 - r) * (RS - r) + pot * GP * (1 - r) - dep * (1 - pot) * GD * r)
            if r < 0.0:
                r = 0.0
            elif r > 1.0:
                r = 1.0
            if TRIG_W:
                trig = (W - fx * ip) > thE if TEX else W > thE      # v7x: W at the exact arrival
            else:
                trig = (V - fx * up) > thE if TEX else V > thE
            if TEX and not VW:
                qa = qa + fx * up                                    # v7x: veto window from t_a
            if ECB > 0:
                trig = (Ws - fx * ips) > thEs if TEX else Ws > thEs  # v10: shaft pool at the (exact) arrival
                if TEX:
                    qs = qs + fx * ups
                vet = qs > thEs
            elif VPK:
                vet = qa > thP                                       # v11 peak veto, own single-bAP peak units
            else:
                vet = qa > thE
            if LPM > 0:                                              # eCB LP: own VDCC low-pass at the (exact) arrival
                trig = (Lp - fx * upL) > thE if TEX else Lp > thE
                if LPM == 1:
                    vet = False
            if c > 0.0 and Ae > 0.0 and trig and not vet:
                d0 = dmin + (d0 - dmin) * math.pow(1.0 - Ae, c)
                if CNT:
                    n_ec += 1.0
            elif CNT and c > 0.0 and Ae > 0.0 and trig:
                n_ve += 1.0
            if (NOM == 1 or NOM == 2) and c > 0.0:                   # v10: NO step at the own arrival
                Nx = (Nn - fx * inN) if TEX else Nn
                if NOM == 1:
                    ex = Nx - thN
                    if ex > 0.0 and Ano > 0.0:
                        dn = dn + c * Ano * ex
                        if CNT:
                            n_no += 1.0
                elif Nx > thN:
                    dn = dn + c * dstep
                    if CNT:
                        n_no += 1.0
                if dn > dnmax:
                    dn = dnmax
            if NOM == 3 and c > 0.0 and pot > 0.0:                   # v10: NO produced by a licensed own arrival
                Nn = Nn + c
                if CNT:
                    n_no += 1.0
            if hm != chV:
                chV = hm
                aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / iscV
            if TRIG_W:
                wb = 1.0 - (Bg if Bg < 1.0 else 1.0)
                if TEX:
                    ip = bV * (wb * s)
                    W = aV * W + ip
                else:
                    W = aV * W + bV * (wb * s)
            else:
                wb = 1.0
            if BGT:
                if hm != chD:
                    chD = hm
                    aD = math.exp(-hm / tauD)
                Bg = Bg * aD
            if TEX:
                up = bV * s
                V = aV * V + up
            else:
                V = aV * V + bV * s
            if LPM > 0:                                              # same input as the trigger pool (wb = 1 - b or 1)
                if hm != chP:
                    chP = hm
                    aP = math.exp(-hm / TLP); bP = TLP * (1.0 - aP) / iscV
                upL = bP * (wb * s)
                Lp = aP * Lp + upL
            xs = 0.0
            if SPN:
                if hm > 0.0 and k + 1 < n:
                    if hm != chS:
                        chS = hm
                        aS = math.exp(-hm / TS); bS = 1e3 / (TS * (1.0 - aS))
                    xs = (E[e0 + k + 1] * 1.0 - aS * e) * bS            # step mean of cai_CR - min_ca_CR, uM
            xh = 0.0
            if SHX:
                xh = X[e0 + k] * 1e-3                                    # nM -> uM
            if ECB > 0:
                ups = (1.0 - aV) * xh
                ips = wb * ups
                Ws = aV * Ws + ips
            if NOM == 1 or NOM == 2:
                if hm != chN:
                    chN = hm
                    aN = math.exp(-hm / tauN)
                inN = tauN * (1.0 - aN) * (xs * 1e-3)                   # mM ms, c* units
                Nn = aN * Nn + inN
            elif NOM == 3:
                if hm != chN:
                    chN = hm
                    aN = math.exp(-hm / tauN); bN = tauN * (1.0 - aN) / KT
                dn = dn + Ano * Nn * bN                                 # exact integral of N over the step, tau_ind units
                if dn > dnmax:
                    dn = dnmax
                Nn = aN * Nn
            if LLP:
                if hm != chL:
                    chL = hm
                    aL = math.exp(-hm / tauL)
                Gl = aL * Gl + (1.0 - aL) * (xs if SRC == 4 else xh)
            if SRC == 1 or SRC == 2:
                if SRC == 1:
                    x = xs
                    tx = clk + hm
                else:
                    x = xh
                    tx = clk
                if WIN > 0.0:
                    if x > thG:
                        tl = tx
                else:
                    Gq = aV * Gq + (1.0 - aV) * x
            elif SRC == 4 and not LLP:
                if xs >= thL:
                    tl = clk + hm
            elif SRC == 5 and not LLP:
                if xh >= thL:
                    tl = clk
            clk += hm
            hpv = hm
        out_rho[p, i] = r
        if NOM > 0:
            dc = d0 + dn
            if dc < dmin:
                dc = dmin
            elif dc > dnmax:
                dc = dnmax
            dm = dp0 + dn
            if dm < dmin:
                dm = dmin
            elif dm > dnmax:
                dm = dnmax
            out_d[p, i, 0] = dc; out_d[p, i, 1] = dm; out_d[p, i, 2] = d0
        else:
            out_d[p, i, 0] = d0; out_d[p, i, 1] = dp0; out_d[p, i, 2] = d0
        if CNT:
            ocnt[p, i, 0] = n_lo; ocnt[p, i, 1] = n_lb; ocnt[p, i, 2] = n_ec; ocnt[p, i, 3] = n_ve; ocnt[p, i, 4] = n_no
            ocnt[p, i, 5] = n_mv

    return kern


@njit(cache=False)
def _open_q(t, s, e, x, src, win, tauE1, isc, thV, TS):
    """Gate quantity at each Vg gate opening (first step k with V_k > thV after V_{k-1} <= thV), as the kernel sees it."""
    n = t.shape[0]
    res = np.empty(n); m = 0
    xs = np.zeros(n); txs = np.zeros(n)
    V = 0.0; Gq = 0.0; clk = 0.0; prev = False
    for k in range(n):
        hm = (t[k + 1] - t[k]) if k + 1 < n else 0.0
        op = V > thV
        if op and not prev:
            if src == 3:
                qv = e[k]
            elif win > 0.0:
                qv = 0.0; j = k - 1
                while j >= 0 and txs[j] >= clk - win:
                    if xs[j] > qv:
                        qv = xs[j]
                    j -= 1
            else:
                qv = Gq
            res[m] = qv; m += 1
        prev = op
        aV = math.exp(-hm / tauE1); bV = tauE1 * (1.0 - aV) / isc
        V = aV * V + bV * s[k]
        xk = 0.0; tx = clk
        if src == 1:
            if hm > 0.0 and k + 1 < n:
                aS = math.exp(-hm / TS); xk = (e[k + 1] - aS * e[k]) * 1e3 / (TS * (1.0 - aS))
            tx = clk + hm
        elif src == 2:
            xk = x[k]
        xs[k] = xk; txs[k] = tx
        Gq = aV * Gq + (1.0 - aV) * xk
        clk += hm
    return res[:m]


@njit(cache=False)
def _lp_peak(x, aw):
    """Peak of W, W <- (1 - aw_k) W + aw_k x_k from W = 0 (the v10 shaft pool of one bAP, unweighted)."""
    w = 0.0; m = 0.0
    for k in range(x.shape[0]):
        w = (1.0 - aw[k]) * w + aw[k] * x[k]
        if w > m:
            m = w
    return m


def _shaft_rest(r, sh):
    """Per-synapse shaft cai (mM) at the last grid sample >= 1 ms before the first pre or post spike."""
    t = np.asarray(r["t"], np.float64)
    ev = np.concatenate([np.asarray(r["prespikes"], float).ravel(), np.asarray(r["postspikes"], float).ravel()])
    k = 0
    if ev.size:
        k = max(int(np.searchsorted(t, ev.min() - 1.0)) - 1, 0)
    return sh[:, k].astype(np.float64)


def _iso_pre(ta, post, nmax=5):
    """Up to nmax isolated own arrival times: no post spike in [t_a - 300, t_a + 30], no other arrival in [t_a - 100, t_a + 30]."""
    u = np.unique(ta)
    if not u.size:
        return u
    gp = np.diff(u, prepend=-np.inf); gn = np.diff(u, append=np.inf)
    ok = (gp > 100.0) & (gn > 30.0)
    if post.size:
        ok &= np.searchsorted(post, u - 300.0, side="left") == np.searchsorted(post, u + 30.0, side="right")
    return u[ok][:nmax]


class GPUModelV11(GPUModelV7X):
    def __init__(self, B, targets, filters, **kw):
        super().__init__(B, targets, filters, **kw)
        P0 = self.P0
        self.f10 = dict(filters)
        self.gate_src = int(P0.get("gate_src", 0))
        self.lic, self.ecb, self.nom = v10_opts(P0)
        self.cnt = bool(int(P0.get("v10_counts", 0)))
        self.llp = bool(int(P0.get("lic_lp", 0))) and self.lic in (4, 5)
        self.tef = bool(_STATE["te1_free"] or int(P0.get("tau_E1_free", 0)))
        self.o11 = v11_opts(P0, max(float(self.veto_T), 0.0))
        self.vw, self.vpk, self.mvd = self.o11["vw"], bool(self.o11["veto_peak"]), bool(self.o11["mvd_mode"])
        self.lp = int(self.o11["ecb_lp_mode"])
        self.k10 = not (self.lic == -1 and self.ecb == 0 and self.nom == 0 and not self.tef and not self.vw and not self.mvd
                        and not self.lp)
        if not self.k10:
            return
        lic, ecb, nom = self.lic, self.ecb, self.nom
        if not self.v5 or self.v5_mode == 3:
            raise ValueError("v10 options need v5_mode 1 or 2")
        if self.vw and ecb != 0:
            raise ValueError("v11 veto window / peak veto act on the spine eCB veto only (ecb_src 0)")
        if self.lp and ecb != 0:
            raise ValueError("v11 ecb_lp_mode acts on the spine eCB trigger only (ecb_src 0)")
        evk = {}                                                             # v11 single-bAP peak VDCC current (veto_peak)
        if lic in (1, 2, 3):
            self.gate_win = float(P0.get("gate_win", 100.0)) if lic in (1, 2) else 0.0
        else:
            self.gate_win = float(P0.get("lic_win", V10_DEFAULTS["lic_win"])) if lic >= 4 else 0.0
            if lic in (4, 5) and self.gate_win <= 0 and not self.llp:
                raise ValueError("v10 lic_src 4 / 5 need lic_win > 0 (or lic_lp 1)")
        self.gate_theta = float(P0.get("gate_theta", 0.0))
        thc = float(P0.get("gate_conv_thV", 3.5751446275588403))
        tauE1 = float(P0.get("tau_E1", 100.0)); isc = float(P0.get("i_scale", 1e-5))
        if self.veto_T <= 0:
            na = int(np.asarray(self.d["aptr"].copy_to_host())[-1])
            self.vq = cuda.to_device(np.zeros(max(na, 1), np.float64))
        eoff = np.asarray(self.d["eoff"].copy_to_host())
        N = self.n_ctrl
        shx = lic in (2, 5) or ecb > 0
        w1t = ecb == 1 and self.tef                                          # W1 table over TE1_GRID
        X = np.zeros(int(eoff[-1]) if shx else 1, np.uint16)
        oq = np.full(N, np.nan); self.burst_sh = np.full(N, np.nan)
        rests, nclip, ev, evr, evp, evw, evwr, keys, vqs = [], 0, {}, {}, {}, {}, {}, [], []
        for r in B.recs:
            t = np.asarray(r["t"], np.float64); sl = r["sl"]; n, T = r["effcai"].shape
            S = np.asarray(r["vdcc"], np.float32).astype(np.float64)
            if self.vpk:
                postv = np.sort(np.asarray(r["postspikes"], float).ravel())
                isov = np.array([not np.any((postv != x) & (postv > x - 100.0) & (postv < x + 20.0)) for x in postv], bool)
            if shx:
                if r.get("shaft_cai") is None:
                    raise ValueError(f"{r['pair']}/{r['proto']}: v10 shaft options need shaft_cai (BatchV2 signals)")
                sh = np.asarray(r["shaft_cai"], np.float64)
                rest = _shaft_rest(r, sh); rests.append(rest)
                du = (sh - rest[:, None]) * 1e3                              # uM above rest
                nm = np.rint(du * 1e3); nclip += int(np.sum(nm > 65535))
                post = np.sort(np.asarray(r["postspikes"], float).ravel())
                iso = np.array([not np.any((post != x) & (post > x - 100.0) & (post < x + 20.0)) for x in post], bool)
                b3 = "3ap_200hz" in r["proto"]
                hmr = np.append(np.diff(t), 0.0); aw = 1.0 - np.exp(-hmr / tauE1)   # shaft pool input weights
                awg = [1.0 - np.exp(-hmr / tg) for tg in TE1_GRID] if w1t else []
                C = np.concatenate([[0.0], np.cumsum(hmr)])[:T]                 # kernel clock at each step
            for ii in range(n):
                i = sl.start + ii
                e = np.asarray(r["effcai"][ii], np.float64)
                keys.append((r["pair"], int(r["syn"][ii])))
                if self.vpk:                                                 # peak step-mean VDCC current of isolated bAPs
                    tav = t[np.asarray(r["arr"][ii]).ravel().astype(np.int64)]
                    Lk = evk.setdefault(keys[-1], ([], []))
                    for x in postv[isov]:
                        a, b = np.searchsorted(t, [x - 1.0, x + 8.0])
                        if b > a:
                            Lk[0 if not np.any((tav > x - 150.0) & (tav < x + 8.0)) else 1].append(float(S[ii, a:b].max()))
                if shx:
                    X[eoff[i]:eoff[i] + T] = np.clip(nm[ii], 0, 65535).astype(np.uint16)
                    xr = np.clip(nm[ii], 0, 65535) * 1e-3                    # what the kernel reads, uM
                    ta = t[np.asarray(r["arr"][ii]).ravel().astype(np.int64)]
                    L = ev.setdefault(keys[-1], []); Lr = evr.setdefault(keys[-1], [])
                    Lw = evw.setdefault(keys[-1], []); Lwr = evwr.setdefault(keys[-1], [])
                    for x in post[iso]:                                      # tier 1 = gpu_v7_rho.bap_ref isolation
                        t1 = not np.any((ta > x - 150.0) & (ta < x + 8.0))
                        if not t1 and (lic != 5 and ecb != 1):
                            continue
                        a, b = np.searchsorted(t, [x - 1.0, x + 30.0])
                        if b > a:
                            (L if t1 else Lr).append(float(du[ii, a:b].max()))
                        Lx = Lw if t1 else Lwr
                        if ecb == 1 and len(Lx) < 5:
                            j = np.searchsorted(post, x, side="right")
                            tend = min(x + 300.0, post[j]) if j < len(post) else x + 300.0
                            b2 = np.searchsorted(t, tend)
                            if b2 > a:
                                Lx.append([float(_lp_peak(xr[a:b2], aw[a:b2]))]
                                          + [float(_lp_peak(xr[a:b2], ag[a:b2])) for ag in awg])
                    if lic == 5:
                        Lp = evp.setdefault(keys[-1], [])
                        for x in _iso_pre(ta, post):
                            a, b = np.searchsorted(t, [x - 1.0, x + 30.0])
                            if b > a:
                                Lp.append(float(du[ii, a:b].max()))
                    if ecb > 0:                                              # shaft veto input per arrival (veto_q order)
                        nz = np.unique(np.asarray(r["arr"][ii]).ravel().astype(np.int64))
                        if self.veto_T > 0:
                            cq = np.concatenate([[0.0], np.cumsum(aw * xr)])
                            je = np.searchsorted(C, C[nz] + self.veto_T, side="right")
                            vqs.append(cq[je] - cq[nz])
                        else:
                            vqs.append(np.zeros(len(nz)))
                    if b3:
                        self.burst_sh[i] = float(du[ii].max())
                else:
                    xr = np.zeros(1)
                if lic in (1, 2, 3):
                    v = _open_q(t, S[ii], e, xr, lic, self.gate_win, tauE1, isc, thc, float(TAU_EFFCA))
                    if v.size:
                        oq[i] = np.median(v)
            if shx:
                r["shaft_cai"] = None
        med = lambda D: np.array([np.median(D[k]) if D.get(k) else np.nan for k in keys])    # per physical synapse
        nw = 1 + (len(TE1_GRID) if w1t else 0)
        medw = lambda D: np.array([np.median(np.asarray(D[k]), axis=0) if D.get(k) else np.full(nw, np.nan) for k in keys])
        self.bap_sh = med(ev) if shx else np.full(N, np.nan)                 # v8 (tier 1 only): calibration print
        nan = np.full(N, np.nan)
        b1 = self.bap_sh; b2 = med(evr) if lic == 5 else nan
        self.bap5 = np.where(np.isfinite(b1), b1, b2) if lic == 5 else nan
        self.pre_sh = med(evp) if lic == 5 else nan
        if ecb == 1:
            w1, w2 = medw(evw), medw(evwr)
            wt = np.where(np.isfinite(w1[:, :1]), w1, w2)
            self.w1_sh = wt[:, 0]; self.w1g = wt[:, 1:] if w1t else None
            tier = np.where(np.isfinite(w1[:, 0]), 1, np.where(np.isfinite(w2[:, 0]), 2, 0))
        else:
            self.w1_sh, self.w1g = nan, None
            tier = np.where(np.isfinite(b1), 1, np.where(np.isfinite(b2), 2, 0)) if lic == 5 else np.zeros(N, int)
        self.tiers = tuple(int(np.sum(tier == k)) for k in (1, 2, 0))
        self._keys10 = keys
        self.open_q = oq
        self.X = cuda.to_device(X) if shx else cuda.to_device(np.zeros(1, np.float32))
        del X
        na = int(np.asarray(self.d["aptr"].copy_to_host())[-1])
        q = np.concatenate(vqs) if vqs else np.zeros(0)
        if ecb > 0:
            assert len(q) == na, (len(q), na)
        self.vqs = cuda.to_device(np.ascontiguousarray(q if len(q) else np.zeros(1), np.float64))
        self.cpr = cuda.to_device(np.ascontiguousarray(self.cpre, np.float64))
        self.cpo = cuda.to_device(np.ascontiguousarray(self.cpost, np.float64))
        self.ltg = cuda.to_device(np.log(TE1_GRID))
        self.tex = bool(self.t_exact)
        if not self.tex:
            self.fa = cuda.to_device(np.zeros(1, np.float64))
        if self.vpk:                                                         # tier 1, else tier 2, else NaN (pooled later)
            p1 = np.array([np.median(evk[k][0]) if k in evk and evk[k][0] else np.nan for k in keys])
            p2 = np.array([np.median(evk[k][1]) if k in evk and evk[k][1] else np.nan for k in keys])
            self.spk = np.where(np.isfinite(p1), p1, p2)
            self.spk_tiers = (int(np.isfinite(p1).sum()), int((~np.isfinite(p1) & np.isfinite(p2)).sum()),
                              int((~np.isfinite(self.spk)).sum()))
        else:
            self.spk = nan
        self._kern10 = _make_kernel_v10(weighted=self.v5_mode == 2, lic=lic, win=self.gate_win, tex=self.tex, ecb=ecb,
                                        nom=nom, cnt=self.cnt, llp=self.llp, tef=self.tef, tv=max(self.veto_T, 0.0),
                                        v11=self.o11)
        _INST.append(self)
        msg = (f"v10 lic_src {lic} (win {self.gate_win} ms, lic_lp {int(self.llp)}) ecb_src {ecb} no_mode {nom} "
               f"t_exact {int(self.tex)} tau_E1 free {int(self.tef)}")
        if self.vw or self.mvd or self.lp:
            msg += f"; v11 {self.o11}"
        if self.vpk:
            msg += f"; v11 single-bAP peak VDCC rows tier 1 / tier 2 / none {self.spk_tiers}"
        if lic in (1, 2, 3):
            msg += (f"; v8 gate_theta {self.gate_theta}: Vg-opening gate quantity per synapse q10/50/90 "
                    f"{np.nanquantile(oq, [.1, .5, .9]).round(4) if np.isfinite(oq).any() else '-'} "
                    f"{UNIT[lic]} ({int(np.isfinite(oq).sum())} of {N} synapses open at theta_V {thc:.4g})")
        if shx:
            rr = np.concatenate(rests) * 1e6
            msg += (f"; shaft rest q10/50/90 {np.quantile(rr, [.1, .5, .9]).round(1)} nM, {nclip} samples clipped at 65.5 uM; "
                    f"uint16 shaft {self.X.nbytes / 1e9:.2f} GB on GPU")
        if ecb > 0:
            msg += f"; shaft veto input per arrival q10/50/90 {np.quantile(q, [.1, .5, .9]).round(4)} uM ({na} arrivals)"
        if lic == 5 or ecb == 1:
            msg += f"; bAP shaft calibration rows tier 1 / tier 2 / none {self.tiers}"
        print(msg, flush=True)

    @staticmethod
    def _kappa(P0):
        if _STATE["kappa"] is None:
            q = np.concatenate([m.open_q for m in _INST]); q = q[np.isfinite(q)]
            thc = float(P0.get("gate_conv_thV", 3.5751446275588403))
            Q50 = float(np.median(q)) if q.size else float("nan")
            kap = float(P0.get("gate_kappa", 0.0)) or Q50 / thc
            _STATE.update(kappa=kap, Q50=Q50, theta_conv=thc, n_open=int(q.size))
            print(f"v8 kappa {kap:.6g} {UNIT.get(int(P0.get('gate_src', 0)), '')} per theta_V unit (Q50 {Q50:.6g} at theta_V "
                  f"{thc:.4g}, {q.size} synapses pooled over {len(_INST)} models)", flush=True)
            if v10_opts(P0)[0] == 2:
                _calib()
        return _STATE["kappa"]

    @staticmethod
    def _finalize():
        """Per-synapse v10 calibrations on the GPU, once, after every model is built: NaN -> pooled median (all models)."""
        if _STATE["fin"]:
            return
        _STATE["fin"] = True
        lic, ecb = _INST[0].lic, _INST[0].ecb
        fill = lambda name: float(np.nanmedian(np.concatenate([getattr(m, name) for m in _INST])))

        def fix(v, f, floor=0.0):
            return np.maximum(np.where(np.isfinite(v), v, f), floor)
        qd = lambda v: np.nanquantile(v, [.1, .5, .9]).round(4) if np.isfinite(v).any() else "-"
        cx = {}
        import fit_v6
        for G, name in fit_v6._M:
            cx[id(G)] = name
        if lic == 4:
            cdir = os.environ.get("CEXP_DIR", "")
            if not cdir:
                raise SystemExit("v10 lic_src 4 needs CEXP_DIR (cexp cpre_cai / cpost_cai)")
            for G, name in fit_v6._M:
                df = fit_v6.join_cexp(G, name, cdir)
                G.cpre_sp = np.maximum(df["cpre_cai"].to_numpy(float) - CA_REST_UM, 0.0)
                G.cpost_sp = np.maximum(df["cpost_cai"].to_numpy(float) - CA_REST_UM, 0.0)
        names = (("bap5", "pre_sh") if lic == 5 else ()) + (("w1_sh",) if ecb == 1 else ()) \
            + (("cpre_sp", "cpost_sp") if lic == 4 else ())
        fb = {k: fill(k) for k in names}
        if ecb == 1 and _INST[0].w1g is not None:
            gfb = np.nanmedian(np.concatenate([m.w1g for m in _INST]), axis=0)
        for m in _INST:
            one = np.zeros(m.n_ctrl)
            nn = {k: f"{int(np.sum(~np.isfinite(getattr(m, k))))}/{m.n_ctrl}" for k in names}
            if lic == 4:
                q1, q2 = fix(m.cpre_sp, fb["cpre_sp"]), fix(m.cpost_sp, fb["cpost_sp"])
            elif lic == 5:
                q1, q2 = fix(m.pre_sh, fb["pre_sh"]), fix(m.bap5, fb["bap5"])
            elif lic == 6:
                q1, q2 = np.asarray(m.cpre, np.float64), np.asarray(m.cpost, np.float64)
            else:
                q1 = q2 = one
            w1 = fix(m.w1_sh, fb["w1_sh"], 1e-6) if ecb == 1 else np.ones(m.n_ctrl)
            if ecb == 1 and m.w1g is not None:
                wg = np.maximum(np.where(np.isfinite(m.w1g), m.w1g, gfb[None, :]), 1e-6)
            else:
                wg = np.ones((1, len(TE1_GRID)))
            m.q10 = (q1, q2, w1)
            m.d10 = tuple(cuda.to_device(np.ascontiguousarray(v, np.float64)) for v in (q1, q2, w1))
            m.w1gd = cuda.to_device(np.ascontiguousarray(wg, np.float64))
            print(f"v10 calib {cx.get(id(m), '')} n {m.n_ctrl}: licence q1 {qd(q1)} q2 {qd(q2)}"
                  f"{' (Cpre_sp, Cpost_sp uM)' if lic == 4 else ' (presh, bapsh uM)' if lic == 5 else ' (c_pre, c_post)' if lic == 6 else ''}"
                  f" | W1_sh {qd(w1) if ecb == 1 else '-'} uM | NaN rows filled {nn} with pooled medians {fb}", flush=True)
        # v11: single-bAP peak VDCC reference (veto_peak), and the P1_i (ecb_ref 2, cexp vdcc_q_post) NaN fraction per pathway
        vpk = any(m.vpk for m in _INST)
        spf = fill("spk") if vpk else 1.0
        for m in _INST:
            if vpk:
                sp = fix(m.spk, spf, 1e-12)
                m.spkd = cuda.to_device(np.ascontiguousarray(sp, np.float64))
            else:
                sp = None
                m.spkd = cuda.to_device(np.ones(1, np.float64))
            try:
                uE = np.asarray(m.u6[1].copy_to_host(), np.float64)
                nanp1 = "-"
                cdir = os.environ.get("CEXP_DIR", "")
                if cdir:
                    df = fit_v6.join_cexp(m, cx.get(id(m), ""), cdir)
                    vq = df["vdcc_q_post"].to_numpy(float)
                    nanp1 = f"{int(np.sum(~np.isfinite(vq)))}/{len(vq)} ({np.mean(~np.isfinite(vq)):.3f})"
                print(f"v11 calib {cx.get(id(m), '')}: P1_i (uE) q10/50/90 {qd(uE)}, cexp vdcc_q_post NaN rows {nanp1} "
                      f"(pooled-median filled)" + (f" | single-bAP peak VDCC I1_i q10/50/90 {qd(sp)} (tier 1/2/none "
                      f"{m.spk_tiers}, NaN -> pooled {spf:.4g}), I1/P1 q10/50/90 {qd(sp / np.maximum(uE, 1e-12))}"
                      if sp is not None else ""), flush=True)
            except Exception as ex:                                          # diagnostics only
                print(f"v11 calib print failed: {type(ex).__name__}: {ex}", flush=True)

    def _params10(self, Ps):
        rows = []
        g = lambda p, k: float(p[k] if k in p else self.f10.get(k, v10_default(k, self.nom)))
        for p in Ps:
            if self.lic in (4, 6):
                lc = (g(p, "a20"), g(p, "a21"))
            elif self.lic == 5:
                lc = (g(p, "kappa_pre"), g(p, "kappa_post"))
            else:
                lc = (0.0, 0.0)
            rows.append([*lc, g(p, "ecb_theta_uM"), g(p, "tau_NO"), g(p, "A_NO"), g(p, "dNO_step"), g(p, "d_NO_max"),
                         g(p, "a_NO_pre"), g(p, "a_NO_post"), g(p, "tau_L"),
                         g(p, "theta_MVD_lo"), g(p, "theta_MVD_hi"), g(p, "veto_peak_k"), g(p, "theta_MVD_abs")])     # v11 slots 20-23
        return np.ascontiguousarray(np.concatenate([self._params5(Ps), np.array(rows, np.float64)], axis=1))

    def rho_dpre_syn(self, td, tp, Ps):
        if not self.k10:
            return super().rho_dpre_syn(td, tp, Ps)
        kap = self._kappa(self.P0) if self.lic in (1, 2, 3) else 0.0
        self._finalize()
        P, N = td.shape
        out_rho = cuda.device_array((P, N), np.float64); out_d = cuda.device_array((P, N, 3), np.float64)
        ocnt = cuda.device_array((P, N, NCNT) if self.cnt else (1, 1, NCNT), np.float64)
        bx, by = self.block
        g = self.d
        self._kern10[(-(-P // bx), -(-N // by)), (bx, by)](
            g["E"], g["S"], g["Hg"], g["eoff"], g["slen"], g["hoff"], g["aptr"], g["astep"], g["acnt"], g["rho0"],
            g["order"], cuda.to_device(np.ascontiguousarray(td, np.float64)),
            cuda.to_device(np.ascontiguousarray(tp, np.float64)), cuda.to_device(self._params10(Ps)), out_rho, out_d,
            *self.u6, self.vq, self.X, float(kap), float(self.gate_theta), self.fa, *self.d10, self.vqs, self.cpr,
            self.cpo, ocnt, self.w1gd, self.ltg, self.spkd)
        if self.cnt and P == 1:
            c = ocnt.copy_to_host()[0]
            tot = c.sum(0); ns = (c > 0).sum(0)
            print(f"v10 counts ({N} synapse rows): licensed pot steps {tot[0]:.0f} (rows {ns[0]}), blocked pot steps "
                  f"{tot[1]:.0f} (rows {ns[1]}), eCB steps {tot[2]:.0f} (rows {ns[2]}), vetoed eCB triggers {tot[3]:.0f} "
                  f"(rows {ns[3]}), NO steps {tot[4]:.0f} (rows {ns[4]})"
                  + (f", MVD-depressing steps {tot[5]:.0f} (rows {ns[5]})" if self.mvd else ""), flush=True)
        return out_rho.copy_to_host(), out_d.copy_to_host()


def _calib():
    """GATE_ALT checklist: single-bAP shaft Delta[Ca] by cexp location, proximal median, Spearman vs vdcc_q_post."""
    try:
        import fit_v6
        from scipy.stats import spearmanr
        cdir = os.environ.get("CEXP_DIR", "/scratch/dhuruva/split1/cexp")
        allb, allq = [], []
        for G, name in fit_v6._M:
            if G not in _INST:
                continue
            df = fit_v6.join_cexp(G, name, cdir)
            u = np.unique([f"{p}:{s}" for p, s in G._keys], return_index=True)[1]     # one row per physical synapse
            b = G.bap_sh[u]; loc = df["loc"].astype(str).to_numpy()[u]; dist = df["dist_um"].to_numpy(float)[u]
            qp = df["vdcc_q_post"].to_numpy(float)[u]
            q = lambda v: np.nanquantile(v, [.1, .5, .9]).round(4) if np.isfinite(v).any() else "-"
            parts = [f"{L} n{int(np.isfinite(b[loc == L]).sum())} {q(b[loc == L])}" for L in sorted(set(loc))]
            prox = b[dist <= 100.0]
            ok = np.isfinite(b) & np.isfinite(qp)
            rs = spearmanr(b[ok], qp[ok])[0] if ok.sum() > 2 else float("nan")
            bu = G.burst_sh
            print(f"v8 calib {name}: single-bAP peak shaft dCa uM q10/50/90 {q(b)} ({int(np.isfinite(b).sum())} of {len(b)} "
                  f"synapses) | by loc {'; '.join(parts)} | proximal <=100 um median "
                  f"{np.nanmedian(prox) if np.isfinite(prox).any() else float('nan'):.4f} (n {int(np.isfinite(prox).sum())}) | "
                  f"Spearman vs vdcc_q_post {rs:.3f} (n {int(ok.sum())}) | 3AP 200 Hz record peak q10/50/90 {q(bu)} "
                  f"(n {int(np.isfinite(bu).sum())})", flush=True)
            allb.append(b[ok]); allq.append(qp[ok])
        if allb:
            b = np.concatenate(allb); qp = np.concatenate(allq)
            print(f"v8 calib pooled: Spearman single-bAP shaft dCa vs vdcc_q_post {spearmanr(b, qp)[0]:.3f} (n {len(b)})",
                  flush=True)
    except Exception as ex:             # diagnostics only, never stops a fit
        print(f"v8 calib failed: {type(ex).__name__}: {ex}", flush=True)


if __name__ == "__main__":
    _s = json.loads(sys.argv[sys.argv.index("--set") + 1]) if "--set" in sys.argv else {}
    if int(_s.get("ecb_ref", 0)) == 2:                  # as gpu_v7_rho.__main__
        os.environ.setdefault("CEXP_DIR", "/scratch/dhuruva/split1/cexp")
        _e = os.environ.get("SCALE_E") or '{"vdcc_q_post": 100000.0}'
        assert json.loads(_e) == {"vdcc_q_post": 1e5}, f"ecb_ref 2 needs SCALE_E vdcc_q_post x 1e5, got {_e}"
        os.environ["SCALE_E"] = _e
    elif int(_s.get("ecb_ref", 0)) == 1:
        assert not os.environ.get("SCALE_E"), "ecb_ref 1 is not combined with SCALE_E"
    src = int(_s.get("gate_src", 0))
    LIC, ECB, NOM = v10_opts(_s)
    LLP = int(_s.get("lic_lp", 0))
    _free = sys.argv[sys.argv.index("--free-filters") + 1].split(",") if "--free-filters" in sys.argv else []
    TEF = "tau_E1" in _free
    _STATE["te1_free"] = TEF
    O11 = v11_opts(_s, max(float(_s.get("veto_T", 0.0)), 0.0))
    MVD, VPK = int(O11["mvd_mode"]), int(O11["veto_peak"])
    if MVD and int(_s.get("ecb_ref", 0)) != 2:
        raise SystemExit("v11 mvd_mode 1 needs ecb_ref 2 (P1_i = 1e5 x cexp vdcc_q_post)")
    V11 = bool(MVD or O11["vw"] or O11["ecb_lp_mode"])
    new = _s.get("lic_src") is not None or ECB > 0 or NOM > 0 or TEF or V11
    KEYS = v10_keys(LIC, ECB, NOM, LLP, MVD, VPK, int(O11["mvd_ref"]))
    import fit_v6
    # v10 parameters freeable through --free-filters (boxes in the fit_v2 FILTER_BOX convention; BOX env still overrides)
    fit_v6.F.FILTER_BOX.update(V10_BOX)
    if NOM == 3:
        fit_v6.F.FILTER_BOX.update(V10_BOX_NO3)
    fit_v6.F.LOG_FILTERS |= V10_LOG
    for _k, _v in V10_DEFAULTS.items():
        fit_v6.MV.DEFAULTS.setdefault(_k, v10_default(_k, NOM))
    if NOM:                                             # MV tau_NO 5 ms is the v2 NO filter; v10 N has its own default
        fit_v6.MV.DEFAULTS["tau_NO"] = v10_default("tau_NO", NOM)
    if KEYS:                                            # fit_v2.unpack forces A_NO = 0 (log slot -inf): restore SET / default
        _u5 = fit_v6.unpack5

        def _unpack10(x):
            a, P = _u5(x)
            for k in KEYS:
                if k not in fit_v6._CFG["free_filters"]:
                    P[k] = float(_s.get(k, v10_default(k, NOM)))
            if int(_s.get("lic_tie", 0)):               # default off: one licence coefficient (a21 = a20 | kappa_pre = kappa_post)
                if LIC in (4, 6):
                    P["a21"] = P["a20"]
                elif LIC == 5:
                    P["kappa_pre"] = P["kappa_post"]
            return a, P
        fit_v6.unpack5 = _unpack10
    TE1_START = float(_s.get("tau_E1_start", 50.0))
    FK = [k for k in KEYS if k in _free]
    if FK or (TEF and TE1_START > 0):
        # DE start of every seed member: a free v10 parameter given in SET starts there (seed pre otherwise; a v2 seed's
        # A_NO 0 / non-positive log value -> v10 default); free tau_E1 starts at tau_E1_start (0 = the seed's value)
        _pack = fit_v6.F.pack

        def _pack10(a, P, cfg):
            P = dict(P)
            for k in FK:
                if k in _s:
                    P[k] = float(_s[k])
                elif k not in P or (k in V10_LOG and not float(P[k]) > 0):
                    P[k] = v10_default(k, NOM)
            if TEF and TE1_START > 0:
                P["tau_E1"] = TE1_START
            return _pack(a, P, cfg)
        fit_v6.F.pack = _pack10
    fit_v6.gpu_v6_rho.GPUModelV6 = GPUModelV11
    if LIC in (2, 5) or ECB > 0:                        # keep the shaft trace in the records (popped after the build)
        _BV = fit_v6.BatchV2
        fit_v6.BatchV2 = lambda *a, signals=("vdcc",), **k: _BV(*a, signals=tuple(signals) + ("shaft_cai",), **k)
    save = sys.argv[sys.argv.index("--save") + 1]
    code = 0
    try:
        fit_v6.main()
    except SystemExit as ex:
        if (src == 0 and not new) or str(ex) != "REPRO DIFF":
            raise
        print("REPRO DIFF expected (v8 gate_src != 0 or a v10 option: a different rule than the seed)", flush=True)
    if (src or new) and os.path.isfile(save + ".json"):
        out = json.load(open(save + ".json"))
        if LIC in (1, 2, 3):
            thV = float(out["v5"]["theta_V"]); gt = float(_s.get("gate_theta", 0.0))
            out["v8"] = dict(gate_src=LIC, t_exact=int(_s.get("t_exact", 0)), gate_win=float(_s.get("gate_win", 100.0)) if LIC in (1, 2) else 0.0, gate_theta=gt,
                             kappa=_STATE.get("kappa"), Q50=_STATE.get("Q50"), theta_conv=_STATE.get("theta_conv"),
                             n_open=_STATE.get("n_open"), theta_G=gt if gt > 0 else thV * float(_STATE.get("kappa") or np.nan),
                             unit=UNIT[LIC])
            print(f"v8: theta_G {out['v8']['theta_G']:.6g} {UNIT[LIC]} -> {save}.json", flush=True)
        if new:
            out["v10"] = dict(lic_src=LIC, ecb_src=ECB, no_mode=NOM, lic_win=float(_s.get("lic_win", V10_DEFAULTS["lic_win"])),
                              lic_lp=LLP, tau_E1_free=int(TEF), tau_E1=out["pre"].get("tau_E1"),
                              params={k: out["pre"].get(k) for k in KEYS}, ca_rest_uM=CA_REST_UM)
            if LIC in (4, 6):                           # the licence coefficients in the json a block (a20, a21 slots)
                out["a"]["a20"] = out["pre"]["a20"]; out["a"]["a21"] = out["pre"]["a21"]
            print(f"v10: {out['v10']} -> {save}.json", flush=True)
            if V11:
                out["v11"] = dict(O11, params={k: out["pre"].get(k) for k in KEYS if k in V11_KEYS})
                print(f"v11: {out['v11']} -> {save}.json", flush=True)
        json.dump(out, open(save + ".json", "w"), indent=1)
    sys.exit(code)
