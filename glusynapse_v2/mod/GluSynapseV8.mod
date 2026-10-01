COMMENT
/**
 * @file GluSynapse.mod
 * Modifications by Dhuruva - Modified Rho Dynamics, added uncaging mode
 * v2 (2026-09-27): presynaptic plasticity pathways (mGluR/eCB pre-LTD, NO pre-LTP),
 *     see plastyfire/glusynapse_v2/DESIGN.md. All v2 amplitudes default to 0,
 *     in which case the mechanism behaves exactly as GluSynapse.
 * V4 (2026-10-01, copy of GluSynapse_v2.mod, live form of the offline rule gpu_v4_rho / rho_v4 / scan_vgate_amp):
 *     t_drive 4: eCB event = upward crossing of -ica_VDCC through K_ca_GB (re-armed below K_ca_GB*ca_hyst_GB),
 *         S += 1 - min(bglu, 1), bglu = own pre arrivals, +1 each, decay tau_glu_GB (model_v2.ca_events/glu_weight).
 *     vamp_mode 1 (C1) / 2 (C2): potentiation also needs Vg > theta_VA_GB, Vg' = -Vg/tau_E1 + (-ica_VDCC)/i_scale.
 *         C2 also drops depression above theta_p. 0 = off (v2 rule).
 *     express_x 1 (validation only): expression (gmax_AMPA, Use) follows a shadow v1 rule rhox with its own
 *         thresholds theta_dx/px and rates gamma_dx/px, without dpre, so the run reproduces the v1 prefire traces
 *         while rho_GB / dpre_GB run passively. 0 = expression follows rho_GB and dpre_GB (v2).
 *     bin_GB > 0 (2026-10-01): pre path on the offline grid (PROCEDURE bin_step), N driven by the bin mean of
 *         -ica_VDCC (extract_v2.step_mean). The continuous path (bin_GB 0) rectifies the instantaneous current and
 *         integrates exactly; it gave dpre +0.01..0.02 above offline in equiv 22131696.
 *     All new options default to off, in which case the mechanism behaves as GluSynapse_v2.
 * V5 (2026-10-01, copy of GluSynapseV4.mod, live form of the offline kernel rho_redesign/gpu_v5_rho.py,
 *     rho_redesign/V5_DESIGN.md). v5_mode_GB 0 (default) = GluSynapseV4 unchanged. The v5 rule is selected by the
 *     globals, with the v4 eCB/NO chains switched off (A_mglu_GB = A_NO_GB = 0):
 *     rho:    Chindemi thresholds theta_d/p_GB per synapse (rho_gamma 1), C1 gate (vamp_mode_GB 1,
 *             theta_VA_GB = theta_V): pot = [effcai > theta_p][Vg > theta_V], dep = [effcai > theta_d], a failed
 *             gate counts as depression. Vg' = -Vg/tau_E1 + (-ica_VDCC)/i_scale (the v4 C1 pool, unchanged).
 *     W:      v5_mode 2 (v5b / v5c): W' = -W/tau_E1 + (1 - min(bglu, 1)) (-ica_VDCC)/i_scale, bglu = own pre
 *             arrivals (+1 each, released or not), decay tau_glu_GB = tau_d_NMDA 70 ms (as v4 t_drive 4).
 *     eCB:    at each own pre arrival (NET_RECEIVE flag 0) with trigger > theta_eCB_GB (theta_eCB_GB <= 0: tied to
 *             theta_VA_GB), dpre <- dpre_min + (dpre - dpre_min)(1 - A_eCB). Trigger = W (v5_mode 2) or Vg
 *             (v5_mode 1, plain v5). A_eCB_GB 1: one gated pre spike sets dpre = dpre_min. The check reads the
 *             trigger before this arrival's bglu increment (the kernel's order; W does not jump at an arrival).
 *             necb_GB counts the steps (diagnostic). No NO arm. Readout unchanged: Use = min(1, U0 (1 + dpre)).
 *             v5_mode 3 (E2) is not implemented here.
 *     Discretisation (deliberate): the offline kernel snaps every own arrival to the first sample of the extracted
 *     grid at or after it (model_v2.arrival_index, 0.25 ms in the spike windows) and advances W, Vg, bglu with an
 *     exact step filter on the sampled -ica_VDCC. Here the arrival, the eCB check and the bglu jump happen at the
 *     exact arrival time and W, Vg, bglu are integrated continuously (CVODE in live runs). The eCB step is applied
 *     at the arrival also when bin_GB > 0 (the bin mode only defers the v4 mGluR step).
 * V7 (2026-10-01, copy of GluSynapseV5.mod, live form of rho_redesign/gpu_v7_rho.py; V5 unchanged at veto_T_GB 0,
 *     uE_GB 1). Two additions, both synapse-local:
 *     uE_GB (RANGE, default 1): per-synapse threshold scale, theta_eCB,i = theta_eCB_GB uE_GB (gpu_v6/v7 thE0 uE[i]);
 *             theta_eCB_GB <= 0 still means theta_VA_GB (unscaled). ecb_ref 2: uE_GB = 1e5 vdcc_q_post (set by
 *             bcl_validation/live_v7.py from the A4 cexp tables).
 *     veto_T_GB > 0 (ms), bAP veto: an own arrival whose trigger fires (W or Vg > theta_eCB,i, read at the arrival
 *             as V5) takes its eCB step only if the own unweighted VDCC charge in the next veto_T ms,
 *             Q = qca(t_a + veto_T) - qca(t_a), qca' = -ica_VDCC/i_scale (the Vg drive, no decay), does not exceed
 *             the same theta_eCB,i (kernel: trig and not (Q > thE)). The window lies after the arrival, so the
 *             decision is a self-event at t_a + veto_T (flag 1000 + qca(t_a); flags 0-21 are V5's), where the step
 *             is applied or counted in nveto_GB. The step is a contraction towards dpre_min, so applying it
 *             veto_T later gives the same final dpre (dpre has no other drive with A_mglu = A_NO = 0); the 25 ms
 *             delay only shifts Use during that window when expression follows dpre (express_x_GB 0).
 *             nvpend_GB = decisions still pending (must be 0 at the end of a run: arrivals within veto_T of
 *             tstop are left undecided). Uses qca_GB, so the veto needs bin_GB 0 (the bin mode resets qca_GB);
 *             live_v7 always runs bin_GB 0. Offline grid: the kernel sums bV s from the arrival's grid sample to
 *             the first sample past t_a + veto_T (<= 0.25 ms longer than the live window in the spike windows).
 * V8 (2026-10-01, copy of GluSynapseV7.mod, live form of the G1 shaft gate of rho_redesign/gpu_v8_rho.py, gate_src 2,
 *     gate_win 100; GATE_IMPL.md / GATE_ALT.md). gate_src_GB 0 (default) = GluSynapseV7 unchanged (no new STATE, the
 *     new WATCHes are not armed, so the CVODE system is V7's).
 *     gate_src_GB 2: the LTP licence is [max over [t - gate_win, t] of dsh > theta_sh_GB], dsh = cai - cai_rest_GB,
 *             cai = Ca of the synapse's own segment (cad pool; GluSynapse never writes ica, so it carries no synaptic
 *             NMDA/VDCC Ca). It REPLACES the Vg > theta_VA licence (gv_GB is not used for pot): potE = pot gsh, and a
 *             failed licence counts as depression (dep (1 - potE) gamma_d), as the kernel (pot = 0 when not lic).
 *             theta_sh_GB <= 0: no licence (pot ungated), as the kernel's thG > 0 condition.
 *             Windowed max: WATCH flags 22 (dsh rises above theta_sh: gsh = 1, ash = 1) and 23 (falls below: ash = 0,
 *             tlsh = t, self-event flag 24 at t + gate_win); flag 24 closes (gsh = 0) only if dsh is still below and
 *             no later fall moved tlsh (t - tlsh >= gate_win), so the licence is open from the first crossing until
 *             gate_win after the last fall. Kernel: lic_k = [t_k - t_last <= gate_win], t_last = last sample j < k with
 *             x_j > theta_G (x = rint(nM) of the 0.25 ms grid shaft trace), i.e. it opens one sample later and closes
 *             up to one sample earlier than live (<= 0.25 ms each in the spike windows; threshold quantisation 0.5 nM).
 *     Rest (kernel _shaft_rest: shaft cai at the last grid sample before min(pre, post spikes) - 1 ms):
 *             t_rest_GB >= 0 (GLOBAL, ms): cai_rest_GB follows cai while t <= t_rest_GB (INITIAL and BREAKPOINT) and is
 *             frozen after it; live_v8 sets t_rest_GB = first pre / post spike of the record - 1 ms.
 *             t_rest_GB < 0: cai_rest_GB (RANGE) is whatever the user set.
 *     Diagnostics (RANGE): nsh_GB = licence openings, tsh_GB = closed licence time (ms; add t - ton_GB if gsh_GB = 1 at
 *             the end), gsh_GB = licence state. gate_src 1 / 3 (spine Ca, c*) and gate_win 0 (G2) are not implemented.
 * @brief Probabilistic synapse featuring long-term plasticity
 * @author king, chindemi, rossert
 * @date 2021-05-19
 * @version 1.0.1
 * @remark Copyright 2005-2023 Blue Brain Project / EPFL
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * 
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 Glutamatergic synapse model featuring:
1) AMPA receptor with a dual-exponential conductance profile.
2) NMDA receptor  with a dual-exponential conductance profile and magnesium
   block as described in Jahr and Stevens 1990.
3) Tsodyks-Markram presynaptic short-term plasticity as Barros et al. 2019.
   Implementation based on the work of Eilif Muller, Michael Reimann and
   Srikanth Ramaswamy (Blue Brain Project, August 2011), who introduced the
   2-state Markov model of vesicle release. The new model is an extension of
   Fuhrmann et al. 2002, motivated by the following constraints:
        a) No consumption on failure
        b) No release until recovery
        c) Same ensemble averaged trace as canonical Tsodyks-Markram using same
           parameters determined from experiment.
   For a pre-synaptic spike or external spontaneous release trigger event, the
   synapse will only release if it is in the recovered state, and with
   probability u (which follows facilitation dynamics). If it releases, it will
   transition to the unrecovered state. Recovery is as a Poisson process with
   rate 1/Dep.
   John Rahmon and Giuseppe Chindemi introduced multi-vesicular release as an
   extension of the 2-state Markov model of vesicle release described above
   (Blue Brain Project, February 2017).
4) NMDAR-mediated calcium current. Fractional calcium current Pf_NMDA from
   Schneggenburger et al. 1993. Fractional NMDAR conductance treated as a
   calcium-only permeable channel with Erev = 40 mV independent of extracellular
   calcium concentration (see Jahr and Stevens 1993). Implemented by Christian
   Rossert and Giuseppe Chindemi (Blue Brain Project, 2016).
5) Spine volume.
6) VDCC.
7) Postsynaptic calcium dynamics.
8) Long-term synaptic plasticity. Calcium-based STDP model based on Graupner and
   Brunel 2012.
Model implementation, optimization and simulation curated by James King (Blue
Brain Project, 2021).
ENDCOMMENT


TITLE Glutamatergic synapse

NEURON {
    THREADSAFE
    POINT_PROCESS GluSynapseV8
    : AMPA Receptor
    GLOBAL tau_r_AMPA, E_AMPA
    RANGE tau_d_AMPA, gmax0_AMPA, gmax_d_AMPA, gmax_p_AMPA, g_AMPA, calcium_current_flag
    : NMDA Receptor
    GLOBAL scale_NMDA, slope_NMDA
    GLOBAL tau_r_NMDA, tau_d_NMDA, E_NMDA
    RANGE gmax_NMDA, g_NMDA, i_AMPA, i_NMDA
    RANGE mg
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    RANGE Use, Dep, Fac, Nrrp
    RANGE Use_d, Use_p
    BBCOREPOINTER rng
    : NMDAR-mediated calcium current
    RANGE ica_NMDA
    : Spine
    RANGE volume_CR, cai_CR, effcai_GB
    : VDCC (R-type)
    GLOBAL ljp_VDCC, vhm_VDCC, km_VDCC, mtau_VDCC, vhh_VDCC, kh_VDCC, htau_VDCC, gca_bar_VDCC
    RANGE ica_VDCC
    : Postsynaptic Ca2+ dynamics
    GLOBAL gamma_ca_CR, tau_ca_CR, min_ca_CR, cao_CR
    : Long-term synaptic plasticity
    GLOBAL rho_star_GB, tau_ind_GB, tau_exp_GB, tau_effca_GB
    GLOBAL gamma_d_GB, gamma_p_GB
    RANGE theta_d_GB, theta_p_GB, rho0_GB, dep_GB, pot_GB
    : Misc
    RANGE vsyn, synapseID, selected_for_report, verbose
    NONSPECIFIC_CURRENT i
    USEION ca READ cai
    : Shaft calcium (from CaDynamics.cad on the parent section)
    RANGE shaft_cai
    : v2 presynaptic plasticity (DESIGN.md sec 3)
    GLOBAL ca_sh_rest_GB, ca_scale_GB
    GLOBAL theta_T_GB, tau_T_GB, A_mglu_GB
    GLOBAL theta_NO_GB, tau_NO_GB, tau_Z_GB, theta_Z_GB, A_NO_GB, z_het_GB
    GLOBAL dpre_min_GB, dpre_max_GB
    GLOBAL t_drive_GB, theta_Te_GB, Te_scale_GB, tau_E1_GB, theta_Tg_GB, theta_V_GB, w_V_GB, tau_b_GB, hyst_V_GB
    GLOBAL pre_drive_GB, i_scale_GB, theta_Ti_GB, theta_NOi_GB
    GLOBAL no_drive_GB, theta_NOe_GB, e_scale_GB, theta_NOc_GB, c_scale_GB, theta_N_GB
    RANGE dpre0_GB, Use_tgt_GB
    : V4: t_drive 4, C1/C2 VDCC gate, shadow expression
    GLOBAL K_ca_GB, ca_hyst_GB, tau_glu_GB
    GLOBAL vamp_mode_GB, theta_VA_GB
    GLOBAL express_x_GB, gamma_dx_GB, gamma_px_GB
    RANGE theta_dx_GB, theta_px_GB, nev_GB, gv_GB, potx_GB, depx_GB
    : V4 bin mode (offline-grid pre path)
    GLOBAL bin_GB, bin_t0_GB, bin_t1_GB
    RANGE Nd_GB, Td_GB, npend_GB, nbin_GB, tTsum_GB
    : V5 shared-pool eCB-LTD
    GLOBAL v5_mode_GB, A_eCB_GB, theta_eCB_GB
    RANGE necb_GB
    : V7 bAP veto and per-synapse eCB threshold scale
    GLOBAL veto_T_GB
    RANGE uE_GB, nveto_GB, nvpend_GB
    : V8 G1 shaft gate
    GLOBAL gate_src_GB, gate_win_GB, theta_sh_GB, t_rest_GB
    RANGE cai_rest_GB, gsh_GB, nsh_GB, tsh_GB, ton_GB
    RANGE conductance
    RANGE next_delay
    BBCOREPOINTER delay_times, delay_weights
    GLOBAL nc_type_param
    GLOBAL minis_single_vesicle
    GLOBAL init_depleted
    : For debugging
    :RANGE sgid, tgid
    RANGE uncaging_mode
}


UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (uS)    = (microsiemens)
    (nS)    = (nanosiemens)
    (pS)    = (picosiemens)
    (umho)  = (micromho)
    (um)    = (micrometers)
    (mM)    = (milli/liter)
    (uM)    = (micro/liter)
    FARADAY = (faraday) (coulomb)
    PI      = (pi)      (1)
    R       = (k-mole)  (joule/degC)
}


PARAMETER {
    celsius                     (degC)
    : AMPA Receptor
    tau_r_AMPA      = 0.2       (ms)        : Tau rise, dual-exponential conductance profile
    tau_d_AMPA      = 1.7       (ms)        : Tau decay, IMPORTANT: tau_r < tau_d
    E_AMPA          = 0         (mV)        : Reversal potential
    gmax0_AMPA      = 1.0       (nS)        : Initial peak conductance
    gmax_d_AMPA     = 1.0       (nS)        : Peak conductance in the depressed state
    gmax_p_AMPA     = 2.0       (nS)        : Peak conductance in the potentitated state
    : NMDA Receptor
    mg              = 1         (mM)        : Extracellular magnesium concentration
    scale_NMDA      = 2.552     (mM)        : Scale of the mg block (Vargas-Caballero and Robinson 2003)
    slope_NMDA      = 0.072     (/mV)       : Slope of the ma block (Vargas-Caballero and Robinson 2003)
    tau_r_NMDA      = 0.29      (ms)        : Tau rise, dual-exponential conductance profile
    tau_d_NMDA      = 70        (ms)        : Tau decay, IMPORTANT: tau_r < tau_d
    E_NMDA          = -3        (mV)        : Reversal potential (Vargas-Caballero and Robinson 2003)
    gmax_NMDA       = 0.55      (nS)        : Peak conductance
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use             = 0.5       (1)         : Initial utilization of synaptic efficacy
    Dep             = 100       (ms)        : Relaxation time constant from depression
    Fac             = 10        (ms)        : Relaxation time constant from facilitation
    Nrrp            = 1         (1)         : Number of release sites for given contact
    Use_d           = 0.2       (1)         : Depressed Use
    Use_p           = 0.8       (1)         : Potentiated Use
    : Spine
    volume_CR       = 0.087     (um3)       : From spine data by Ruth Benavides-Piccione (unpublished)
    : VDCC (R-type)
    gca_bar_VDCC    = 0.0744    (nS/um2)    : Density spines: 20 um-2 (Sabatini 2000), unitary conductance VGCC 3.72 pS (Bartol 2015)
    ljp_VDCC        = 0         (mV)
    vhm_VDCC        = -5.9      (mV)        : v 1/2 for act, Magee and Johnston 1995 (corrected for m*m)
    km_VDCC         = 9.5       (mV)        : act slope, Magee and Johnston 1995 (corrected for m*m)
    vhh_VDCC        = -39       (mV)        : v 1/2 for inact, Magee and Johnston 1995
    kh_VDCC         = -9.2      (mV)        : inact, Magee and Johnston 1995
    mtau_VDCC       = 1         (ms)        : max time constant (guess)
    htau_VDCC       = 27        (ms)        : max time constant 100*0.27
    : Postsynaptic Ca2+ dynamics
    gamma_ca_CR     = 0.04      (1)         : Percent of free calcium (not buffered), Sabatini et al 2002: kappa_e = 24+-11 (also 14 (2-31) or 22 (18-33))
    tau_ca_CR       = 12        (ms)        : Rate of removal of calcium, Sabatini et al 2002: 14ms (12-20ms)
    min_ca_CR       = 70e-6     (mM)        : Sabatini et al 2002: 70+-29 nM, per AP: 1.1 (0.6-8.2) uM = 1100 e-6 mM = 1100 nM
    cao_CR          = 2.0       (mM)        : Extracellular calcium concentration in slices
    : Long-term synaptic plasticity
    rho_star_GB     = 0.5       (1)
    tau_ind_GB      = 70        (s)
    tau_exp_GB      = 100       (s)
    tau_effca_GB    = 200       (ms)
    gamma_d_GB      = 100       (1)
    gamma_p_GB      = 450       (1)
    theta_d_GB      = 0.006     (us/liter)
    theta_p_GB      = 0.001     (us/liter)
    rho0_GB         = 0         (1)
    : Misc
    synapseID       = 0
    verbose         = 0
    selected_for_report = 0
    conductance     = 0.0
    nc_type_param = 5
    minis_single_vesicle = 0   :// 0 -> no limit (old behavior)
    init_depleted = 0          :// 0 -> init full (old behavior)
    :sgid = -1
    :tgid = -1
    calcium_current_flag = 1

    uncaging_mode       = 0    : 0=normal synapse, 1=uncaging mode
    : v2 presynaptic plasticity. A_mglu_GB = A_NO_GB = 0 -> identical to GluSynapse
    ca_sh_rest_GB   = 6.5e-5    (mM)        : resting shaft Ca (CaDynamics minCai)
    ca_scale_GB     = 1e-3      (mM)        : normalisation of the shaft Ca drive
    theta_T_GB      = 1e-4      (mM)        : VGCC-Ca threshold of the pre-LTD (PLC/eCB) trace
    tau_T_GB        = 50        (ms)        : pre-LTD trace time constant
    A_mglu_GB       = 0         (1)         : pre-LTD step per presynaptic spike (mGluR x T)
    theta_NO_GB     = 5e-4      (mM)        : dendritic Ca-spike threshold for NO synthesis
    tau_NO_GB       = 5         (ms)        : NO trace time constant (tau_Z > tau_NO keeps NO causal)
    tau_Z_GB        = 30        (ms)        : presynaptic activity trace (NO coincidence window)
    theta_Z_GB      = 0         (1)         : pre-rate threshold of NO action, tanh(pos(Z - theta_Z)) (v2.2)
    A_NO_GB         = 0         (1)         : NO-driven presynaptic rate
    z_het_GB        = 0         (1)         : NO depresses terminals with tanh(Z) < z_het (Tong 2021)
    dpre_min_GB     = -0.8      (1)         : soft lower bound of dpre
    dpre_max_GB     = 1.0       (1)         : soft upper bound of dpre
    dpre0_GB        = 0         (1)         : initial presynaptic modifier
    pre_drive_GB    = 0         (1)         : T/N drive: 0 = shaft Ca (cai), 1 = spine VDCC influx (-ica_VDCC)
    i_scale_GB      = 1e-6      (nA)        : normalisation of the VDCC drive
    theta_Ti_GB     = 1e-7      (nA)        : VDCC threshold of the pre-LTD trace (pre_drive 1)
    theta_NOi_GB    = 1e-7      (nA)        : VDCC threshold of NO synthesis (pre_drive 1)
    no_drive_GB     = 0         (1)         : N drive: 0 = as pre_drive, 1 = effcai (v2.3), 2 = spine Ca cai_CR (v2.4)
    theta_NOe_GB    = 0.05      (us/liter)  : effcai threshold of NO synthesis (no_drive 1)
    e_scale_GB      = 0.01      (us/liter)  : normalisation of the effcai drive (no_drive 1)
    theta_NOc_GB    = 0         (mM)        : spine Ca (cai_CR - min_ca_CR) threshold of NO synthesis (no_drive 2, v2.4)
    c_scale_GB      = 1e-3      (mM)        : normalisation of the spine Ca drive (no_drive 2)
    theta_N_GB      = 0         (1)         : threshold on N for NO action, tanh(pos(N - theta_N)) (v2.4)
    : t_drive 3 (synapse-local eCB drive, D1; offline model_v2.py t_drive 3). Off at t_drive_GB = 0, theta_Tg_GB = 0.
    t_drive_GB      = 0         (1)         : T drive: 0 = pre_drive (Ca or VDCC), 3 = local-V event count
    theta_Te_GB     = 0.5       (1)         : threshold on the event trace S in the T drive
    Te_scale_GB     = 1         (1)         : normalisation of the T drive, u_T = pos(S - theta_Te)/Te_scale (1/ms)
    tau_E1_GB       = 100       (ms)        : event trace S decay
    theta_Tg_GB     = 0         (1)         : gate threshold at a pre spike, tanh(pos(T - theta_Tg)); 0 = old tanh(T)
    theta_V_GB      = 3         (mV)        : event = upward crossing of v - vbar above this
    w_V_GB          = 15        (ms)        : events within w_V after this synapse's own pre arrival are vetoed
    tau_b_GB        = 300       (ms)        : time constant of the slow local baseline vbar
    hyst_V_GB       = 2         (mV)        : re-arm when v - vbar falls below theta_V - hyst_V
    : V4 t_drive 4 (model_v2.DEFAULTS K_ca, CA_HYST, tau_d_NMDA)
    K_ca_GB         = 2.469e-8  (nA)        : -ica_VDCC event threshold
    ca_hyst_GB      = 0.5       (1)         : re-arm below K_ca_GB*ca_hyst_GB
    tau_glu_GB      = 70        (ms)        : decay of the own glutamate-bound state bglu
    : V4 C1/C2 gate (rho_v4.vamp). Off at vamp_mode_GB = 0 or theta_VA_GB <= 0
    vamp_mode_GB    = 0         (1)         : 0 off, 1 = C1, 2 = C2
    theta_VA_GB     = 0         (1)         : threshold on Vg (ms, -ica_VDCC/i_scale filtered with tau_E1)
    : V4 shadow expression (validation only)
    express_x_GB    = 0         (1)         : 1 = gmax_AMPA, Use follow rhox_GB (v1 rule, no dpre)
    gamma_dx_GB     = 101.5     (1)
    gamma_px_GB     = 199.773931 (1)
    theta_dx_GB     = 0.006     (us/liter)
    theta_px_GB     = 0.001     (us/liter)
    : V4 bin mode (pre path on the offline grid, model_v2.feature_T / feature_K / dpre_final). 0 = continuous v2 path
    bin_GB          = 0         (ms)        : grid step (extract_v2 decim 0.25); N driven by the bin MEAN of -ica_VDCC
    bin_t0_GB       = 0         (ms)        : grid start (prefire t = 0; full runs: the shifted induction start)
    bin_t1_GB       = 1e9       (ms)        : last grid event (full runs: the snap); no pre plasticity outside
    : V5 (gpu_v5_rho). Off at v5_mode_GB = 0 or A_eCB_GB = 0
    v5_mode_GB      = 0         (1)         : 0 off (V4), 1 = v5 (trigger Vg), 2 = v5b / v5c (trigger W)
    A_eCB_GB        = 0         (1)         : eCB step per gated own pre arrival (1 = to dpre_min)
    theta_eCB_GB    = 0         (1)         : threshold on the trigger pool (ms, as theta_VA_GB); <= 0: theta_VA_GB
    : V7 (gpu_v7_rho). veto_T_GB 0 and uE_GB 1 = V5
    veto_T_GB       = 0         (ms)        : bAP veto window after each own arrival (0 = no veto)
    uE_GB           = 1         (1)         : per-synapse scale of theta_eCB_GB (ecb_ref 2: 1e5 vdcc_q_post)
    : V8 (gpu_v8_rho G1). gate_src_GB 0 = V7
    gate_src_GB     = 0         (1)         : LTP licence: 0 = Vg (V7), 2 = shaft dCa windowed max
    gate_win_GB     = 100       (ms)        : licence window after the last shaft dCa above theta_sh_GB
    theta_sh_GB     = 1.5e-4    (mM)        : shaft dCa threshold (0.15 uM)
    t_rest_GB       = -1        (ms)        : cai_rest_GB follows cai while t <= t_rest_GB; < 0: user-set cai_rest_GB
}


VERBATIM
/**
 * This Verbatim block is needed to generate random numbers from a uniform
 * distribution U(0, 1).
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef NRN_VERSION_GTEQ_8_2_0
#include "nrnran123.h"
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

#ifndef CORENEURON_BUILD
extern int ifarg(int iarg);

extern void* vector_arg(int iarg);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
#endif
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif


ENDVERBATIM


ASSIGNED {
    g_AMPA          (uS)    : AMPA Receptor
    g_NMDA          (uS)    : NMDA Receptor
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    rng                     : Random Number Generator
    usingR123               : TEMPORARY until mcellran4 completely deprecated
    ica_NMDA        (nA)    : NMDAR-mediated calcium current
    ica_VDCC        (nA)    : VDCC (R-type)
    : Long-term synaptic plasticity
    dep_GB          (1)
    Use_tgt_GB      (1)     : v2 target of Use_GB
    pot_GB          (1)
    : Shaft calcium (ion variable, auto-populated by NEURON from CaDynamics.cad)
    cai             (mM)
    shaft_cai       (mM)
    : Misc
    v               (mV)
    vsyn            (mV)
    i               (nA)

    : stuff for delayed connections
    delay_times
    delay_weights
    next_delay (ms)
    i_AMPA (nA)
    i_NMDA (nA)
    tlast_GB (ms)   : t_drive 3: time of this synapse's last pre arrival
    armed_GB (1)    : t_drive 3: 1 = the next upward crossing is an event
    armca_GB (1)    : t_drive 4: 1 = the next -ica_VDCC crossing of K_ca is an event
    nev_GB   (1)    : t_drive 4: number of events (diagnostic)
    gv_GB    (1)    : C1/C2 gate, 1 = Vg_GB > theta_VA_GB
    potx_GB  (1)    : shadow rule
    depx_GB  (1)    : shadow rule
    : V4 bin mode
    Nd_GB    (1)    : N on the grid (N[k], left-point lowpass of the bin-mean drive)
    Td_GB    (1)    : T on the grid (T[k])
    gprev_GB (1)    : tanh(pos(N[k-1] - theta_N)) tanh(pos(Z[k-1] - theta_Z))
    uTprev_GB (1)   : pos(S[k-1] - theta_Te)/Te_scale
    npend_GB (1)    : own arrivals since the last grid event (mGluR steps applied at the next grid sample)
    binon_GB (1)    : 0 before bin_t0, 1 running, 2 after bin_t1
    nbin_GB  (1)    : number of mGluR steps applied (diagnostic)
    tTsum_GB (1)    : sum of tanh(pos(T - theta_Tg)) over the applied steps (diagnostic)
    : V5
    necb_GB  (1)    : number of eCB steps applied (diagnostic)
    : V7
    nveto_GB (1)    : number of triggered arrivals vetoed (diagnostic)
    nvpend_GB (1)   : veto decisions pending (self-events not yet delivered)
    : V8
    cai_rest_GB (mM) : pre-stimulus shaft cai of this synapse
    gsh_GB   (1)    : shaft licence, 1 = open
    ash_GB   (1)    : 1 = dsh currently above theta_sh_GB
    tlsh_GB  (ms)   : time of the last fall of dsh below theta_sh_GB
    ton_GB   (ms)   : time the licence last opened
    nsh_GB   (1)    : number of licence openings (diagnostic)
    tsh_GB   (ms)   : licence open time of the closed intervals (diagnostic)
}

STATE {
    : AMPA Receptor
    A_AMPA      (1)
    B_AMPA      (1)
    gmax_AMPA   (nS)
    : NMDA Receptor
    A_NMDA      (1)
    B_NMDA      (1)
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use_GB      (1)
    : VDCC (R-type)
    m_VDCC      (1)
    h_VDCC      (1)
    : Postsynaptic Ca2+ dynamics
    cai_CR      (mM)        <1e-6>
    : Long-term synaptic plasticity
    rho_GB      (1)
    effcai_GB   (us/liter)  <1e-3>
    : v2 presynaptic plasticity
    dpre_GB     (1)
    T_GB        (1)
    N_GB        (1)
    Z_GB        (1)
    : t_drive 3
    S_GB        (1)
    vbar_GB     (mV)
    : V4
    bglu_GB     (1)
    Vg_GB       (1)
    rhox_GB     (1)
    qca_GB      (1)     <1e-4>  : V4 bin mode: charge of -ica_VDCC/i_scale since the last grid event (ms)
    : V5
    W_GB        (1)             : (1 - min(bglu, 1))-weighted copy of Vg, the eCB trigger pool (v5_mode 2)
}

INITIAL{
    : AMPA Receptor
    A_AMPA      = 0
    B_AMPA      = 0
    gmax_AMPA   = gmax0_AMPA
    : NMDA Receptor
    A_NMDA      = 0
    B_NMDA      = 0
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    if (uncaging_mode) {
        Use_GB  = 1
    } else {
        Use_GB  = Use
    }
    : Postsynaptic Ca2+ dynamics
    cai_CR      = min_ca_CR
    : Long-term synaptic plasticity
    rho_GB      = rho0_GB
    effcai_GB   = 0
    dpre_GB     = dpre0_GB
    T_GB        = 0
    N_GB        = 0
    Z_GB        = 0
    S_GB        = 0
    vbar_GB     = v
    tlast_GB    = -1e9
    armed_GB    = 1
    bglu_GB     = 0
    Vg_GB       = 0
    rhox_GB     = rho0_GB
    armca_GB    = 1
    nev_GB      = 0
    gv_GB       = 0
    potx_GB     = 0
    depx_GB     = 0
    qca_GB      = 0
    Nd_GB       = 0
    Td_GB       = 0
    gprev_GB    = 0
    uTprev_GB   = 0
    npend_GB    = 0
    binon_GB    = 0
    nbin_GB     = 0
    tTsum_GB    = 0
    W_GB        = 0
    necb_GB     = 0
    nveto_GB    = 0
    nvpend_GB   = 0
    gsh_GB      = 0
    ash_GB      = 0
    tlsh_GB     = -1e9
    ton_GB      = 0
    nsh_GB      = 0
    tsh_GB      = 0
    if (t_rest_GB >= 0) { cai_rest_GB = cai }
    Use_tgt_GB  = Use
    dep_GB      = 0
    pot_GB      = 0
    : Delayed connection
    next_delay = -1

    : Initialize watchers
    net_send(0, 1)
    : V4 bin mode grid (flag 21)
    if (bin_GB > 0) { net_send(bin_t0_GB, 21) }

}

PROCEDURE setup_delay_vecs() {
VERBATIM
#ifndef CORENEURON_BUILD
    void** vv_delay_times = (void**)(&_p_delay_times);
    void** vv_delay_weights = (void**)(&_p_delay_weights);
    *vv_delay_times = (void*)NULL;
    *vv_delay_weights = (void*)NULL;
    if (ifarg(1)) {
        *vv_delay_times = vector_arg(1);
    }
    if (ifarg(2)) {
        *vv_delay_weights = vector_arg(2);
    }
#endif
ENDVERBATIM
}


BREAKPOINT {
    LOCAL Eca_syn, mggate, Pf_NMDA, gca_bar_abs_VDCC, gca_VDCC
    SOLVE state METHOD euler
    : AMPA Receptor
    g_AMPA = (1e-3)*gmax_AMPA*(B_AMPA - A_AMPA)
    i_AMPA = g_AMPA*(v-E_AMPA)
    : NMDA Receptor
    mggate = 1 / (1 + exp(-slope_NMDA*v) * (mg/scale_NMDA))
    g_NMDA = (1e-3)*gmax_NMDA*mggate*(B_NMDA - A_NMDA)
    i_NMDA = g_NMDA*(v - E_NMDA)
    : NMDAR-mediated calcium current
    Pf_NMDA  = (4*cao_CR) / (4*cao_CR + (1/1.38) * 120 (mM)) * 0.6
    ica_NMDA = Pf_NMDA*g_NMDA*(v-40.0)
    : VDCC (R-type), assuming sphere for spine head
    gca_bar_abs_VDCC = gca_bar_VDCC * 4(um2)*PI*(3(1/um3)/4*volume_CR*1/PI)^(2/3)
    gca_VDCC = (1e-3) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC
    Eca_syn = nernst(cai_CR, cao_CR, 2)
    ica_VDCC = gca_VDCC*(v-Eca_syn)*calcium_current_flag
    : Update synaptic voltage (for recording convenience)
    vsyn = v
    : Update current
    i = i_AMPA + i_NMDA + ica_VDCC
    : Copy shaft calcium for recording
    shaft_cai = cai
    : V8 pre-stimulus shaft rest (frozen after t_rest_GB)
    if (gate_src_GB > 1.5 && t_rest_GB >= 0 && t <= t_rest_GB) { cai_rest_GB = cai }
}


DERIVATIVE state {
    LOCAL minf_VDCC, hinf_VDCC, rhoE, potE, depE, cont
    : V4 bin mode: dpre changes only at the grid events (flag 21) and arrivals
    cont = 1
    if (bin_GB > 0) { cont = 0 }
    : V4: expression state (rho_GB, or the shadow rhox_GB when express_x_GB = 1)
    if (express_x_GB > 0.5) { rhoE = rhox_GB } else { rhoE = rho_GB }
    : AMPA Receptor
    A_AMPA'      = - A_AMPA/tau_r_AMPA
    B_AMPA'      = - B_AMPA/tau_d_AMPA
    gmax_AMPA'   = (gmax_d_AMPA + rhoE*(gmax_p_AMPA - gmax_d_AMPA) - gmax_AMPA) / ((1e3)*tau_exp_GB)
    : NMDA Receptor
    A_NMDA'      = - A_NMDA/tau_r_NMDA
    B_NMDA'      = - B_NMDA/tau_d_NMDA
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    if (uncaging_mode) {
        Use_GB'  = 0
    } else {
        if (express_x_GB > 0.5) {
            Use_tgt_GB = Use_d + rhox_GB*(Use_p - Use_d)
        } else {
            Use_tgt_GB = (Use_d + rho_GB*(Use_p - Use_d)) * (1 + dpre_GB)
        }
        if (Use_tgt_GB > 1) { Use_tgt_GB = 1 }
        Use_GB'  = (Use_tgt_GB - Use_GB) / ((1e3)*tau_exp_GB)
    }
    : V4 C1 / C2 gate (scan_vgate_amp.rho_rec): pot needs Vg > theta_VA; C2 also drops dep above theta_p
    : V8 gate_src 2: the shaft licence replaces the Vg licence (gpu_v8_rho kernel: pot = 0 when thG > 0 and not lic)
    potE = pot_GB
    depE = dep_GB
    if (gate_src_GB > 1.5) {
        if (theta_sh_GB > 0) { potE = pot_GB*gsh_GB }
    } else {
        if (vamp_mode_GB > 0.5 && theta_VA_GB > 0) { potE = pot_GB*gv_GB }
    }
    if (vamp_mode_GB > 1.5) { depE = dep_GB*(1 - pot_GB) }
    : VDCC (R-type)
    minf_VDCC    = 1 / (1 + exp(((vhm_VDCC - ljp_VDCC) - v) / km_VDCC))
    hinf_VDCC    = 1 / (1 + exp(((vhh_VDCC - ljp_VDCC) - v) / kh_VDCC))
    m_VDCC'      = (minf_VDCC-m_VDCC)/mtau_VDCC
    h_VDCC'      = (hinf_VDCC-h_VDCC)/htau_VDCC
    : Long-term synaptic plasticity
    effcai_GB'   = - effcai_GB/tau_effca_GB + (cai_CR - min_ca_CR)
    
    rho_GB'      = ( - rho_GB*(1 - rho_GB)*(rho_star_GB - rho_GB)
                 + potE*gamma_p_GB*(1 - rho_GB)
                 - depE*(1 - potE)*gamma_d_GB*rho_GB ) / ((1e3)*tau_ind_GB)
    : V4 shadow v1 rule (express_x_GB); runs always, only read when express_x_GB = 1
    rhox_GB'     = ( - rhox_GB*(1 - rhox_GB)*(rho_star_GB - rhox_GB)
                 + potx_GB*gamma_px_GB*(1 - rhox_GB)
                 - depx_GB*(1 - potx_GB)*gamma_dx_GB*rhox_GB ) / ((1e3)*tau_ind_GB)
    : V4 own glutamate-bound state and filtered VDCC charge
    bglu_GB'     = - bglu_GB/tau_glu_GB
    Vg_GB'       = - Vg_GB/tau_E1_GB + vdcc_drive()
    qca_GB'      = vdcc_drive()
    : V5 eCB trigger pool (v5_mode 2): same filter as Vg, input weighted by 1 - min(bglu, 1)
    W_GB'        = - W_GB/tau_E1_GB + ecb_drive()
    : Postsynaptic Ca2+ dynamics
    : v2 presynaptic plasticity. Drive: shaft Ca (cai, CaDynamics) or spine VDCC influx, see pre_drive
    T_GB'        = - T_GB/tau_T_GB + T_drive()
    N_GB'        = - N_GB/tau_NO_GB + no_drive()
    Z_GB'        = - Z_GB/tau_Z_GB
    S_GB'        = - S_GB/tau_E1_GB
    vbar_GB'     = (v - vbar_GB)/tau_b_GB
    dpre_GB'     = cont*A_NO_GB*tanh(pos(N_GB - theta_N_GB))*( pos(tanh(pos(Z_GB - theta_Z_GB)) - z_het_GB)*(dpre_max_GB - dpre_GB)
                                     - pos(z_het_GB - tanh(pos(Z_GB - theta_Z_GB)))*(dpre_GB - dpre_min_GB) )
                   / ((1e3)*tau_ind_GB)
    cai_CR'      = - (1e-9)*(ica_NMDA+ica_VDCC)*gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY)
                   - (cai_CR - min_ca_CR)/tau_ca_CR
}


NET_RECEIVE (weight, u, tsyn (ms), recovered, unrecovered, nc_type) {
    : nc_type: 0=presynaptic netcon, 1=spontmini, 2=replay
    LOCAL p_rec, released, tp, factor, rec, thE, trig
    INITIAL {
        weight = 1
        if (uncaging_mode) {
            u = 1
            recovered = Nrrp
            unrecovered = 0
        } else {
            u = 0
            if (init_depleted){
                recovered = 0
                unrecovered = Nrrp
            } else {
                recovered = Nrrp
                unrecovered = 0
            }
        }
        tsyn = 0 (ms)
        if (nc_type == 0) {   : pre-synaptic netcon
    VERBATIM
            // setup self events for delayed connections to change weights
            IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
            IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
            if (vv_delay_times && vector_capacity(vv_delay_times)>=1) {
                double* deltm_el = vector_vec(vv_delay_times);
                int delay_times_idx;
                next_delay = 0;
                for(delay_times_idx = 0; delay_times_idx < vector_capacity(vv_delay_times); ++delay_times_idx) {
                    double next_delay_t = deltm_el[delay_times_idx];
    ENDVERBATIM
                    net_send(next_delay_t, 10)  : use flag 10 to avoid interfering with GluSynapse logic
    VERBATIM
                }
            }
    ENDVERBATIM
        }
    }


    if(verbose > 0){ UNITSOFF printf("Time = %g ms, incoming spike at synapse %g\n", t, synapseID) UNITSON }
    if(flag == 0) {
        if(weight <= 0){
            : Do not perform any calculations if the synapse (netcon) is deactivated.
            : This avoids drawing from the random stream
            : WARNING In this model *weight* is only used to activate/deactivate the
            :         synapse. The conductance is stored in gmax_AMPA and gmax_NMDA.
            if(verbose > 0){ printf("Inactive synapse, weight = %g\n", weight) }
        } 
        else {
                if (uncaging_mode) {
                    : Uncaging mode
                    u = 1
                    recovered = Nrrp
                    released = Nrrp
                } 
                else {
                    : Flag 0: Regular spike
                    if(verbose > 0){ printf("Flag 0, Regular spike\n") }
                    : Update facilitation variable as Eq. 2 in Fuhrmann et al. 2002
                    u = Use_GB + u*(1 - Use_GB)*exp(-(t - tsyn)/Fac)
                    if ( verbose > 0 ) { printf("\tVesicle release probability = %g\n", u) }
                    : Recovery
                    p_rec = 1 - exp(-(t - tsyn)/Dep)
                    if ( verbose > 0 ) { printf("\tVesicle recovery probability = %g\n", p_rec) }
                    if ( verbose > 0 ) { printf("\tVesicle available before recovery = %g\n", recovered) }
                    recovered = recovered + brand(unrecovered, p_rec)
                    if ( verbose > 0 ) { printf("\tVesicles available after recovery = %g\n", recovered) }
                    : Release
                    rec = recovered  : Make a copy so we can change it for single vesicle minis w/o messing with recovered
                    : Consider only a single recovered vesicle for minis (if minis_single_vesicle flag is set to 1)
                    if (rec > 1 && minis_single_vesicle && nc_type == 1) { rec = 1 }
                    released = brand(rec, u)
                    if ( verbose > 0 ) { printf("\tReleased %g vesicles out of %g\n", released, recovered) }
                    : Update vesicle pool
                    recovered = recovered - released
                    unrecovered = Nrrp - recovered
                    if ( verbose > 0 ) { printf("\tFinal vesicle count, Recovered = %g, Unrecovered = %g, Nrrp = %g\n", recovered, unrecovered, Nrrp) }
                }
                
                : Update AMPA variables
                tp = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)  : Time to peak
                factor = 1 / (-exp(-tp/tau_r_AMPA)+exp(-tp/tau_d_AMPA))  : Normalization factor
                A_AMPA = A_AMPA + released/Nrrp*factor
                B_AMPA = B_AMPA + released/Nrrp*factor
                : Update NMDA variables
                tp = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)  : Time to peak
                factor = 1 / (-exp(-tp/tau_r_NMDA)+exp(-tp/tau_d_NMDA))  : Normalization factor
                A_NMDA = A_NMDA + released/Nrrp*factor
                B_NMDA = B_NMDA + released/Nrrp*factor
                : v2: mGluR/eCB pre-LTD and presynaptic trace, once per presynaptic spike
                if (!uncaging_mode) {
                    if (bin_GB > 0) {
                        : V4 bin mode: the mGluR step waits for the next grid sample (arrival_index, first sample
                        : at or after the arrival), after the NO increment of the preceding step (dpre_final)
                        if (binon_GB > 0.5 && binon_GB < 1.5) { npend_GB = npend_GB + 1 }
                    } else {
                        dpre_GB = dpre_GB - A_mglu_GB*tanh(pos(T_GB - theta_Tg_GB))*(dpre_GB - dpre_min_GB)
                    }
                    Z_GB = Z_GB + 1
                    tlast_GB = t
                    : V5 eCB-LTD at this own arrival, read before its bglu increment (gpu_v5_rho step order)
                    if (v5_mode_GB > 0.5 && A_eCB_GB > 0) {
                        : V7: theta_eCB,i = theta_eCB_GB uE_GB (gpu_v7_rho thE0 uE[i]; <= 0: theta_VA_GB)
                        if (theta_eCB_GB > 0) { thE = theta_eCB_GB*uE_GB } else { thE = theta_VA_GB }
                        if (v5_mode_GB > 1.5) { trig = W_GB } else { trig = Vg_GB }
                        if (trig > thE) {
                            if (veto_T_GB > 0) {
                                : V7 bAP veto: decided at t + veto_T on the VDCC charge since this arrival
                                nvpend_GB = nvpend_GB + 1
                                net_send(veto_T_GB, 1000 + qca_GB)
                            } else {
                                dpre_GB = dpre_min_GB + (dpre_GB - dpre_min_GB)*(1 - A_eCB_GB)
                                necb_GB = necb_GB + 1
                            }
                        }
                    }
                    : V4 t_drive 4: own arrival (released or not, as model_v2.glu_weight)
                    bglu_GB = bglu_GB + 1
                }
                : Update tsyn
                : tsyn knows about all spikes, not only those that released
                : i.e. each spike can increase the u, regardless of recovered state
                :      and each spike trigger an evaluation of recovery
                tsyn = t
        }
    } else if(flag == 1) {
        : Flag 1, Initialize watchers
        if(verbose > 0){ printf("Flag 1, Initialize watchers\n") }
        WATCH (effcai_GB > theta_d_GB) 2
        WATCH (effcai_GB < theta_d_GB) 3
        WATCH (effcai_GB > theta_p_GB) 4
        WATCH (effcai_GB < theta_p_GB) 5
        : t_drive 3 event detector, see the comment above FUNCTION T_drive
        if (t_drive_GB > 2.5 && t_drive_GB < 3.5) {
            WATCH (v - vbar_GB > theta_V_GB) 6
            WATCH (v - vbar_GB < theta_V_GB - hyst_V_GB) 7
        }
        : V4 t_drive 4 event detector on -ica_VDCC (flags 8, 9)
        if (t_drive_GB > 3.5) {
            WATCH (ica_VDCC < - K_ca_GB) 8
            WATCH (ica_VDCC > - K_ca_GB*ca_hyst_GB) 9
        }
        : V4 C1 / C2 gate (flags 11, 12)
        if (vamp_mode_GB > 0.5) {
            WATCH (Vg_GB > theta_VA_GB) 11
            WATCH (Vg_GB < theta_VA_GB) 12
        }
        : V4 shadow v1 rule (flags 13-16)
        if (express_x_GB > 0.5) {
            WATCH (effcai_GB > theta_dx_GB) 13
            WATCH (effcai_GB < theta_dx_GB) 14
            WATCH (effcai_GB > theta_px_GB) 15
            WATCH (effcai_GB < theta_px_GB) 16
        }
        : V8 shaft gate (flags 22, 23; closing self-event 24)
        if (gate_src_GB > 1.5) {
            WATCH (cai - cai_rest_GB > theta_sh_GB) 22
            WATCH (cai - cai_rest_GB < theta_sh_GB) 23
        }
    } else if(flag == 2) {
        : Flag 2, Activate depression mechanisms
        if(verbose > 0){ printf("Flag 2, Activate depression mechanisms\n") }
        dep_GB = 1
    } else if(flag == 3) {
        : Flag 3, Deactivate depression mechanisms
        if(verbose > 0){ printf("Flag 3, Deactivate depression mechanisms\n") }
        dep_GB = 0
    } else if(flag == 4) {
        : Flag 4, Activate potentiation mechanisms
        if(verbose > 0){ printf("Flag 4, Activate potentiation mechanisms\n") }
        pot_GB = 1
    } else if(flag == 5) {
        : Flag 5, Deactivate potentiation mechanisms
        if(verbose > 0){ printf("Flag 5, Deactivate potentiation mechanisms\n") }
        pot_GB = 0
    } else if(flag == 6) {
        : Flag 6, t_drive 3: upward crossing of v - vbar over theta_V. One event per excursion (armed),
        : dropped if it falls within w_V after this synapse's own pre arrival (the veto still uses the arm)
        if (armed_GB > 0.5) {
            armed_GB = 0
            if (t - tlast_GB >= w_V_GB) { S_GB = S_GB + 1 }
        }
    } else if(flag == 7) {
        : Flag 7, t_drive 3: v - vbar fell below theta_V - hyst_V, re-arm
        armed_GB = 1
    } else if(flag == 8) {
        : Flag 8, t_drive 4: upward crossing of -ica_VDCC through K_ca, one event per excursion, weight 1 - min(bglu, 1)
        if (armca_GB > 0.5) {
            armca_GB = 0
            nev_GB = nev_GB + 1
            if (bglu_GB < 1) { S_GB = S_GB + 1 - bglu_GB }
        }
    } else if(flag == 9) {
        : Flag 9, t_drive 4: -ica_VDCC fell below K_ca*ca_hyst, re-arm
        armca_GB = 1
    } else if(flag == 11) {
        gv_GB = 1
    } else if(flag == 12) {
        gv_GB = 0
    } else if(flag == 13) {
        depx_GB = 1
    } else if(flag == 14) {
        depx_GB = 0
    } else if(flag == 15) {
        potx_GB = 1
    } else if(flag == 16) {
        potx_GB = 0
    } else if(flag == 21) {
        : V4 bin mode grid event at t_k = bin_t0 + k bin (see PROCEDURE bin_step)
        bin_step()
        if (binon_GB < 1.5) { net_send(bin_GB, 21) }
    } else if(flag == 22) {
        : V8: shaft dCa rose above theta_sh, the licence opens (or stays open)
        ash_GB = 1
        if (gsh_GB < 0.5) {
            gsh_GB = 1
            nsh_GB = nsh_GB + 1
            ton_GB = t
        }
    } else if(flag == 23) {
        : V8: shaft dCa fell below theta_sh, the licence stays open for gate_win
        ash_GB = 0
        tlsh_GB = t
        net_send(gate_win_GB, 24)
    } else if(flag == 24) {
        : V8: close unless dsh rose again or fell again later (a stale self-event)
        if (ash_GB < 0.5 && gsh_GB > 0.5 && t - tlsh_GB >= gate_win_GB - 1e-9) {
            gsh_GB = 0
            tsh_GB = tsh_GB + t - ton_GB
        }
    } else if(flag > 500) {
        : V7 veto decision for a triggered own arrival at t - veto_T_GB; flag = 1000 + qca_GB at that arrival
        nvpend_GB = nvpend_GB - 1
        if (theta_eCB_GB > 0) { thE = theta_eCB_GB*uE_GB } else { thE = theta_VA_GB }
        if (qca_GB - (flag - 1000) > thE) {
            nveto_GB = nveto_GB + 1
        } else {
            dpre_GB = dpre_min_GB + (dpre_GB - dpre_min_GB)*(1 - A_eCB_GB)
            necb_GB = necb_GB + 1
        }
    } else if(flag == 10) {
        : Flag 10, Handle delayed connection weight changes
    VERBATIM
        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {
            double* weights_v = vector_vec(vv_delay_weights);
            double next_delay_weight = weights_v[(int)next_delay];
    ENDVERBATIM
            weight = conductance * next_delay_weight
            next_delay = next_delay + 1
    VERBATIM
        }
    ENDVERBATIM
    }
}

FUNCTION nernst(ci(mM), co(mM), z) (mV) {
    nernst = (1000) * R * (celsius + 273.15) / (z*FARADAY) * log(co/ci)
    if(verbose > 1) { UNITSOFF printf("nernst:%g R:%g temperature (c):%g \n", nernst, R, celsius) UNITSON }
}

PROCEDURE setRNG() {
    VERBATIM
    #ifndef CORENEURON_BUILD
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
    #endif
    ENDVERBATIM
}


PROCEDURE clearRNG() {
VERBATIM
    #ifndef CORENEURON_BUILD
    if (usingR123) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
    } else {
        void** pv = (void**)(&_p_rng);
        if (*pv) {
            *pv = (void*)0;
        }
    }
    #endif
ENDVERBATIM
}


FUNCTION urand() {
    VERBATIM
    double value;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(RANDCAST _p_rng);
        #endif
    } else {
        value = 0.0;
    }
    _lurand = value;
    ENDVERBATIM
}

FUNCTION brand(n, p) {
    LOCAL result, count, success
    success = 0
    FROM count = 0 TO (n - 1) {
        result = urand()
        if(result <= p) {
            success = success + 1
        }
    }
    brand = success
}


FUNCTION bbsavestate() {
    bbsavestate = 0
    VERBATIM
    #ifndef CORENEURON_BUILD
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval;
        #ifndef NRN_VERSION_GTEQ_8_2_0
        double *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        #endif
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(RANDCAST _p_rng);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(RANDCAST _p_rng, (long)(xval[0]));
                }
            }
        }
    #endif
    ENDVERBATIM
}


VERBATIM
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
    IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

    // make sure offset array non-null
    if (iArray) {
        // get handle to random123 instance
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        // get location for storing ids
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        // retrieve/store identifier seeds
        nrnran123_getids3(*pv, ia, ia+1, ia+2);
        // retrieve/store stream sequence
        char which;
        nrnran123_getseq(*pv, ia+3, &which);
        ia[4] = (int)which;
    }

    // increment integer offset (2 identifier), no double data
    *ioffset += 5;
    *doffset += 0;

    // serialize connection delay vectors
    if (vv_delay_times && vv_delay_weights &&
       (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {
        if (iArray) {
            uint32_t* di = ((uint32_t*)iArray) + *ioffset;
            // store vector sizes for deserialization
            di[0] = vector_capacity(vv_delay_times);
            di[1] = vector_capacity(vv_delay_weights);
        }
        if (dArray) {
            double* delay_times_el = vector_vec(vv_delay_times);
            double* delay_weights_el = vector_vec(vv_delay_weights);
            double* x_i = dArray + *doffset;
            int delay_vecs_idx;
            int x_idx = 0;
            for(delay_vecs_idx = 0; delay_vecs_idx < vector_capacity(vv_delay_times); ++delay_vecs_idx) {
                 x_i[x_idx++] = delay_times_el[delay_vecs_idx];
                 x_i[x_idx++] = delay_weights_el[delay_vecs_idx];
            }
        }
        // reserve space for connection delay data on serialization buffer
        *doffset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);
    } else {
        if (iArray) {
            uint32_t* di = ((uint32_t*)iArray) + *ioffset;
            di[0] = 0;
            di[1] = 0;
        }
    }
    // reserve space for delay vectors (may be 0)
    *ioffset += 2;
}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure it's not previously set
    assert(!_p_rng);
    assert(!_p_delay_times && !_p_delay_weights);

    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
    // make sure non-zero identifier seeds
    if (ia[0] != 0 || ia[1] != 0 || ia[2] != 0) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        // get new stream
        *pv = nrnran123_newstream3(ia[0], ia[1], ia[2]);
        // restore sequence
        nrnran123_setseq(*pv, ia[3], (char)ia[4]);
    }
    // increment intger offset (2 identifiers), no double data
    *ioffset += 5;

    int delay_times_sz = iArray[5];
    int delay_weights_sz = iArray[6];
    *ioffset += 2;

    if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {
        double* x_i = dArray + *doffset;

        // allocate vectors
        _p_delay_times = (double*)vector_new1(delay_times_sz);
        _p_delay_weights = (double*)vector_new1(delay_weights_sz);

        double* delay_times_el = vector_vec((IvocVect*)_p_delay_times);
        double* delay_weights_el = vector_vec((IvocVect*)_p_delay_weights);

        // copy data
        int x_idx;
        int vec_idx = 0;
        for(x_idx = 0; x_idx < delay_times_sz + delay_weights_sz; x_idx += 2) {
            delay_times_el[vec_idx] = x_i[x_idx];
            delay_weights_el[vec_idx++] = x_i[x_idx+1];
        }
        *doffset += delay_times_sz + delay_weights_sz;
    }
}
ENDVERBATIM


FUNCTION T_drive() {
    : drive of the T (eCB) trace. t_drive 3 (synapse-local, D1): u_T = pos(S - theta_Te)/Te_scale, S counting
    : local-V events. Detection: WATCH on v - vbar in NET_RECEIVE (flags 6, 7). WATCH is evaluated every fixed
    : step and, under CVODE, root-found, so the event time is exact there and the first step after the crossing
    : at fixed dt (the offline detector's grid sample). A BREAKPOINT check with net_send would be evaluated
    : only at accepted steps under CVODE, and delayed by one event queue pass at fixed step.
    : Re-arm hysteresis hyst_V. vbar_GB' = (v - vbar_GB)/tau_b_GB is a STATE (euler / CVODE alike).
    UNITSOFF
    if (t_drive_GB > 2.5) {
        T_drive = pos(S_GB - theta_Te_GB)/Te_scale_GB
    } else {
        T_drive = pre_drive(theta_T_GB, theta_Ti_GB)
    }
    UNITSON
}

FUNCTION pre_drive(thc (mM), thi (nA)) {
    : drive of the T and N traces (1/ms after dividing by the scale; UNITSOFF for the mixed branches)
    UNITSOFF
    if (pre_drive_GB > 0.5) {
        pre_drive = pos(-ica_VDCC - thi)/i_scale_GB
    } else {
        pre_drive = pos(cai - ca_sh_rest_GB - thc)/ca_scale_GB
    }
    UNITSON
}

FUNCTION no_drive() {
    : drive of the N (NO) trace: effcai above theta_NOe (no_drive 1), else the T-pathway form
    UNITSOFF
    if (no_drive_GB > 1.5) {
        no_drive = pos(cai_CR - min_ca_CR - theta_NOc_GB)/c_scale_GB
    } else if (no_drive_GB > 0.5) {
        no_drive = pos(effcai_GB - theta_NOe_GB)/e_scale_GB
    } else {
        no_drive = pre_drive(theta_NO_GB, theta_NOi_GB)
    }
    UNITSON
}

FUNCTION vdcc_drive() {
    : V4 C1 / C2: drive of Vg, -ica_VDCC / i_scale (1/ms with tau_E1 in ms, as scan_vgate_amp bV)
    UNITSOFF
    vdcc_drive = -ica_VDCC/i_scale_GB
    UNITSON
}

FUNCTION ecb_drive() {
    : V5 drive of W (gpu_v5_rho TRIG_W): (1 - min(bglu, 1)) (-ica_VDCC)/i_scale in v5_mode 2, else 0
    UNITSOFF
    if (v5_mode_GB > 1.5 && bglu_GB < 1) {
        ecb_drive = (1 - bglu_GB)*(-ica_VDCC)/i_scale_GB
    } else {
        ecb_drive = 0
    }
    UNITSON
}

FUNCTION pos(x) {
    : rectifier [x]+
    if (x > 0) { pos = x } else { pos = 0 }
}

PROCEDURE bin_step() {
    : V4 bin mode, the offline pre path on a uniform grid of step bin_GB (model_v2, extract_v2 --decim 0.25):
    :   vdcc[k-1] = mean of -ica_VDCC over [t_{k-1}, t_k) (extract_v2.step_mean, charge-preserving), so the N drive
    :   is pos(mean - theta_NOi)/i_scale: the bin mean is rectified, not the instantaneous current.
    :   N[k] = aN N[k-1] + bN uN[k-1], T[k] = aT T[k-1] + bT uT[k-1] (model_v2._lowpass_grid / _coef);
    :   S, Z jump at the first sample at or after the event / arrival (_impulse_grid: the continuous states here).
    :   K: left-point, d <- dmax - (dmax - d) exp(-r g[k-1] bin) (dpre_final composes exactly over the steps);
    :   then one mGluR step per arrival in (t_{k-1}, t_k] with tanh(pos(T[k] - theta_Tg)) (feature_T at arr).
    LOCAL aN, bN, aT, bT, uN, n
    UNITSOFF
    if (binon_GB < 0.5) {
        : k = 0 (x[0] = 0)
        binon_GB = 1
        Nd_GB = 0
        Td_GB = 0
        gprev_GB = 0
        npend_GB = 0
    } else {
        if (bin_GB/tau_NO_GB < 0.1) {
            aN = 1 - bin_GB/tau_NO_GB
            bN = bin_GB
        } else {
            aN = exp(-bin_GB/tau_NO_GB)
            bN = tau_NO_GB*(1 - aN)
        }
        if (bin_GB/tau_T_GB < 0.1) {
            aT = 1 - bin_GB/tau_T_GB
            bT = bin_GB
        } else {
            aT = exp(-bin_GB/tau_T_GB)
            bT = tau_T_GB*(1 - aT)
        }
        uN = pos(qca_GB/bin_GB - theta_NOi_GB/i_scale_GB)
        dpre_GB = dpre_max_GB - (dpre_max_GB - dpre_GB)*exp(-A_NO_GB*gprev_GB*bin_GB/((1e3)*tau_ind_GB))
        Nd_GB = aN*Nd_GB + bN*uN
        Td_GB = aT*Td_GB + bT*uTprev_GB
        n = npend_GB
        while (n > 0.5) {
            dpre_GB = dpre_GB - A_mglu_GB*tanh(pos(Td_GB - theta_Tg_GB))*(dpre_GB - dpre_min_GB)
            nbin_GB = nbin_GB + 1
            tTsum_GB = tTsum_GB + tanh(pos(Td_GB - theta_Tg_GB))
            n = n - 1
        }
        npend_GB = 0
        gprev_GB = tanh(pos(Nd_GB - theta_N_GB))*tanh(pos(Z_GB - theta_Z_GB))
    }
    uTprev_GB = pos(S_GB - theta_Te_GB)/Te_scale_GB
    qca_GB = 0
    if (t + bin_GB > bin_t1_GB + 1e-6) { binon_GB = 2 }
    UNITSON
}
