#!/usr/bin/env python3
"""
CICR Model — Simplified Two-Stage IP3R Trigger + RyR Amplifier

Replaces the full Li-Rinzel IP3R model with a minimal biophysically-motivated
formulation. No gating ODEs, no stiffness — all cnexp-solvable in NEURON.

Two-stage ER release cascade:
  Stage 1 (IP3R): J_IP3R = V_IP3R * ip3_act * gradient       (trigger)
  Stage 2 (RyR):  J_RyR  = V_RyR * Hill(ca_cicr, K_RyR, 3) * gradient  (amplifier)

where ip3_act = IP3 / (IP3 + K_IP3)  (saturating activation, like Li-Rinzel m_inf)

RyR is gated by ca_cicr (ER-released calcium only), NOT total ca_i.
This ensures bAPs alone cannot trigger CICR — IP3R must release Ca first.

Key biology (Hong & Ross 2007):
  - IP3 is REQUIRED for release (bAPs alone don't trigger it)
  - SERCA loads ER during bAPs (priming, lasts 1-2s)
  - ER loading state (gradient) is the slow memory variable
  - IP3 provides the fast coincidence signal from presynaptic input

Priming timescale: tau_priming = gamma_ER / V_leak ~ 1-2 seconds
  V_leak bounds [5e-5, 2e-4] /ms give tau_priming in [0.6s, 2.5s]

Refractory inactivation gate (h_ref):
  h_inf = K_h^2 / (K_h^2 + ca_cicr^2)
  dh/dt = (h_inf - h) / tau_ref
  Both J_IP3R and J_RyR are multiplied by h_ref.
  This produces a refractory period after CICR: channels inactivate when
  ca_cicr rises, then recover with tau_ref (~665ms, Caya-Bissonnette 2023).

State variables (carry): IP3, Ca_ER, ca_cicr, h_ref, effcai, rho

Fixed biophysical constants:
  K_IP3 = 0.001 mM (1 µM, IP3R half-activation for IP3, from Li-Rinzel d1)
  K_RyR = 0.0005 mM (0.5 µM, RyR2 literature)
  gamma_ER = 0.1255 (ER/cytosol volume ratio, morphological)
  CA_ER_MAX = 2.0 mM (finite ER calcium capacity)

Free parameters: 21 (10 plasticity + 11 CICR/scaling)
  K_SERCA is now a free parameter (was fixed at 0.8934 mM).
  Neuronal SERCA2b has K_d ~ 0.5-10 µM, much lower than the Lytton 1992 value.
  tau_ref, K_h_ref: refractory gate parameters (Caya-Bissonnette 2023).
"""

import numpy as np
import jax
import jax.numpy as jnp
from cicr_common import CICRModel, compute_effcai_piecewise_linear_jax

# Fixed biophysical constants
CAI_REST = 70e-6       # Resting cytosolic calcium (mM)
TRACE_TAU_REF = 200.0  # Reference tau for threshold scaling

# IP3R activation (from Li-Rinzel d1, fixed)
K_IP3 = 0.001    # IP3R half-activation for IP3 (mM) = 1 µM

# RyR2 activation (literature values, fixed)
K_RYR = 0.0005   # RyR2 half-activation for Ca (mM) = 0.5 µM
N_RYR = 3         # Hill coefficient (cooperative binding)

# Compartment constants (fixed)
GAMMA_ER = 0.1255   # ER/cytosol volume ratio (morphological)

# Initial conditions and capacity
CA_ER_0 = 0.4       # Initial ER calcium (mM) — resting steady state
CA_ER_MAX = 2.0     # Maximum ER calcium (mM) — finite ER capacity


def _debug_sim(cai_1syn, t, cp, cq, is_apical, dp):
    """Numpy debug simulation — Simplified IP3R + RyR two-stage CICR with refractory gate."""
    delta_IP3 = float(dp["delta_IP3"]); tau_IP3 = float(dp["tau_IP3"])
    V_IP3R = float(dp["V_IP3R"]); V_RyR = float(dp["V_RyR"])
    V_SERCA = float(dp["V_SERCA"]); K_SERCA = float(dp["K_SERCA"])
    V_leak = float(dp["V_leak"])
    tau_ext = float(dp["tau_extrusion"]); tau_eff = float(dp["tau_eff"])
    tau_ref = float(dp.get("tau_ref", 665.0))
    K_h_ref = float(dp.get("K_h_ref", 0.0005))
    gamma_d = float(dp.get("gamma_d_GB_GluSynapse", 150.0))
    gamma_p = float(dp.get("gamma_p_GB_GluSynapse", 150.0))

    if is_apical:
        theta_d = dp["a20"]*cp + dp["a21"]*cq
        theta_p = dp["a30"]*cp + dp["a31"]*cq
    else:
        theta_d = dp["a00"]*cp + dp["a01"]*cq
        theta_p = dp["a10"]*cp + dp["a11"]*cq
    theta_p = max(theta_p, theta_d + 0.01)

    _nves_raw = dp.get("_nves", None)
    if _nves_raw is None or np.ndim(_nves_raw) == 0:
        nves = np.zeros(len(t))
    else:
        nves = np.asarray(_nves_raw, dtype=np.float64)

    n = len(cai_1syn); cai_total = np.copy(cai_1syn)
    ip3_out = np.zeros(n); ca_er_out = np.zeros(n)
    ca_cicr_out = np.zeros(n); ryr_open_out = np.zeros(n)
    h_ref_out = np.zeros(n)
    effcai_out = np.zeros(n); rho_out = np.zeros(n)

    # Initial state
    IP3 = 0.0; Ca_ER = CA_ER_0; ca_cicr = 0.0; h_ref = 1.0; effcai = 0.0
    rho = float(dp.get("rho0", 0.5))
    rho_out[0] = rho; ca_er_out[0] = CA_ER_0; h_ref_out[0] = 1.0

    for i in range(n - 1):
        dt = t[i+1] - t[i]
        if dt <= 0:
            ip3_out[i+1] = IP3; ca_er_out[i+1] = Ca_ER
            ca_cicr_out[i+1] = ca_cicr; ryr_open_out[i+1] = 0.0
            h_ref_out[i+1] = h_ref
            effcai_out[i+1] = effcai; rho_out[i+1] = rho
            cai_total[i+1] = cai_1syn[i+1] + ca_cicr
            continue

        ca_ext = max(0.0, cai_1syn[i] - CAI_REST)
        ca_i = ca_ext + ca_cicr

        # 1. IP3 dynamics (cnexp: dIP3/dt = -IP3/tau_IP3, with bolus)
        IP3 = IP3 + delta_IP3 * nves[i]
        IP3 = max(0.0, IP3 * np.exp(-dt / tau_IP3))

        # 2. Stage 1 — IP3R trigger (saturating activation, gated by h_ref)
        ip3_act = IP3 / (IP3 + K_IP3 + 1e-30)
        gradient = max(0.0, Ca_ER - ca_i)
        J_ip3r = V_IP3R * ip3_act * h_ref * gradient

        # 3. Stage 2 — RyR amplifier (Hill function of ca_cicr, gated by h_ref)
        cicr3 = ca_cicr ** N_RYR
        ryr_open = cicr3 / (K_RYR**N_RYR + cicr3 + 1e-30)
        J_ryr = V_RyR * ryr_open * h_ref * gradient

        # 4. SERCA pump (Hill n=2, driven by ca_ext ONLY, not ca_cicr)
        #    SERCA loads ER from VDCC calcium. Released ca_cicr gets extruded,
        #    not recycled — this gives single-shot refractory behavior.
        J_serca = V_SERCA * (ca_ext**2) / (K_SERCA**2 + ca_ext**2 + 1e-12)

        # 5. Passive ER leak
        J_leak = V_leak * gradient

        # 6. ER store dynamics (capped at CA_ER_MAX)
        dCa_ER = (J_serca - J_ip3r - J_ryr - J_leak) / GAMMA_ER
        Ca_ER = min(CA_ER_MAX, max(0.0, Ca_ER + dt * dCa_ER))

        # 7. Cytosolic CICR pool
        dca_cicr = J_ip3r + J_ryr + J_leak - J_serca - (ca_cicr / tau_ext)
        ca_cicr = max(0.0, ca_cicr + dt * dca_cicr)

        # 7b. Refractory inactivation gate (h_ref)
        h_inf = K_h_ref**2 / (K_h_ref**2 + ca_cicr**2 + 1e-30)
        dh = (h_inf - h_ref) / tau_ref
        h_ref = min(1.0, max(0.0, h_ref + dt * dh))

        # 8. Effective calcium integrator
        d_e = np.exp(-dt / tau_eff)
        effcai = effcai * d_e + (ca_ext + ca_cicr) * tau_eff * (1.0 - d_e)

        # 9. Plasticity (Graupner-Brunel)
        pot = 1.0 if effcai > theta_p else 0.0
        dep = 1.0 if effcai > theta_d else 0.0
        drho = (-rho*(1-rho)*(0.5-rho) + gamma_p*(1-rho)*pot - gamma_d*rho*dep) / 70000.0
        rho = min(1.0, max(0.0, rho + dt * drho))

        cai_total[i+1] = cai_1syn[i+1] + ca_cicr
        ip3_out[i+1] = IP3; ca_er_out[i+1] = Ca_ER
        ca_cicr_out[i+1] = ca_cicr; ryr_open_out[i+1] = ryr_open
        h_ref_out[i+1] = h_ref
        effcai_out[i+1] = effcai; rho_out[i+1] = rho

    return {"cai_total": cai_total, "priming": ip3_out, "ca_er": ca_er_out,
            "ca_cicr": ca_cicr_out, "ca_ryr": ryr_open_out, "ca_ip3r": ca_cicr_out,
            "h_ref": h_ref_out, "effcai": effcai_out, "rho": rho_out}


class SimpleCICRModel(CICRModel):
    DESCRIPTION = "Simplified IP3R Trigger + RyR Amplifier (no gating ODEs)"
    FIT_PARAMS = [
        # --- Plasticity Parameters (10) ---
        ("gamma_d_GB_GluSynapse", 50.0, 250.0), ("gamma_p_GB_GluSynapse", 50.0, 300.0),
        ("a00", 1.0, 5.0), ("a01", 1.0, 5.0), ("a10", 1.0, 5.0), ("a11", 1.0, 5.0),
        ("a20", 1.0, 5.0), ("a21", 1.0, 5.0), ("a30", 1.0, 10.0), ("a31", 1.0, 5.0),

        # --- CICR Parameters (11) ---
        # IP3 production (physiological: 0.5-50 µM per vesicle)
        ("delta_IP3", 0.0005, 0.05),      # IP3 bolus per vesicle (mM)
        ("tau_IP3", 200.0, 1000.0),        # IP3 decay time constant (ms)
        # Release conductances
        ("V_IP3R", 1e-5, 0.1),             # IP3R release conductance (/ms)
        ("V_RyR", 1e-5, 0.1),              # RyR release conductance (/ms)
        # SERCA pump
        ("V_SERCA", 1e-5, 0.1),            # SERCA max pump rate (mM/ms)
        ("K_SERCA", 0.0005, 0.01),         # SERCA half-activation (mM) = 0.5-10 µM
        # Leak
        ("V_leak", 6e-5, 2e-4),            # Passive ER leak (/ms) — sets priming duration
        # Cytosolic dynamics
        ("tau_extrusion", 10.0, 100.0),     # Cytosolic CICR clearance (ms) — Caya τ_c=23.5ms
        ("tau_effca_GB_GluSynapse", 50.0, 500.0),           # Effective calcium time constant (ms)
        # Refractory inactivation gate (Caya-Bissonnette 2023)
        ("tau_ref", 200.0, 1500.0),         # Refractory recovery time (ms) — Caya τ_r=665ms
        ("K_h_ref", 0.0001, 0.005),         # Inactivation half-point (mM) = 0.1-5 µM
    ]
    DEFAULT_PARAMS = {
        # Plasticity (from previous ReLU fit)
        "gamma_d_GB_GluSynapse": 244.27, "gamma_p_GB_GluSynapse": 201.89,
        "a00": 1.38, "a01": 2.65, "a10": 4.07, "a11": 4.20,
        "a20": 3.97, "a21": 1.47, "a30": 3.01, "a31": 2.44,
        # CICR (physiological defaults)
        "delta_IP3": 0.005, "tau_IP3": 330.0,    # 5 µM bolus
        "V_IP3R": 5e-4, "V_RyR": 2e-4,
        "V_SERCA": 0.001,
        "K_SERCA": 0.003,       # 3 µM — neuronal SERCA2b range
        "V_leak": 1e-4,         # tau_priming ~ 1.3s
        "tau_extrusion": 23.5, "tau_effca_GB_GluSynapse": 130.0,  # Caya τ_c = 23.5ms
        # Refractory gate (Caya-Bissonnette 2023)
        "tau_ref": 665.0,       # Refractory recovery (ms) — Caya τ_r = 665ms
        "K_h_ref": 0.0005,      # Inactivation half-point (mM) = 0.5 µM
    }

    def unpack_params(self, x):
        gb = {'gamma_d': x[0], 'gamma_p': x[1],
              'a00': x[2], 'a01': x[3], 'a10': x[4], 'a11': x[5],
              'a20': x[6], 'a21': x[7], 'a30': x[8], 'a31': x[9]}
        cicr = {'delta_IP3': x[10], 'tau_IP3': x[11],
                'V_IP3R': x[12], 'V_RyR': x[13],
                'V_SERCA': x[14], 'K_SERCA': x[15],
                'V_leak': x[16],
                'tau_extrusion': x[17], 'tau_effca_GB_GluSynapse': x[18],
                'tau_ref': x[19], 'K_h_ref': x[20]}
        return gb, cicr

    def get_init_fn(self):
        def init_fn(cai_first, rho0):
            # IP3, Ca_ER, ca_cicr, h_ref, effcai, rho
            return (0.0, CA_ER_0, 0.0, 1.0, 0.0, rho0)
        return init_fn

    def get_debug_sim_fn(self):
        return _debug_sim

    def get_step_factory(self):
        def scan_factory(params, syn_params):
            gb, cicr = params
            c_pre_opt, c_post_opt, is_apical, cai_pre, t_pre, cai_post, t_post = syn_params
            c_pre = jnp.max(compute_effcai_piecewise_linear_jax(cai_pre, t_pre, tau_effca=cicr['tau_effca_GB_GluSynapse']))
            c_post = jnp.max(compute_effcai_piecewise_linear_jax(cai_post, t_post, tau_effca=cicr['tau_effca_GB_GluSynapse']))

            theta_d = jnp.where(is_apical, gb['a20']*c_pre + gb['a21']*c_post, gb['a00']*c_pre + gb['a01']*c_post)
            theta_p = jnp.maximum(
                jnp.where(is_apical, gb['a30']*c_pre + gb['a31']*c_post, gb['a10']*c_pre + gb['a11']*c_post),
                theta_d + 0.01)

            # Pre-compute fixed RyR constants
            k_ryr_n = K_RYR ** N_RYR

            # Pre-compute fixed h_ref constants
            k_h2 = cicr['K_h_ref'] ** 2

            def scan_step(carry, inputs):
                IP3, Ca_ER, ca_cicr, h_ref, effcai, rho = carry
                cai_raw, dt, nves_i = inputs

                # External calcium (from VDCC/NMDA)
                ca_ext = jnp.maximum(0.0, cai_raw - CAI_REST)

                # Total cytosolic calcium
                ca_i = ca_ext + ca_cicr

                # 1. IP3 dynamics (cnexp decay + bolus)
                IP3_post = IP3 + cicr['delta_IP3'] * nves_i
                IP3_new = jnp.maximum(0.0, IP3_post * jnp.exp(-dt / cicr['tau_IP3']))

                # 2. ER-cytosol gradient (driving force for release)
                gradient = jnp.maximum(0.0, Ca_ER - ca_i)

                # 3. Stage 1 — IP3R trigger (saturating activation, gated by h_ref)
                ip3_act = IP3_new / (IP3_new + K_IP3 + 1e-30)
                J_ip3r = cicr['V_IP3R'] * ip3_act * h_ref * gradient

                # 4. Stage 2 — RyR amplifier (algebraic Hill, gated by h_ref)
                cicr_n = ca_cicr ** N_RYR
                ryr_open = cicr_n / (k_ryr_n + cicr_n + 1e-30)
                J_ryr = cicr['V_RyR'] * ryr_open * h_ref * gradient

                # 5. SERCA pump (driven by ca_ext ONLY — loads ER from VDCC Ca)
                #    Released ca_cicr gets extruded, not recycled to ER.
                J_serca = cicr['V_SERCA'] * (ca_ext**2) / (cicr['K_SERCA']**2 + ca_ext**2 + 1e-12)

                # 6. Passive ER leak (sets priming decay timescale)
                J_leak = cicr['V_leak'] * gradient

                # 7. ER store dynamics (capped at CA_ER_MAX)
                dCa_ER = (J_serca - J_ip3r - J_ryr - J_leak) / GAMMA_ER
                Ca_ER_new = jnp.clip(Ca_ER + dt * dCa_ER, 0.0, CA_ER_MAX)

                # 8. Cytosolic CICR pool
                dca_cicr = J_ip3r + J_ryr + J_leak - J_serca - (ca_cicr / cicr['tau_extrusion'])
                ca_cicr_new = jnp.maximum(0.0, ca_cicr + dt * dca_cicr)

                # 8b. Refractory inactivation gate
                h_inf = k_h2 / (k_h2 + ca_cicr_new**2 + 1e-30)
                dh = (h_inf - h_ref) / cicr['tau_ref']
                h_ref_new = jnp.where(dt > 0, jnp.clip(h_ref + dt * dh, 0.0, 1.0), h_ref)

                # 9. Effective calcium integrator
                d_eff = jnp.exp(-dt / cicr['tau_effca_GB_GluSynapse'])
                effcai_new = effcai * d_eff + (ca_ext + ca_cicr_new) * cicr['tau_effca_GB_GluSynapse'] * (1.0 - d_eff)

                # 10. Plasticity (Graupner-Brunel)
                pot = jnp.where(effcai_new > theta_p, 1.0, 0.0)
                dep = jnp.where(effcai_new > theta_d, 1.0, 0.0)
                drho = (-rho*(1.0-rho)*(0.5-rho) + gb['gamma_p']*(1.0-rho)*pot - gb['gamma_d']*rho*dep) / 70000.0
                rho_new = jnp.where(dt > 0, jnp.clip(rho + dt * drho, 0.0, 1.0), rho)

                return (IP3_new, Ca_ER_new, ca_cicr_new, h_ref_new, effcai_new, rho_new), None
            return scan_step
        return scan_factory

if __name__ == "__main__":
    SimpleCICRModel().run()
