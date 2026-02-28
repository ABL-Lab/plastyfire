#!/usr/bin/env python3
"""
Dual-Hill CICR variant (analog calcium amplitude tracking, 22 parameters).
Equivalent to fit_params_g_CICR.py but using the shared CICRModel base class.
"""

import numpy as np
from numba import njit
from cicr_common import CICRModel, _compute_effcai_batch, TAU_EFFCA_GB, MIN_CA_CR


# ══════════════════════════════════════════════════════════════════
# Dual-Hill CICR numba kernel
# ══════════════════════════════════════════════════════════════════

@njit(cache=True)
def _apply_cicr_batch_dual_hill(cai, t, c_pre, c_post,
                                 tau_IP3, K_P, n_prime, a_evt, b_evt,
                                 a_er, b_er, g_serca, g_release,
                                 sat_cap, K_Ca, m_ca, max_dt_sub=0.1):
    """Dual-Hill CICR with analog calcium amplitude tracking."""
    K_serca = 0.1
    fc = 0.83
    fe = 0.17
    caAvg = 0.0017
    MAX_DT_SUB = max_dt_sub

    n_syn = cai.shape[0]
    n_time = cai.shape[1]
    cai_total = np.copy(cai)

    for s in range(n_syn):
        ca_event_thresh_s = a_evt * c_pre[s] + b_evt * c_post[s]
        gamma_er = a_er * c_pre[s] + b_er * c_post[s]
        ca_er = (caAvg - fc * cai[s, 0]) / fe
        if ca_er < 0.0:
            ca_er = 0.0
        Priming_val = 0.0
        ca_cicr = 0.0
        was_above = 0
        event_peak_s = 0.0

        for i in range(n_time - 1):
            dt_full = t[i + 1] - t[i]
            if dt_full <= 0.0:
                cai_total[s, i + 1] = cai[s, i + 1] + ca_cicr
                continue

            # 1. Analog Calcium Amplitude Tracking
            is_above = 1 if cai[s, i] > ca_event_thresh_s else 0
            event_amp = 0.0
            if is_above == 1 and was_above == 0:
                event_peak_s = cai[s, i]
            elif is_above == 1:
                if cai[s, i] > event_peak_s:
                    event_peak_s = cai[s, i]
            elif was_above == 1:
                raw_amp = event_peak_s - ca_event_thresh_s
                event_amp = (raw_amp / sat_cap) if raw_amp < sat_cap else 1.0
                event_peak_s = 0.0
            was_above = is_above

            Priming_val = Priming_val * np.exp(-dt_full / tau_IP3) + event_amp

            # 2. Dual-Hill Gate Evaluation
            P_safe = max(0.0, Priming_val)
            P_gate = (P_safe**n_prime) / (K_P**n_prime + P_safe**n_prime + 1e-12)

            ca_raw_safe = max(0.0, cai[s, i])
            Ca_gate = (ca_raw_safe**m_ca) / (K_Ca**m_ca + ca_raw_safe**m_ca + 1e-12)

            # 3. Subcycled ER Dynamics
            n_sub = max(1, int(np.ceil(dt_full / MAX_DT_SUB)))
            dt_sub = dt_full / n_sub

            for _ in range(n_sub):
                ca_cyt = cai[s, i] + ca_cicr
                ca_cicr_um = 1000.0 * ca_cicr
                if ca_cicr > gamma_er:
                    J_serca_raw = g_serca * (ca_cicr_um**2) / (K_serca**2 + ca_cicr_um**2)
                    J_serca = min(J_serca_raw, ca_cicr / dt_sub)
                else:
                    J_serca = 0.0
                if ca_er > ca_cyt:
                    J_release = g_release * P_gate * Ca_gate * (ca_er - ca_cyt)
                else:
                    J_release = 0.0
                J_net = J_release - J_serca
                ca_er = ca_er - (fc / fe) * J_net * dt_sub
                if ca_er < 0.0: ca_er = 0.0
                ca_cicr = ca_cicr + J_net * dt_sub
                if ca_cicr < 0.0: ca_cicr = 0.0

            cai_total[s, i + 1] = cai[s, i + 1] + ca_cicr

    return cai_total


# ══════════════════════════════════════════════════════════════════
# DualHillCICRModel subclass
# ══════════════════════════════════════════════════════════════════

class DualHillCICRModel(CICRModel):
    DESCRIPTION = "Dual-Hill CICR"
    SLURM_SCRIPT_NAME = "submit_fitting_g_cicr.sh"
    RESULTS_SUFFIX = "_g_cicr"

    FIT_PARAMS = [
        ("gamma_d_GB_GluSynapse", 50.0, 200.0),
        ("gamma_p_GB_GluSynapse", 150.0, 300.0),
        ("a00", 1.0, 5.0), ("a01", 1.0, 5.0),
        ("a10", 1.0, 5.0), ("a11", 1.0, 5.0),
        ("a20", 1.0, 10.0), ("a21", 1.0, 5.0),
        ("a30", 1.0, 10.0), ("a31", 1.0, 5.0),
        ("tau_IP3", 100.0, 5000.0),
        ("K_P", 0.5, 10.0),
        ("n_prime", 1.0, 8.0),
        ("a_evt", 0.0, 0.05), ("b_evt", 0.0, 0.05),
        ("a_er", 0.0, 0.05), ("b_er", 0.0, 0.05),
        ("g_serca_cicr", 0.01, 50.0),
        ("g_release_cicr", 0.001, 10.0),
        ("sat_cap", 0.001, 0.5),
        ("K_Ca", 0.0001, 0.005),
        ("m_ca", 1.0, 8.0),
    ]

    DEFAULT_PARAMS = {
        "gamma_d_GB_GluSynapse": 101.5,
        "gamma_p_GB_GluSynapse": 216.2,
        "a00": 1.002, "a01": 1.954,
        "a10": 1.159, "a11": 2.483,
        "a20": 1.127, "a21": 2.456,
        "a30": 5.236, "a31": 1.782,
        "tau_IP3": 800.0,
        "K_P": 2.5,
        "n_prime": 3.0,
        "a_evt": 0.005, "b_evt": 0.005,
        "a_er": 0.005, "b_er": 0.005,
        "g_serca_cicr": 1.9565,
        "g_release_cicr": 0.5,
        "sat_cap": 0.05,
        "K_Ca": 0.0005,
        "m_ca": 3.0,
    }

    def apply_cicr(self, pd_item, x, max_dt_sub=0.1):
        return _apply_cicr_batch_dual_hill(
            pd_item["cai"], pd_item["t"], pd_item["c_pre"], pd_item["c_post"],
            x[10], x[11], x[12], x[13], x[14], x[15], x[16],
            x[17], x[18], x[19], x[20], x[21],
            max_dt_sub,
        )

    def unpack_gb_and_cicr(self, x):
        gb = {"gamma_d": x[0], "gamma_p": x[1],
              "a00": x[2], "a01": x[3], "a10": x[4], "a11": x[5],
              "a20": x[6], "a21": x[7], "a30": x[8], "a31": x[9]}
        cicr = {"tau_IP3": x[10], "K_P": x[11], "n_prime": x[12],
                "a_evt": x[13], "b_evt": x[14], "a_er": x[15], "b_er": x[16],
                "g_serca": x[17], "g_release": x[18], "sat_cap": x[19],
                "K_Ca": x[20], "m_ca": x[21]}
        return gb, cicr

    def warmup_cicr(self):
        _dummy_cai = np.zeros((1, 100), dtype=np.float64)
        _dummy_t = np.linspace(0, 100, 100, dtype=np.float64)
        _apply_cicr_batch_dual_hill(
            _dummy_cai, _dummy_t, np.zeros(1), np.zeros(1),
            800.0, 2.5, 3.0, 0.005, 0.005, 0.005, 0.005,
            1.9565, 0.5, 0.05, 0.0005, 3.0,
        )


if __name__ == "__main__":
    DualHillCICRModel().run()
