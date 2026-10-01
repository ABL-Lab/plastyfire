"""T15 gate: offline model_v2 dpre vs the NEURON GluSynapse_v2 mod on the toy section.

Records shaft cai from the NEURON run, feeds it and the spike arrivals to model_v2, compares final dpre.
Run at the NEURON dt (0.025 ms) and decimated to 0.25 ms (the extraction grid).
Drive 0 (shaft cai) with 30 ms clamp steps; drive 1 (spine VDCC influx) with 2 ms AP-like steps, fed
as -ica_VDCC step means (extract_v2.step_mean) on the coarse grids.

    python glusynapse_v2/tests/test_offline_gate.py
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, ".."))
import test_v2_single_synapse as S
import model_v2 as MV

CASES = {
    "mglu only":  (-40.0, -10.0, dict(A_mglu_GB=0.05)),
    "NO only":    (+10.0,   0.0, dict(A_NO_GB=500.0, tau_Z_GB=20.0)),
    "both":       (-10.0,   0.0, dict(A_mglu_GB=0.05, A_NO_GB=500.0, tau_Z_GB=20.0,
                                      theta_T_GB=2e-4, ca_scale_GB=1e-2)),
}

VDCC = dict(pre_drive_GB=1, i_scale_GB=1e-6, theta_Ti_GB=1e-7, theta_NOi_GB=1e-7)
CASES_VDCC = {
    "vdcc mglu":  (-10.0, 10.0, dict(VDCC, A_mglu_GB=0.05, tau_T_GB=20.0)),
    "vdcc NO":    (+5.0,  10.0, dict(VDCC, A_NO_GB=500.0, tau_NO_GB=5.0, tau_Z_GB=20.0)),
    "vdcc NO thZ": (+5.0, 10.0, dict(VDCC, A_NO_GB=2000.0, tau_NO_GB=5.0, tau_Z_GB=20.0, theta_Z_GB=0.3)),
    # v2.3: NO driven by effcai (point samples, as extract_v2), mGluR by VDCC
    "eff NO+mglu": (-5.0, 10.0, dict(VDCC, no_drive_GB=1, theta_NOe_GB=0.0, e_scale_GB=0.01, A_NO_GB=500.0,
                                     tau_NO_GB=5.0, tau_Z_GB=20.0, A_mglu_GB=0.05, tau_T_GB=20.0)),
    "eff NO thr": (+5.0, 10.0, dict(VDCC, no_drive_GB=1, theta_NOe_GB=2e-3, e_scale_GB=1e-3, A_NO_GB=500.0,
                                    tau_NO_GB=5.0, tau_Z_GB=20.0, theta_Z_GB=0.3)),
    # v2.4: NO driven by spine Ca, recovered offline from effcai (batch_v2.cacr_from_effcai), threshold on N
    "cacr NO thN": (+5.0, 10.0, dict(VDCC, no_drive_GB=2, theta_NOc_GB=1e-4, c_scale_GB=1e-3, theta_N_GB=2.0,
                                     A_NO_GB=60.0, tau_NO_GB=20.0, tau_Z_GB=40.0)),
}


def _signal(r, par):
    return -r["v2.ica_VDCC"] if par.get("pre_drive_GB") else r["shaft_cai"]


def _features(r, par, keep, n, tg, arr, P):
    """(tT, K) from the coarse signals: T from _signal, N from effcai when no_drive_GB is set."""
    tT = MV.feature_T(_coarse(_signal(r, par), keep, n, par)[None], tg, arr, P)
    if par.get("no_drive_GB") == 2:                     # as extraction: effcai point samples -> cacr
        # the toy runs the mod default tau_effca_GB (200 ms); extracted data use the cooker 278.3 ms
        # (constants.TAU_EFFCA_GB), which cacr_from_effcai defaults to
        import neuron
        from batch_v2 import cacr_from_effcai
        nsig = cacr_from_effcai(r["v2.effcai_GB"][keep][None], tg, tau=neuron.h.tau_effca_GB_GluSynapse_v2)[0]
    elif par.get("no_drive_GB"):
        nsig = r["v2.effcai_GB"][keep]
    else:
        nsig = _coarse(_signal(r, par), keep, n, par)
    return tT, MV.feature_K(nsig[None], tg, arr, P)


def _coarse(x, keep, n, par):
    if par.get("pre_drive_GB"):
        from extract_v2 import step_mean
        return step_mean(x, keep, np.append(keep[1:], n))
    return x[keep]


def main():
    S._snapshot_defaults()
    worst = 0.0
    for name, (dt_pp, v, par) in {**CASES, **CASES_VDCC}.items():
        pre, steps = S.pairing(dt_pp, v=v, dur=2.0 if par.get("pre_drive_GB") else 30.0)
        r = S.run(pre, steps, par)
        P = {k[:-3]: val for k, val in {**S.DEFAULTS, **par}.items()}
        P.update(ca_sh_rest=6.5e-5, ca_scale=par.get("ca_scale_GB", 1e-3))
        # mod trace in NEURON is recorded at t = k*dt starting 0; arrivals at pre + NetCon delay 0.1
        # (the mod's GLOBAL ca_scale default is set in the test only if passed)
        import neuron; h = neuron.h
        sig = _signal(r, par)
        for dec in (1, 10):
            keep = np.arange(0, len(sig), dec); dt = 0.025 * dec
            tg = MV.grid(len(keep), dt)
            arr = MV.arrival_index(pre, [0.1], tg)
            tT, K = _features(r, par, keep, len(sig), tg, arr, P)
            d = MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"])[0]
            ref = r["v2.dpre_GB"][-1]
            err = abs(d - ref); worst = max(worst, err / max(abs(ref), 1e-3)) if dec == 1 else worst
            print(f"{name:10s} dt {dt:5.3f}  mod {ref:+.5f}  offline {d:+.5f}  |err| {err:.2e}")
    # windowed grid (Ebner extraction): 0.25 ms in [spike - 50, spike + 150], 5 ms elsewhere
    sys.path.insert(0, os.path.join(HERE, "..", "..", "analytical_method"))
    from extract import window_index
    for name, (dt_pp, v, par) in {**CASES, **CASES_VDCC}.items():
        pre, steps = S.pairing(dt_pp, v=v, dur=2.0 if par.get("pre_drive_GB") else 30.0)
        r = S.run(pre, steps, par)
        P = {k[:-3]: val for k, val in {**S.DEFAULTS, **par}.items()}
        sig = _signal(r, par)
        keep = window_index(len(sig), 0.025, pre, 10, 200, 50.0, 150.0)
        tg = keep * 0.025
        arr = MV.arrival_index(pre, [0.1], tg)
        tT, K = _features(r, par, keep, len(sig), tg, arr, P)
        d = MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"])[0]
        ref = r["v2.dpre_GB"][-1]
        print(f"{name:10s} windowed ({len(keep)} pts)  mod {ref:+.5f}  offline {d:+.5f}  |err| {abs(d - ref):.2e}")
    print("PASS" if worst < 0.02 else "FAIL", f"worst rel err at 0.025 ms: {worst:.3%}")

if __name__ == "__main__":
    main()
