"""Fast unit tests for the v2 offline stack (no NEURON). Run after every change:

    python glusynapse_v2/tests/test_units.py            # all
    python glusynapse_v2/tests/test_units.py dpre       # tests whose name contains "dpre"

NEURON-level gates live in test_v2_single_synapse.py (mod behaviour) and test_offline_gate.py
(offline model vs mod).
"""
import glob, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(HERE); ROOT = os.path.dirname(V2)
sys.path.insert(0, V2); sys.path.insert(0, os.path.join(ROOT, "analytical_method"))
os.environ.setdefault("ANALYTICAL_BASIS_DIR", os.path.join(ROOT, "basis_results_edges_sabrina_n120_delta"))
import model_v2 as MV   # noqa: E402

MARKRAM = os.path.join(V2, "extracted", "markram_delta-cooker")
V1_EXTRACTED = os.path.join(ROOT, "analytical_method", "extracted_sabrina_n120_delta_d0p25")
PAIR = "180351-198084"


# ------------------------------------------------------------------ helpers
def synthetic(n_syn=3, dt=0.025, T=3000.0, seed=0):
    """Shaft Ca with bAP-like bumps (1-8 uM) and pre spikes at 5 Hz, some before and some after bumps."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    ca = np.full((n_syn, len(t)), 6.5e-5)
    pre = np.arange(100.0, T - 200, 200.0)
    for i in range(n_syn):
        for p in pre:
            on = p + rng.uniform(-40, 40)
            amp = rng.uniform(1e-3, 8e-3)
            m = (t >= on)
            ca[i, m] += amp * np.exp(-(t[m] - on) / 30.0) * (t[m] - on < 300)
    return t, ca, pre


def bruteforce(t, ca, arrivals, P):
    """Direct Euler integration of the mod's v2 equations (reference)."""
    P = {**MV.DEFAULTS, **P}
    dt = t[1] - t[0]
    n = ca.shape[0]
    T = np.zeros(n); N = np.zeros(n); Z = np.zeros(n); d = np.full(n, P["dpre0"])
    arr_idx = [set(np.searchsorted(t, a - 1e-9)) for a in arrivals]
    for k in range(len(t)):
        for i in range(n):
            if k in arr_idx[i]:
                d[i] -= P["A_mglu"] * np.tanh(T[i]) * (d[i] - P["dpre_min"])
                Z[i] += 1
        x = ca[:, k] - P["ca_sh_rest"]
        dd = P["A_NO"] * np.tanh(N) * np.tanh(np.maximum(Z - P["theta_Z"], 0)) * (P["dpre_max"] - d) / (1e3 * MV.TAU_IND)
        T += dt * (-T / P["tau_T"] + np.maximum(x - P["theta_T"], 0) / P["ca_scale"])
        N += dt * (-N / P["tau_NO"] + np.maximum(x - P["theta_NO"], 0) / P["ca_scale"])
        Z += dt * (-Z / P["tau_Z"])
        d += dt * dd
    return d


def offline(t, ca, pre, delay, P):
    P = {**MV.DEFAULTS, **P}
    arr = MV.arrival_index(pre, delay, t)
    tT, K = MV.pre_features(ca, t, arr, P)
    return MV.dpre_final(tT, K, P["A_mglu"], P["A_NO"], P["dpre_min"], P["dpre_max"], P["dpre0"])


# ------------------------------------------------------------------ model_v2
def test_dpre_vs_bruteforce():
    t, ca, pre = synthetic(dt=0.1, T=2000.0)
    delay = np.array([0.5, 1.0, 2.0])
    worst = 0.0
    for P in (dict(A_mglu=0.1, theta_T=1e-3, ca_scale=1e-2),
              dict(A_NO=300.0, theta_NO=2e-3, ca_scale=1e-2),
              dict(A_mglu=0.2, A_NO=300.0, theta_T=5e-4, theta_NO=3e-3, ca_scale=5e-3, tau_Z=15.0),
              dict(A_NO=300.0, theta_NO=2e-3, ca_scale=1e-2, theta_Z=0.4)):
        ref = bruteforce(t, ca, [pre + dl for dl in delay], P)
        got = offline(t, ca, pre, delay, P)
        err = np.max(np.abs(got - ref)); worst = max(worst, err)
        assert err < 2e-3 * max(1.0, np.max(np.abs(ref))), f"{P}: ref {ref} got {got}"
        assert np.max(np.abs(ref)) > 1e-3, f"{P}: test does nothing"
    return f"max |offline - brute force| {worst:.1e}"


def test_dpre_zero_amplitudes():
    t, ca, pre = synthetic()
    d = offline(t, ca, pre, [1.0, 1.0, 1.0], dict(dpre0=0.0))
    assert np.all(d == 0.0)
    d = offline(t, ca, pre, [1.0, 1.0, 1.0], dict(dpre0=0.3))
    assert np.allclose(d, 0.3)
    return "A = 0 keeps dpre at dpre0"


def test_dpre_bounds():
    t, ca, pre = synthetic()
    lo = offline(t, ca, pre, [1.0] * 3, dict(A_mglu=1.0, theta_T=0.0, ca_scale=1e-4))
    hi = offline(t, ca, pre, [1.0] * 3, dict(A_NO=1e6, theta_NO=0.0, ca_scale=1e-4))
    assert np.all(lo >= MV.DEFAULTS["dpre_min"] - 1e-12) and np.all(lo < -0.7), lo
    assert np.all(hi <= MV.DEFAULTS["dpre_max"] + 1e-12) and np.all(hi > 0.9), hi
    return f"saturates inside [{lo.min():.3f}, {hi.max():.3f}]"


def test_mglu_needs_post_before_pre():
    """Ca bump 20 ms BEFORE the pre spike depresses; bump 20 ms AFTER does not (tau_T short)."""
    dt = 0.025; t = np.arange(0, 2000, dt); pre = np.arange(100.0, 1900, 200.0)
    def run(offset):
        ca = np.full((1, len(t)), 6.5e-5)
        for p in pre:
            m = (t >= p + offset) & (t < p + offset + 10); ca[0, m] += 3e-3
        return offline(t, ca, pre, [0.0], dict(A_mglu=0.1, theta_T=5e-4, tau_T=15.0, ca_scale=1e-2))[0]
    before, after = run(-20.0), run(+20.0)
    assert before < -0.05 and abs(after) < 0.01 * abs(before), (before, after)
    return f"post-before-pre {before:+.3f}, pre-before-post {after:+.4f}"


def test_no_needs_pre_before_post():
    """NO pathway: pre 10 ms before a big Ca event potentiates much more than pre 50 ms after it,
    in the causal regime tau_Z > tau_NO (fit_v2 rule 4). With tau_Z < tau_NO and a saturated tanh(N)
    the order reverses, which is why the rule exists."""
    dt = 0.025; t = np.arange(0, 2000, dt); pre = np.arange(100.0, 1900, 200.0)
    def run(offset):
        ca = np.full((1, len(t)), 6.5e-5)
        for p in pre:
            m = (t >= p + offset) & (t < p + offset + 20); ca[0, m] += 5e-3
        return offline(t, ca, pre, [0.0], dict(A_NO=500.0, theta_NO=2e-3, tau_NO=5.0, tau_Z=30.0,
                                              ca_scale=1e-2))[0]
    pre_first, post_first = run(+10.0), run(-50.0)
    assert pre_first > 0.02 and pre_first > 5 * max(post_first, 1e-9), (pre_first, post_first)
    return f"pre-first {pre_first:+.3f}, post-first {post_first:+.4f}"


def _burst_signals(t, pre, offset, n_ap=3, isi=20.0):
    """A 3-AP burst per pairing starting at pre + offset: slow shaft Ca (tau 80 ms, as in the T5 traces)
    and the spine VDCC influx (-ica, a 1 ms pulse of 2e-6 nA per AP)."""
    ca = np.full((1, len(t)), 6.5e-5); vd = np.zeros((1, len(t)))
    for p in pre:
        for k in range(n_ap):
            a = p + offset + k * isi; m = t >= a
            ca[0, m] += 1e-3 * np.exp(-(t[m] - a) / 80.0)
            vd[0, (t >= a) & (t < a + 1.0)] += 2e-6
    return ca, vd


def test_vdcc_drive_keeps_timing():
    """Why pre_drive 1 exists: with the slow shaft Ca drive, NO overlap barely tells pre-before-burst from
    pre-after-burst; with the VDCC pulse drive it does, and the mGluR trace only sees post-before-pre."""
    dt = 0.025; t = np.arange(0, 2000, dt); pre = np.arange(100.0, 1900, 200.0)
    NO = dict(A_NO=500.0, tau_NO=5.0, tau_Z=30.0)
    shaft = dict(NO, theta_NO=2e-4, ca_scale=1e-3)
    vdcc = dict(NO, pre_drive=1, theta_NOi=1e-7, i_scale=1e-6)
    out = {}
    for off in (+10.0, -50.0):                                   # pre 10 ms before / 10 ms after the burst
        ca, vd = _burst_signals(t, pre, off)
        out[off] = offline(t, ca, pre, [0.0], shaft)[0], offline(t, vd, pre, [0.0], vdcc)[0]
    sel_shaft = out[10.0][0] / max(out[-50.0][0], 1e-9)
    sel_vdcc = out[10.0][1] / max(out[-50.0][1], 1e-9)
    # residual post-first overlap is the tau_NO = 5 ms tail of the last AP, 10 ms before the pre spike
    assert sel_shaft < 1.5 and sel_vdcc > 4 and out[10.0][1] > 0.02, (out, sel_shaft, sel_vdcc)
    # mGluR (T at arrival): only when the burst came first
    Pm = dict(pre_drive=1, theta_Ti=1e-7, i_scale=1e-6, tau_T=20.0)
    tt = {}
    for off in (+10.0, -50.0):
        _, vd = _burst_signals(t, pre, off)
        tt[off] = MV.feature_T(vd, t, MV.arrival_index(pre, [0.0], t), Pm).mean()
    assert tt[-50.0] > 0.3 and tt[10.0] < 1e-2, tt          # pre-first: tail of the burst 150 ms earlier
    return (f"NO pre-first/post-first: shaft {sel_shaft:.1f}x, vdcc {sel_vdcc:.0f}x; "
            f"tanh(T) at arrival post-first {tt[-50.0]:.2f}, pre-first {tt[10.0]:.1e}")


def test_theta_z_needs_high_pre_rate():
    """v2.2: with theta_Z the NO pathway needs a high presynaptic rate. Same pairing (5 pre spikes, each
    10 ms before a VDCC pulse), bursts every 1 s: at 50 Hz Z sums past theta_Z, at 10 Hz it does not.
    theta_Z = 0 gives NO at both rates (v2.1)."""
    dt = 0.025; t = np.arange(0, 6000, dt)
    def run(f, th):
        pre = np.concatenate([b + np.arange(5) * 1000.0 / f for b in np.arange(100.0, 5000, 1000.0)])
        vd = np.zeros((1, len(t)))
        for p in pre:
            vd[0, (t >= p + 10) & (t < p + 11)] += 2e-6
        return offline(t, vd, pre, [0.0], dict(pre_drive=1, theta_NOi=1e-7, i_scale=1e-6, A_NO=500.0,
                                               tau_NO=5.0, tau_Z=30.0, theta_Z=th))[0]
    hi, lo, lo0 = run(50.0, 1.5), run(10.0, 1.5), run(10.0, 0.0)
    assert hi > 0.01 and lo < 1e-6 < lo0, (hi, lo, lo0)
    return f"theta_Z 1.5: 50 Hz {hi:+.3f}, 10 Hz {lo:+.1e}; theta_Z 0 at 10 Hz {lo0:+.3f}"


def test_fit_tie_and_levels():
    """fit_v2 --tie copies the shared threshold; log-spaced levels hit both box ends."""
    import fit_v2 as F
    cfg = dict(model="v2", filters={"pre_drive": 1}, free_filters=["theta_Ti", "tau_NO"], filter_steps=5,
               tie={"theta_NOi": "theta_Ti"})
    for q in (0.0, 0.5, 1.0):
        a, P = F.unpack(np.r_[1., 1., 1., 1., -2., 1., q, q], cfg)
        assert P["theta_NOi"] == P["theta_Ti"], P
    lo, hi = F.FILTER_BOX["tau_NO"]
    assert np.isclose(F.unpack(np.r_[1., 1., 1., 1., -2., 1., 0, 0], cfg)[1]["tau_NO"], lo)
    assert np.isclose(F.unpack(np.r_[1., 1., 1., 1., -2., 1., 1, 1], cfg)[1]["tau_NO"], hi)
    return "theta_NOi tied to theta_Ti; log levels span the box"


def test_rho_kernel_numba():
    """numba rho_full_all (serial and parallel) == per-record numpy rho_full to 1e-12, also from non-binary
    rho0 (no fixed-point skipping), and faster than the old stacked numpy loop."""
    from batch_v2 import BatchV2
    B = BatchV2([MARKRAM], protocols=["10Hz_-10ms", "10Hz_10ms", "10Hz_30ms"], pairs={PAIR}, verbose=False)
    a = dict(a00=1.5, a01=1.3, a10=1.8, a11=2.2); a.update(a20=a["a00"], a21=a["a01"], a30=a["a10"], a31=a["a11"])
    th = [B.thetas(r, a) for r in B.recs]
    td = np.concatenate([t[0] for t in th]); tp = np.concatenate([t[1] for t in th])
    ref = np.concatenate([B.rho_full(r, *t) for r, t in zip(B.recs, th)])
    B._stack(); B.rho_full_all(td, tp)                                   # compile
    t0 = time.time(); got = B.rho_full_all(td, tp); t_nb = time.time() - t0
    par = B.rho_full_all(td, tp, parallel=True)
    assert np.abs(got - ref).max() < 1e-12 and np.abs(par - ref).max() < 1e-12, np.abs(got - ref).max()
    tp_low = td * 0.5                                                   # theta_p < theta_d (rule 1 broken)
    ref_low = np.concatenate([B.rho_full(r, t[0], t[0] * 0.5) for r, t in zip(B.recs, th)])
    assert np.abs(B.rho_full_all(td, tp_low) - ref_low).max() < 1e-12, "theta_p < theta_d"
    r0 = B._rho0.copy(); B._rho0[:] = 0.37
    for r in B.recs: r["rho0"] = np.full(len(r["syn"]), 0.37)
    ref2 = np.concatenate([B.rho_full(r, *t) for r, t in zip(B.recs, th)])
    got2 = B.rho_full_all(td, tp); B._rho0 = r0
    assert np.abs(got2 - ref2).max() < 1e-12, np.abs(got2 - ref2).max()
    E, H = B._E.T, B._H[B._rec].T                                        # old layout (T, N)
    t0 = time.time(); rho = r0.copy()
    for k in range(E.shape[0]):
        e = E[k]; dep = e > td; pot = e > tp
        rho += H[k] * (-rho * (1 - rho) * (0.5 - rho) + pot * 101.5 * (1 - rho) - dep * (1 - pot) * 199.77 * rho)
        np.clip(rho, 0.0, 1.0, out=rho)
    t_np = time.time() - t0
    return f"|numba - numpy| {np.abs(got - ref).max():.0e}; {len(td)} syn: numpy loop {t_np:.2f}s, numba {t_nb:.3f}s ({t_np / t_nb:.0f}x)"


def test_bap_gate():
    """V2_BAP_GATE: synapses with c_post below the cutoff keep rho0 (rho_full, rho_fast, numba), the rest are
    unchanged vs the ungated rule; off (None) is the Chindemi rule."""
    import batch_v2
    from batch_v2 import BatchV2
    B = BatchV2([MARKRAM], protocols=["10Hz_10ms", "10Hz_-10ms"], pairs={PAIR}, verbose=False)
    a = dict(a00=1.5, a01=1.3, a10=1.8, a11=2.2); a.update(a20=a["a00"], a21=a["a01"], a30=a["a10"], a31=a["a11"])
    cp = np.concatenate([r["c_post"] for r in B.recs])
    cut = float(np.median(cp))                                         # gate about half the synapses
    off = cp < cut
    assert 0 < off.sum() < len(off)
    free = [np.concatenate([f(r, *B.thetas(r, a)) for r in B.recs]) for f in (B.rho_full, B.rho_fast)]
    batch_v2.BAP_GATE = cut
    try:
        th = [B.thetas(r, a) for r in B.recs]
        gated = [np.concatenate([f(r, *t) for r, t in zip(B.recs, th)]) for f in (B.rho_full, B.rho_fast)]
        td = np.concatenate([t[0] for t in th]); tp = np.concatenate([t[1] for t in th])
        B._stack(); nb = B.rho_full_all(td, tp)
    finally:
        batch_v2.BAP_GATE = None
    rho0 = np.concatenate([r["rho0"] for r in B.recs])
    for g, f in zip(gated, free):
        assert np.array_equal(g[off], rho0[off]), "gated synapses moved"
        assert np.abs(g[~off] - f[~off]).max() < 1e-12, "ungated synapses changed"
    assert np.abs(nb - gated[0]).max() < 1e-12, "numba != numpy with inf thresholds"
    moved = (np.abs(free[0] - rho0) > 1e-6)
    return f"{off.sum()}/{len(off)} gated; {moved[off].sum()} of them would have moved ungated"


def test_step_mean_preserves_charge():
    """extract_v2.step_mean keeps the integral of a 1 ms pulse on the 0.25 ms and windowed grids."""
    from extract_v2 import step_mean
    dt = 0.025; t = np.arange(0, 400, dt); x = ((t >= 100.3) & (t < 101.3)) * 2e-6
    for keep in (np.arange(0, len(t), 10), np.r_[np.arange(0, 3000, 200), np.arange(3000, 6000, 10),
                                                 np.arange(6000, len(t), 200)]):
        ends = np.append(keep[1:], len(t))
        m = step_mean(x, keep, ends).astype(np.float64)
        assert abs(np.sum(m * (ends - keep)) * dt - x.sum() * dt) < 1e-12
    return "charge preserved on uniform and windowed grids"


def test_drive0_unchanged():
    """pre_drive 0 (default) gives the same features as before the switch existed; drive 1 without a
    vdcc signal is refused by BatchV2."""
    t, ca, pre = synthetic(); delay = np.zeros(len(ca))
    arr = MV.arrival_index(pre, delay, t)
    a = MV.pre_features(ca, t, arr, dict(theta_T=2e-4))
    b = MV.pre_features(ca, t, arr, dict(theta_T=2e-4, pre_drive=0, i_scale=5.0, theta_Ti=9.0))
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    from batch_v2 import _sig
    try:
        _sig(dict(pair="p", proto="q", shaft_cai=ca, vdcc=None), dict(pre_drive=1)); raise AssertionError
    except ValueError:
        pass
    return "drive-1 keys ignored under drive 0; missing vdcc refused"


def test_windowed_grid():
    """Windowed grid (fine near spikes, coarse in quiet gaps) gives the uniform-grid answer."""
    from extract import window_index
    dt = 0.025; t = np.arange(0, 12000, dt)
    pre = np.array([100.0, 150.0, 5100.0, 5150.0, 10100.0])
    ca = np.full((2, len(t)), 6.5e-5)
    for p in pre:
        m = (t >= p - 10) & (t < p + 200); ca[:, m] += 4e-3 * np.exp(-(t[m] - p + 10) / 40.0)
    P = dict(A_mglu=0.1, A_NO=300.0, theta_T=5e-4, theta_NO=2e-3, ca_scale=1e-2)
    uni = offline(t[::10], ca[:, ::10], pre, [0.5, 1.5], P)
    keep = window_index(len(t), dt, pre, 10, 200, 50.0, 1000.0)
    win = offline(t[keep], ca[:, keep], pre, [0.5, 1.5], P)
    assert len(keep) < len(t) // 10 and np.allclose(uni, win, atol=1e-6), (uni, win, len(keep))
    return f"{len(keep)} vs {len(t)//10} samples, max diff {np.abs(uni - win).max():.1e}"


def test_arrival_index():
    t = np.arange(0, 10, 0.25)
    a = MV.arrival_index([1.0, 2.1], [0.0, 0.5], t)
    assert a.tolist() == [[4, 9], [6, 11]], a        # 1.0 -> 1.0; 2.1 -> 2.25; 1.5; 2.6 -> 2.75
    return "first grid point at or after prespike + delay"


# ------------------------------------------------------------------ readout
def test_pairbasis_matches_v1():
    import model as M
    from batch_v2 import PairBasis
    d = np.load(os.path.join(MARKRAM, f"{PAIR}__10Hz_-10ms.npz"))
    pb = PairBasis(PAIR, d["syn"]); df = M.load_basis(*PAIR.split("-"))
    rng = np.random.default_rng(0)
    for _ in range(20):
        r0 = rng.integers(0, 2, len(d["syn"])).astype(float); r1 = rng.random(len(d["syn"]))
        m1, s1 = M.epsp_from_basis(df, r1); m2, s2 = pb.epsp(r1, np.zeros(len(r1)))
        assert abs(m1 - m2) < 1e-9 and abs(s1 - s2) < 1e-9
        assert abs(M.epsp_ratio(df, r0, r1) - pb.ratio(r0, r1, np.zeros(len(r1)))) < 1e-9
    return "PairBasis == model.epsp_from_basis / epsp_ratio at dpre = 0 (20 random states)"


def test_readout_use_scaling():
    from batch_v2 import PairBasis
    d = np.load(os.path.join(MARKRAM, f"{PAIR}__10Hz_-10ms.npz"))
    pb = PairBasis(PAIR, d["syn"]); n = len(d["syn"]); z = np.zeros(n)
    e0, _ = pb.epsp(z, np.zeros(n))
    half, _ = pb.epsp(z, np.full(n, -0.5))
    assert abs(half / e0 - 0.5) < 1e-9, half / e0              # linear in Use while Use < 1
    big, _ = pb.epsp(np.ones(n), np.full(n, 5.0))               # Use capped at 1
    cap = np.sum((pb.e0 * pb.w + pb.delta) / pb.Use_p)
    assert abs(big - cap) < 1e-9
    return "Use ratio scaling linear; capped at Use = 1"


# ------------------------------------------------------------------ targets
def test_targets():
    import csv
    from targets import load_targets, EBNER_CSV
    t = load_targets()
    assert all(v[1] > 0 for v in t.values()), "every target needs an SEM"
    assert ("nevian_3ap_50hz_dt-50ms", "mglu_block") in t and ("10Hz_-10ms", "control") in t
    m, s, n, _ = t[("nevian_3ap_50hz_dt-50ms", "mglu_block")]
    w = np.array([1 / .16 ** 2, 1 / .14 ** 2, 1 / .07 ** 2])
    assert abs(m - (w @ [1.06, 1.09, 1.11]) / w.sum()) < 1e-9 and n == 16, (m, n)   # MCPG+AM251+U73122
    if os.path.isfile(EBNER_CSV):
        for r in csv.DictReader(open(EBNER_CSV)):
            if not r["sem"]:
                continue
            loc = r["synapse_location"].split()[0]
            k = (r["protocol_id"] if loc == "all" else f'{r["protocol_id"]}@{loc}', "control")
            assert abs(t[k][0] - float(r["mean_ratio"])) < 1e-9, k
    p = load_targets(("paired_l5",))
    assert not any(k[0].startswith(("nevian", "letzkus")) or "@" in k[0] for k in p), "paired_l5: L5-L5 unitary only"
    assert ("10Hz_-50ms", "control") not in p, "tail anchors are not data"
    assert p[("sjostrom_20hz_dt-10ms", "mglu_block")][:3] == (1.02, 0.07, 7)
    assert {c for _, c in p} == {"control", "mglu_block", "nmdar_block"}
    assert ("sjostrom_burst5x20hz_r50_dt-200ms", "control") in p and len(p) == 24
    q = load_targets(("paired_l5", "sjostrom07"))
    assert len(q) == 29 and q[("sjostrom07_step200ms_pair", "no_block")][2] == 14
    assert abs(q[("sjostrom07_step200ms_pair", "no_block")][0] - 1.36) < 0.01, "text: pooled NO block 136%"
    assert {c for _, c in q} == {"control", "mglu_block", "nmdar_block", "no_block"}
    return f"{len(t)} targets, csv consistent, pharmacology combined; paired_l5 {len(p)}, + sjostrom07 {len(q)}"


def test_targets_l23():
    from targets import load_targets, EBNER_CSV
    q = load_targets(("paired_l23l5",))
    assert not any(k[0].startswith(("nevian", "markram")) or k[0].startswith("10Hz") for k in q), "unitary L2/3->L5 only"
    assert q[("sjostrom_50hz_dt+10ms", "control")][:3] == (1.06, 0.09, 19), "L2/3->L5 value, not the L5-L5 1.57"
    assert {c for _, c in q} == {"control", "nmdar_block"} and all(v[1] > 0 for v in q.values())
    if os.path.isfile(EBNER_CSV):
        assert q[("letzkus_1ap_dt+10ms", "control")][:3] == (0.72, 0.03, 15)      # '@pooled' dropped
        assert q[("letzkus_3ap_200hz_dt-10ms@distal", "control")][:3] == (1.42, 0.09, 7)
        assert q[("sjostrom_50hz_dt+10ms@distal", "control")][:3] == (0.86, 0.09, 8)
        assert len(q) == 9, len(q)
    assert load_targets(("paired_l23l23",)) == {}
    for g in (("paired_l5", "paired_l23l5"), ("paired_l23l5", "paired_l23l23")):
        try:
            load_targets(g); raise AssertionError(f"{g} must not combine")
        except ValueError:
            pass
    assert len(load_targets(("paired_l5",))) == 24, "paired_l5 unchanged"
    return f"paired_l23l5 {len(q)}, paired_l23l23 0, pathway groups exclusive"


# ------------------------------------------------------------------ batch
def _tmp_v1_dir(pairs, protos):
    import tempfile
    d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR", "/tmp"))
    for p in pairs:
        for pr in protos:
            f = os.path.join(V1_EXTRACTED, f"{p}__{pr}.npz")
            if os.path.isfile(f):
                os.symlink(f, os.path.join(d, os.path.basename(f)))
    return d


def test_batch_v1_equivalence():
    import json, batch as BA
    from batch_v2 import BatchV2
    a = json.load(open(os.path.join(ROOT, "fit_results", "delta-cooker.json")))["params"]
    pairs = sorted({os.path.basename(f).split("__")[0] for f in glob.glob(os.path.join(V1_EXTRACTED, "*.npz"))})[:3]
    d = _tmp_v1_dir(pairs, ["10Hz_-10ms", "10Hz_10ms"])
    B1 = BA.Batch(extracted=d); B2 = BatchV2([d], verbose=False)
    r1 = B1.rho_fast(*B1.thetas(a))
    r2 = np.concatenate([B2.rho_fast(r, *B2.thetas(r, a)) for r in B2.recs])
    f1 = B1.rho_full(*B1.thetas(a))
    f2 = np.concatenate([B2.rho_full(r, *B2.thetas(r, a)) for r in B2.recs])
    assert np.abs(r1 - r2).max() < 1e-4 and np.abs(f1 - f2).max() < 1e-4
    fa = B2.rho_full_all(np.concatenate([B2.thetas(r, a)[0] for r in B2.recs]),
                         np.concatenate([B2.thetas(r, a)[1] for r in B2.recs]))
    assert np.abs(fa - f2).max() < 1e-6, f"rho_full_all vs per-record: {np.abs(fa - f2).max()}"   # float32 steps
    _, df1 = B1.curve(r1); df2 = B2.record_ratios(a, {})
    assert np.allclose(np.sort(df1["ratio"].values), np.sort(df2["ratio"].values), atol=1e-4)
    return f"rho fast {np.abs(r1-r2).max():.1e}, full {np.abs(f1-f2).max():.1e}, ratios equal"


def test_batch_conditions():
    import json
    from batch_v2 import BatchV2
    a = json.load(open(os.path.join(ROOT, "fit_results", "delta-cooker.json")))["params"]
    B = BatchV2([MARKRAM], protocols=["10Hz_-10ms"], pairs={PAIR}, verbose=False)
    P = dict(A_mglu=0.05, A_NO=100.0, theta_NO=2e-3)
    df = B.record_ratios(a, P, ("control", "mglu_block", "post_nmdar", "nmdar_block")).set_index("condition")
    r0 = B.recs[0]; same = B.basis(r0).ratio(r0["rho0"], r0["rho0"], np.zeros(len(r0["syn"])))
    assert abs(df.loc["nmdar_block", "ratio"] - same) < 1e-12 and df.loc["nmdar_block", "dpre"] == 0.0   # APV: no change
    only_no = B.record_ratios(a, dict(P, A_mglu=0.0), ("control",)).set_index("condition")
    assert abs(df.loc["mglu_block", "ratio"] - only_no.loc["control", "ratio"]) < 1e-12
    assert abs(df.loc["post_nmdar", "rho"] - B.recs[0]["rho0"].mean()) < 1e-12
    assert df.loc["control", "dpre"] < df.loc["mglu_block", "dpre"]
    from targets import load_targets
    T = {k: v for k, v in load_targets().items() if k[0].startswith("letzkus")}
    B2 = BatchV2([MARKRAM], protocols=["10Hz_-10ms"], pairs={PAIR}, verbose=False)
    for r in B2.recs:
        r["proto"] = "letzkus_3ap_200hz_dt+10ms"          # relabel one record to exercise @loc handling
    _, res = B2.objective(a, {}, T)
    assert list(res["target"]) == ["letzkus_3ap_200hz_dt+10ms@proximal"], list(res["target"])
    none = B.record_ratios(a, {}, ("control",)).iloc[0]
    assert none["dpre"] == 0.0
    return "mglu_block = A_mglu 0; post_nmdar freezes rho; nmdar_block = no change; A = 0 -> dpre 0"


def test_no_drive_effcai():
    """v2.3 (no_drive 1): NO is driven by effcai above theta_NOe, so a burst-sized effcai rise gives NO and a
    single-pairing rise does not; the T pathway keeps its own signal; post_nmdar switches NO off."""
    import json
    from batch_v2 import BatchV2, _sig
    dt = 0.025; t = np.arange(0, 1000, dt); pre = np.array([100.0]); arr = MV.arrival_index(pre, [0.0], t)
    base = 0.01
    def eff(rise):                                           # effcai-like: step up after the pre spike, 200 ms decay
        e = np.full(len(t), base); m = t >= 105
        e[m] += rise * np.exp(-(t[m] - 105) / 200.0); return e[None, :]
    P = dict(no_drive=1, theta_NOe=base + 0.05, e_scale=0.01, tau_NO=10.0, tau_Z=60.0, theta_Z=0.0)
    k_burst = MV.feature_K(eff(0.08), t, arr, P).sum(); k_single = MV.feature_K(eff(0.04), t, arr, P).sum()
    assert k_single == 0.0 and k_burst > 1.0, (k_single, k_burst)
    r = dict(effcai=eff(0.08).astype(np.float32), vdcc=np.zeros((1, len(t))), shaft_cai=None)
    assert _sig(r, dict(no_drive=1, pre_drive=1), "NO") is r["effcai"]
    assert _sig(r, dict(no_drive=1, pre_drive=1), "T") is r["vdcc"]
    assert _sig(r, dict(no_drive=0, pre_drive=1), "NO") is r["vdcc"]
    a = json.load(open(os.path.join(ROOT, "fit_results", "delta-cooker.json")))["params"]
    B = BatchV2([MARKRAM], protocols=["10Hz_10ms"], pairs={PAIR}, verbose=False)
    P2 = dict(A_mglu=0.05, A_NO=300.0, no_drive=1, theta_NOe=0.02, theta_Z=0.0, tau_Z=60.0)
    df = B.record_ratios(a, P2, ("control", "post_nmdar")).set_index("condition")
    off = B.record_ratios(a, dict(P2, A_NO=0.0), ("control",)).set_index("condition")
    assert abs(df.loc["post_nmdar", "dpre"] - off.loc["control", "dpre"]) < 1e-12
    assert df.loc["control", "dpre"] > df.loc["post_nmdar", "dpre"] + 0.01, df
    old = B.record_ratios(a, dict(P2, no_drive=0, pre_drive=0, theta_NO=2e-4), ("control", "post_nmdar"))
    assert abs(old.dpre.iloc[0] - old.dpre.iloc[1]) < 1e-12              # no_drive 0: unchanged behaviour
    return (f"K burst {k_burst:.1f} ms vs single {k_single}; post_nmdar dpre {df.loc['post_nmdar', 'dpre']:+.3f} "
            f"(= A_NO 0) vs control {df.loc['control', 'dpre']:+.3f}")


def test_no_block():
    """T30 no_block (L-NAME / cPTIO): dpre and ratio equal control with A_NO = 0; rho and the eCB term are kept."""
    import json
    from batch_v2 import BatchV2
    a = json.load(open(os.path.join(ROOT, "fit_results", "delta-cooker.json")))["params"]
    B = BatchV2([MARKRAM], protocols=["10Hz_10ms"], pairs={PAIR}, verbose=False)
    P2 = dict(A_mglu=0.05, A_NO=300.0, no_drive=1, theta_NOe=0.02, theta_Z=0.0, tau_Z=60.0)
    df = B.record_ratios(a, P2, ("control", "no_block")).set_index("condition")
    off = B.record_ratios(a, dict(P2, A_NO=0.0), ("control",)).set_index("condition")
    for k in ("dpre", "ratio", "rho"):
        assert abs(df.loc["no_block", k] - off.loc["control", k]) < 1e-12, k
    assert df.loc["control", "dpre"] > df.loc["no_block", "dpre"] + 0.01, df
    return f"no_block dpre {df.loc['no_block', 'dpre']:+.3f} (= A_NO 0) vs control {df.loc['control', 'dpre']:+.3f}"


def test_t_drive_ecb():
    """T25 (t_drive 1): T is driven by effcai / c_post above theta_Te, so post-before-pre activity (a bAP
    burst) depresses the next pre spike and a pre spike before the burst is untouched; a rise below theta_Te
    gives nothing; the depression saturates at the fitted floor dpre_min."""
    from batch_v2 import _sig
    dt = 0.025; t = np.arange(0, 1000, dt)
    def eff(rise, t0=100.0):                                 # effcai / c_post: step up at t0, 278 ms decay
        e = np.zeros(len(t)); m = t >= t0; e[m] = rise * np.exp(-(t[m] - t0) / 278.0); return e[None, :]
    pre = np.array([50.0, 150.0, 600.0]); arr = MV.arrival_index(pre, [0.0], t)
    P = dict(t_drive=1, theta_Te=0.5, tau_T=50.0)
    tT = MV.feature_T(eff(1.5), t, arr, P)[0]
    assert tT[0] == 0.0 and tT[1] > 0.5 and 0.0 < tT[2] < tT[1], tT        # pre-before: 0; after: depressed
    assert MV.feature_T(eff(0.4), t, arr, P).max() == 0.0                    # below theta_Te
    assert np.all(MV.feature_T(eff(1.5), t, arr, dict(P, t_drive=0))[0] == MV.feature_T(eff(1.5), t, arr,
                  dict(P, t_drive=0, theta_Te=9.0))[0])                      # t_drive 0 ignores theta_Te
    r = dict(effcai=np.full((2, 5), 0.2, np.float32), c_post=np.array([0.1, 0.4]), vdcc=np.ones((2, 5)))
    x = _sig(r, dict(t_drive=1, pre_drive=1), "T")
    assert np.allclose(x[:, 0], [2.0, 0.5]) and _sig(r, dict(t_drive=1, pre_drive=1), "NO") is r["vdcc"]
    d = MV.dpre_final(np.ones((1, 40)), np.zeros((1, 41)), 0.3, 0.0, dpre_min=-0.35)
    assert abs(d[0] + 0.35) < 1e-6, d                                       # floor = dpre_min
    return f"tanh(T) at pre spikes {np.round(tT, 3)}; floor {d[0]:+.3f}"


def test_t_drive_post_aps():
    """T25 (t_drive 2): cell-level eCB from post APs (S +1 per AP, tau_E1), supralinear through theta_Te, gated
    by theta_Tg. Reproduces the Sjostrom 2003 window: LTD for one AP 10-50 ms before the pre spike, none at
    100-200 ms or for pre-before-post, but LTD for a 5-AP 20 Hz burst ending 120-200 ms before."""
    from batch_v2 import _sig, post_impulses
    dt = 0.1; t = np.arange(0, 1000, dt)
    P = dict(t_drive=2, tau_E1=100.0, theta_Te=0.8, tau_T=40.0, theta_Tg=0.5)
    def gate(posts, pre):
        r = dict(postspikes=np.array(posts), t=t, syn=np.arange(3))
        sig = _sig(r, P, "T"); assert sig.shape == (3, len(t)) and sig.sum() == 3 * len(posts)
        g = MV.feature_T(sig, t, MV.arrival_index(np.array([pre]), np.zeros(3), t), P)
        assert np.all(g == g[0]), "cell-level: same gate at every synapse"
        return float(g[0, 0])
    one = lambda d: gate([100.0], 100.0 - d)
    burst = lambda d: gate([100.0, 150.0, 200.0, 250.0, 300.0], 300.0 - d)
    ltd = {"-10": one(-10), "-25": one(-25), "-50": one(-50), "B-120": burst(-120), "B-200": burst(-200)}
    none = {"+10": one(10), "+25": one(25), "-100": one(-100), "-120": one(-120), "-200": one(-200)}
    assert min(ltd.values()) > 0.05 and max(none.values()) == 0.0, (ltd, none)
    assert np.allclose(post_impulses([0.05, 0.05], np.array([0.0, 0.1, 0.2])), [0, 2, 0])
    g0 = MV.feature_T(np.ones((2, len(t))), t, np.array([[500], [900]]), dict(tau_T=40.0))
    assert np.all(g0 == np.tanh(MV._lowpass_grid(np.maximum(np.ones((2, len(t))) - MV.DEFAULTS["ca_sh_rest"]
                  - MV.DEFAULTS["theta_T"], 0) / MV.DEFAULTS["ca_scale"], t, 40.0)[[0, 1], [500, 900]][:, None]))
    return "LTD gates " + " ".join(f"{k} {v:.2f}" for k, v in ltd.items()) + "; none: all 0"


def test_t_drive3_veto():
    """v_events / veto: one event per excursion (hysteresis), a slow drift is not an event, and events within
    w ms after the synapse's own arrival are dropped (boundary: dropped inside [arr, arr + w), kept at arr + w)."""
    dt = 0.1; t = np.arange(0, 2000, dt); v = np.full(len(t), -70.0)
    v += 10.0 * t / 2000.0                                   # 10 mV drift over 2 s: v_bar follows, no event
    for t0 in (300.0, 305.0, 900.0):                         # 300 and 305 ms bumps merge into one excursion
        v += 8.0 * np.exp(-np.maximum(t - t0, 0) / 3.0) * (t >= t0)
    ev = MV.v_events(v, dt)
    assert len(ev) == 2 and abs(ev[0] - 300.0) < 1.0 and abs(ev[1] - 900.0) < 1.0, ev   # 300+305 merge (hysteresis)
    quiet = MV.v_events(-70.0 + 10.0 * t / 2000.0, dt); assert len(quiet) == 0
    e = np.array([100.0, 110.0, 114.9, 115.0, 200.0, np.nan])
    k = MV.veto(e, [100.0], 15.0)
    assert np.allclose(k, [115.0, 200.0]), k               # 100, 110, 114.9 inside [100, 115) dropped; 115 kept
    assert np.allclose(MV.veto(e, [], 15.0), e[:5]) and len(MV.veto([np.nan], [1.0], 15.0)) == 0
    return "events %s; veto keeps %s" % (np.round(ev, 1).tolist(), k.tolist())


def _vrec(vev, postspikes=(), prespikes=(300.0,), n=3, dt=0.1, T=1000.0):
    t = np.arange(0, T, dt)
    return dict(pair="p", proto="q", syn=np.arange(n), t=t, vev=np.array(vev, float), prespikes=np.array(prespikes),
                postspikes=np.array(postspikes), delay=np.zeros(n))


def test_t_drive3_window():
    """t_drive 3 on synthetic local V events reproduces the Sjostrom 2003 window like test_t_drive_post_aps: the
    bAP of each post AP is a local event, the pre EPSP (own arrival) is vetoed."""
    from batch_v2 import _sig
    P = dict(t_drive=3, tau_E1=100.0, theta_Te=0.8, tau_T=40.0, theta_Tg=0.5)
    def gate(posts, pre):
        # event at each post AP (+1 ms bAP delay) and at the pre arrival + 1 ms (EPSP, vetoed)
        ev = np.array(sorted([p + 1.0 for p in posts] + [pre + 1.0]))
        r = _vrec([ev] * 3, prespikes=[pre]); sig = _sig(r, P, "T")
        assert sig.shape == (3, len(r["t"])) and sig.sum() in (0, 3 * len(posts)), sig.sum()   # +10: bAP inside w
        g = MV.feature_T(sig, r["t"], MV.arrival_index(np.array([pre]), np.zeros(3), r["t"]), P)
        return float(g[0, 0])
    one = lambda d: gate([100.0], 100.0 - d)
    burst = lambda d: gate([100.0, 150.0, 200.0, 250.0, 300.0], 300.0 - d)
    ltd = {"-10": one(-10), "-25": one(-25), "-50": one(-50), "B-120": burst(-120), "B-200": burst(-200)}
    none = {"+10": one(10), "+25": one(25), "-100": one(-100), "-120": one(-120), "-200": one(-200)}
    assert min(ltd.values()) > 0.05 and max(none.values()) == 0.0, (ltd, none)
    # no veto -> the EPSP itself is a T drive: one extra event at the pre arrival
    assert _sig(_vrec([[301.0]] * 3), dict(P, w_V=0.0), "T").sum() == 3
    assert _sig(_vrec([[301.0]] * 3), P, "T").sum() == 0
    return "LTD gates " + " ".join(f"{k} {v:.2f}" for k, v in ltd.items()) + "; none: all 0"


def test_t_drive3_equals_t_drive2():
    """When every synapse sees each post AP as a local event (same sample) and the EPSP is vetoed, t_drive 3 gives
    exactly the t_drive 2 gate."""
    from batch_v2 import _sig
    posts = [100.0, 150.0, 200.0]; pre = 260.0
    r = _vrec([np.array(posts + [pre + 0.5])] * 3, postspikes=posts, prespikes=[pre])
    arr = MV.arrival_index(np.array([pre]), np.zeros(3), r["t"])
    Pb = dict(tau_E1=100.0, theta_Te=0.8, tau_T=40.0, theta_Tg=0.1)
    g2 = MV.feature_T(_sig(r, dict(Pb, t_drive=2), "T"), r["t"], arr, dict(Pb, t_drive=2))
    g3 = MV.feature_T(_sig(r, dict(Pb, t_drive=3), "T"), r["t"], arr, dict(Pb, t_drive=3))
    assert g2.max() > 0.05 and np.allclose(g2, g3, atol=1e-12), (g2, g3)
    return f"gate {g3[0, 0]:.3f} == t_drive 2"


def _crec(cev, prespikes=(300.0,), postspikes=(), n=3, dt=0.1, T=1000.0, delay=0.0):
    t = np.arange(0, T, dt)
    return dict(pair="p", proto="q", syn=np.arange(n), t=t, cev=np.array(cev, float), prespikes=np.array(prespikes),
                postspikes=np.array(postspikes), delay=np.full(n, delay))


def _brute_t4(x, dt, arr, P):
    """Direct sample-by-sample reference for t_drive 4: crossing detector with hysteresis on x = -ica_VDCC, own
    glutamate-bound state b (jump 1 at each own arrival, exp decay tau_d_NMDA, cap 1), S += 1 - b at each event,
    S decays with tau_E1, u = pos(S - theta_Te), T low-passed (model_v2 Euler/exact coefficients). Returns T at
    each arrival sample."""
    n = len(x); t = dt * np.arange(n); K = P["K_ca"] * P["K_mult"]
    ai = {int(np.searchsorted(t, a - 1e-9)) for a in arr}
    aE = MV._coef(dt, P["tau_E1"])[0]; ab = np.exp(-dt / P["tau_d_NMDA"])
    aT, bT = MV._coef(dt, P["tau_T"])
    S = b = T = 0.0; armed = True; out = {}
    for k in range(n):
        if k > 0:
            S *= aE; T = aT * T + bT * max(Sprev - P["theta_Te"], 0.0) / P["Te_scale"]; b *= ab
        if k in ai:
            out[k] = T; b += 1.0                         # T read at the arrival, before the own release
        if armed and x[k] > K:
            armed = False; S += 1.0 - min(b, 1.0)
        elif not armed and x[k] < K * MV.CA_HYST:
            armed = True
        Sprev = S
    return np.array([out[k] for k in sorted(out)])


def test_t_drive4_weight():
    """ca_events / glu_weight: one event per excursion (hysteresis at K/2), a sub-K bump is not an event, and the
    weight is 1 - b: 1 with no own release, exp(-dt/70) decayed after a release, 0 at the release itself, capped
    at 1 for stacked releases."""
    dt = 0.025; t = np.arange(0, 400, dt); K = MV.DEFAULTS["K_ca"]
    x = np.zeros(len(t))
    for t0, a in ((50.0, 3 * K), (50.5, 3 * K), (200.0, 0.9 * K), (300.0, 1.5 * K)):   # 50 and 50.5 merge (no re-arm)
        x += a * np.exp(-np.maximum(t - t0, 0) / 1.0) * (t >= t0)
    ev = MV.ca_events(x, dt)
    assert len(ev) == 2 and abs(ev[0] - 50.0) < 0.1 and abs(ev[1] - 300.0) < 0.1, ev
    assert len(MV.ca_events(0.5 * K * np.ones(10), dt)) == 0
    e, w = MV.glu_weight([10.0, 40.0, 70.0, 100.0, np.nan], [40.0], 70.0)
    assert np.allclose(w, [1.0, 0.0, 1 - np.exp(-30 / 70), 1 - np.exp(-60 / 70)]), w        # own release at 40
    e, w = MV.glu_weight([50.0], [10.0, 20.0, 30.0, 40.0], 70.0)
    assert w[0] == 0.0, w                                # four stacked releases: b capped at 1, weight 0, not negative
    return "events %s; weights %s" % (np.round(ev, 2).tolist(), np.round(w, 3).tolist())


def test_t_drive4_events_fast():
    """The block-skipping ca_events equals the plain sample-by-sample detector (random spiky traces, length not a
    multiple of the block, events at block edges, trace starting above K)."""
    K = MV.DEFAULTS["K_ca"]; rng = np.random.default_rng(1); ntot = 0
    def naive(x, K, h=0.5):
        armed = True; o = []
        for k, v in enumerate(x):
            if armed and v > K: armed = False; o.append(k)
            elif (not armed) and v < K * h: armed = True
        return np.array(o, float)
    for n in (50000, 20011, 4096 * 5 + 7):
        x = np.abs(rng.normal(0, 0.1 * K, n)) + 0.3 * K * (rng.random() < 0.5)
        for _ in range(30):
            k0 = rng.integers(0, n); x[k0:k0 + rng.integers(1, 60)] += rng.uniform(0.2, 6) * K
        x[:5] = 3 * K
        for m in (0.5, 1.0, 2.0):
            a = MV.ca_events(x, 1.0, K=K * m, block=256); b = naive(x, K * m)
            assert np.array_equal(a, b), (n, m, a[:5], b[:5])
            ntot += len(b)
    return f"{ntot} events identical"


def test_t_drive4_brute_force():
    """The event solution (batch_v2._sig -> impulses -> feature_T) equals the sample-by-sample brute force with the
    detector, b, S, T integrated directly, for K_ca, K_ca/2 (cev_lo) and 2 K_ca (cev_hi)."""
    from batch_v2 import _sig
    dt = 0.1; t = np.arange(0, 1500, dt); K = MV.DEFAULTS["K_ca"]; rng = np.random.default_rng(5)
    x = np.zeros((2, len(t)))
    for i in range(2):
        for t0 in np.sort(rng.uniform(20, 1300, 12)):
            x[i] += rng.uniform(0.3, 6) * K * np.exp(-np.maximum(t - t0, 0) / rng.uniform(1, 4)) * (t >= t0)
    pre = np.array([120.0, 150.0, 700.0, 900.0, 930.0]); delay = 0.7
    P = dict(t_drive=4, tau_E1=100.0, theta_Te=0.6, tau_T=40.0, theta_Tg=0.05)
    arr_idx = MV.arrival_index(pre, [delay] * 2, t)
    n_ev = []
    for m, key in ((1.0, "cev"), (0.5, "cev_lo"), (2.0, "cev_hi")):
        cev = [MV.ca_events(x[i], dt, K=K * m) for i in range(2)]
        nm = max(len(c) for c in cev); A = np.full((2, nm), np.nan)
        for i, c in enumerate(cev):
            A[i, :len(c)] = c
        r = _crec([[1.0]] * 2, prespikes=pre, n=2, dt=dt, T=1500.0, delay=delay)
        r["cev"] = r["cev_lo"] = r["cev_hi"] = A
        Pm = dict(P, K_mult=m)
        g = MV.feature_T(_sig(r, Pm, "T"), r["t"], arr_idx, Pm)
        Pf = {**MV.DEFAULTS, **Pm}
        for i in range(2):
            ref = _brute_t4(x[i], dt, pre + delay, Pf)
            # brute force works on the arrival sample; feature_T reads T there too
            Tb = np.tanh(np.maximum(ref - P["theta_Tg"], 0.0))
            assert np.allclose(g[i], Tb, atol=1e-12, rtol=1e-9), (m, i, g[i], Tb)
        n_ev.append(nm)
    assert g.max() > 0.01, "vacuous"
    return f"event solution == brute force at K x0.5/1/2 (max events {n_ev})"


def test_t_drive4_window():
    """t_drive 4 on synthetic Ca events reproduces the Sjostrom 2003 window: a bAP event 10-50 ms before the pre
    arrival (b = 0) drives T, a bAP after the arrival is down-weighted by 1 - b (pre-before-post +10 gives
    weight 1 - exp(-10/70) of the LTD unit... which is under theta_Tg), the own EPSP event at the arrival is
    weighted 0, and trains of 5 bAPs build T."""
    from batch_v2 import _sig
    P = dict(t_drive=4, tau_E1=100.0, theta_Te=0.8, tau_T=40.0, theta_Tg=0.5)
    def gate(posts, pre, epsp=True):
        ev = np.array(sorted([p + 1.0 for p in posts] + ([pre + 0.0] if epsp else [])))    # bAP +1 ms; EPSP Ca at arrival
        r = _crec([ev] * 3, prespikes=[pre]); sig = _sig(r, P, "T")
        g = MV.feature_T(sig, r["t"], MV.arrival_index(np.array([pre]), np.zeros(3), r["t"]), P)
        return float(g[0, 0]), sig
    one = lambda d: gate([100.0], 100.0 - d)[0]
    burst = lambda d: gate([100.0, 150.0, 200.0, 250.0, 300.0], 300.0 - d)[0]
    ltd = {"-10": one(-10), "-25": one(-25), "-50": one(-50), "B-120": burst(-120), "B-200": burst(-200)}
    none = {"+10": one(10), "+25": one(25), "-100": one(-100), "-120": one(-120), "-200": one(-200)}
    assert min(ltd.values()) > 0.05 and max(none.values()) == 0.0, (ltd, none)
    # own-glutamate weighting: an event 20 ms after the own arrival has weight 1 - exp(-20/70), at the arrival 0
    r = _crec([[320.0, 300.0]] * 3, prespikes=[300.0]); sig = _sig(r, P, "T")
    assert np.allclose(sig.sum(1), 1 - np.exp(-20 / 70)), sig.sum(1)
    assert np.allclose(_sig(_crec([[100.0]] * 3), P, "T").sum(1), 1.0)         # no own release before: full unit
    return "LTD gates " + " ".join(f"{k} {v:.2f}" for k, v in ltd.items()) + "; none: all 0"


def test_t_drive4_equals_t_drive2():
    """When every synapse sees each post AP as a Ca event (same sample) and its own EPSP event is at the arrival
    (weight 0), t_drive 4 gives exactly the t_drive 2 gate (events before any own release have weight 1)."""
    from batch_v2 import _sig
    posts = [100.0, 150.0, 200.0]; pre = 260.0
    r = _crec([np.array(posts + [pre])] * 3, postspikes=posts, prespikes=[pre])
    arr = MV.arrival_index(np.array([pre]), np.zeros(3), r["t"])
    Pb = dict(tau_E1=100.0, theta_Te=0.8, tau_T=40.0, theta_Tg=0.1)
    g2 = MV.feature_T(_sig(r, dict(Pb, t_drive=2), "T"), r["t"], arr, dict(Pb, t_drive=2))
    g4 = MV.feature_T(_sig(r, dict(Pb, t_drive=4), "T"), r["t"], arr, dict(Pb, t_drive=4))
    assert g2.max() > 0.05 and np.allclose(g2, g4, atol=1e-12), (g2, g4)
    return f"gate {g4[0, 0]:.3f} == t_drive 2"


def test_cacr_and_theta_n():
    """v2.4: spine Ca recovered from effcai is exact on a uniform grid; theta_N (threshold on N) lets a
    two-spike Ca burst drive NO while a single Ca transient does not; theta_N 0 keeps the old tanh(N)."""
    from batch_v2 import cacr_from_effcai
    sys.path.insert(0, os.path.join(ROOT, "analytical_method"))
    from extract import effcai_from_cai
    from constants import MIN_CA_CR
    dt = 0.25; t = np.arange(0, 1000, dt); rng = np.random.default_rng(0)
    cai = MIN_CA_CR + np.abs(rng.normal(0, 1e-3, (2, len(t))))
    eff = effcai_from_cai(cai, dt).astype(np.float32)
    u = cacr_from_effcai(eff, t)
    err = np.abs(u[:, :-1] - (cai[:, :-1] - MIN_CA_CR)).max()
    assert err < 2e-6, err                                   # float32 effcai: ~1e-6 of the 1e-3 signal
    pre = np.array([100.0]); arr = MV.arrival_index(pre, [0.0], t)
    def ca(n_ap):                                            # 2e-3 mM transients, 12 ms decay, 20 ms apart
        x = np.zeros(len(t))
        for k in range(n_ap):
            m = t >= 105 + 20 * k; x[m] += 2e-3 * np.exp(-(t[m] - 105 - 20 * k) / 12.0)
        return x[None]
    P = dict(no_drive=2, theta_NOc=5e-4, c_scale=1e-3, tau_NO=50.0, tau_Z=100.0, theta_Z=0.0, theta_N=12.0)
    k1, k2 = MV.feature_K(ca(1), t, arr, P).sum(), MV.feature_K(ca(2), t, arr, P).sum()
    assert k1 == 0.0 and k2 > 1.0, (k1, k2)
    a0 = MV.feature_K(ca(1), t, arr, dict(P, theta_N=0.0)); assert a0.sum() > 1.0
    old = MV.feature_K(ca(1) + 1e-4, t, arr, dict(theta_NO=1e-4, tau_Z=30.0))
    new = MV.feature_K(ca(1) + 1e-4, t, arr, dict(theta_NO=1e-4, tau_Z=30.0, theta_N=0.0))
    assert np.array_equal(old, new)
    return f"cacr recovery err {err:.1e} mM; K single {k1} vs double {k2:.1f} ms at theta_N 12"


def test_feature_cache_precompute():
    """Split T/K caches and the parallel precompute give exactly pre_features."""
    from batch_v2 import BatchV2
    B = BatchV2([MARKRAM], protocols=["10Hz_-10ms", "10Hz_10ms"], pairs={PAIR}, verbose=False)
    Ps = [dict(theta_T=th, tau_NO=tn, tau_Z=30.0) for th in (2e-4, 1e-3) for tn in (3.0, 10.0)]
    n = B.precompute(Ps, workers=2)
    assert n == 2 + 2, n                                   # 2 T keys + 2 K keys, not 4 x 2
    for P in Ps:
        for (tT, K), r in zip(B.features(P), B.recs):
            tT0, K0 = MV.pre_features(r["shaft_cai"], r["t"], r["arr"], P)
            assert np.array_equal(tT, tT0) and np.array_equal(K, K0)
    return f"{n} passes cover {len(Ps)} filter sets; identical to pre_features"


if __name__ == "__main__":
    pick = sys.argv[1] if len(sys.argv) > 1 else ""
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn) and pick in name:
            t0 = time.time()
            try:
                msg = fn()
                print(f"PASS {name} ({time.time() - t0:.1f}s): {msg}", flush=True)
            except Exception as e:
                fails += 1; print(f"FAIL {name}: {type(e).__name__}: {e}", flush=True)
    sys.exit(1 if fails else 0)
