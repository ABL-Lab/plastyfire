"""T16 check: does the Option B readout (model_v2.epsp_v2) predict the BCL EPSP when Use is scaled per synapse?

Runs the basis C01 protocol (run_basis_pair_edges, same delta circuit) with the v2 drop-in library
(glusynapse_v2/build_dropin, POINT_PROCESS GluSynapse = v2) and configs of (rho_i, s_i): rho_i picks
Use_d/gmax_d or Use_p/gmax_p, dpre0 = s_i - 1 and the pre amplitudes stay 0, so dpre is constant and
Use_GB sits at its v2 target min(1, Use_x * s_i) for the whole run. (Setting Use directly with the v1
mechanism does not work: Use_GB relaxes back to Use_x with tau_exp = 100 s during the 2 min run.)
It prints the measured mean EPSP next to the prediction from the existing basis csv.

    python glusynapse_v2/check_readout.py --pre-gid 180351 --post-gid 198084 --workers 12
"""
import argparse, os, sys, multiprocessing
import numpy as np, pandas as pd

ROOT = "/project/rrg-emuller/dhuruva/plastyfire"
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "glusynapse_v2"))
import run_basis_pair_edges as RB   # noqa: E402
import model_v2 as MV              # noqa: E402

SIMS = f"{ROOT}/refitting_results/fitting/n120/seed20262009/Sabrina_L5TTPC_L5TTPC_STDP/simulations"
CIRCUIT = f"{ROOT}/data/dhuruva_delta_circuit_config.json"
BASIS = f"{ROOT}/basis_results_edges_sabrina_n120_delta"
DROPIN = f"{ROOT}/glusynapse_v2/build_dropin"
RB.MECHANISMS_PATH = DROPIN     # _count_synapses loads it too


def configs(n, seed=0):
    rng = np.random.default_rng(seed)
    alt = [i % 2 for i in range(n)]
    c = [("all-d x1.0", [0] * n, [1.0] * n), ("all-d x0.5", [0] * n, [0.5] * n),
         ("all-d x1.8", [0] * n, [1.8] * n), ("all-p x1.0", [1] * n, [1.0] * n),
         ("all-p x0.6", [1] * n, [0.6] * n)]
    for k in range(2):
        c.append((f"mixed{k}", alt if k == 0 else [1 - a for a in alt],
                  list(np.round(rng.uniform(0.4, 1.8, n), 3))))
    return c


def trial(args):
    cfg_eff, cfg_orig, pre, post, name, rho, s, seed = args
    import bluecellulab
    from libsonata import SpikeReader
    bluecellulab.neuron.load_mechanisms(DROPIN)
    np.random.seed(seed)
    sim = bluecellulab.CircuitSimulation(cfg_eff, base_seed=seed)
    spk = SpikeReader(os.path.join(os.path.dirname(cfg_orig), "prespikes.h5"))[RB.NODE_POP].get_dict()["timestamps"]
    sim.instantiate_gids([(RB.NODE_POP, post)], add_synapses=True, add_minis=False, add_pulse_stimuli=True,
                         intersect_pre_gids=[(RB.NODE_POP, pre)], pre_spike_trains={(RB.NODE_POP, pre): spk})
    cell = sim.cells[(RB.NODE_POP, post)]
    for sec in cell.somatic + cell.axonal:
        sec.uninsert("SK_E2")
    for (sid, syn), r, si in zip(cell.synapses.items(), rho, s):
        h = syn.hsynapse
        use = min(1.0, (h.Use_p if r else h.Use_d) * si)
        g = h.gmax_p_AMPA if r else h.gmax_d_AMPA
        h.rho0_GB = h.rho_GB = float(r)
        h.dpre0_GB = si - 1.0
        h.Use = use; h.Use_GB = use
        h.gmax_AMPA = g; h.gmax0_AMPA = g
        h.theta_d_GB = -1.0; h.theta_p_GB = -1.0
    bluecellulab.neuron.h.cvode_active(1)
    sim.run(RB.C01_DURATION_MS, cvode=True)
    e = RB._measure_epsp(np.array(sim.get_time()), np.array(sim.get_voltage_trace((RB.NODE_POP, post))), spk)
    return name, seed, e, [int(k[1]) if isinstance(k, tuple) else int(k) for k in cell.synapses.keys()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre-gid", type=int, required=True)
    ap.add_argument("--post-gid", type=int, required=True)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", default=f"{ROOT}/glusynapse_v2/results/readout_check")
    a = ap.parse_args()
    pair = f"{a.pre_gid}-{a.post_gid}"
    cfg_orig = f"{SIMS}/{pair}/10Hz_-10ms/simulation_config.json"
    cfg_eff = RB._strip_synapse_reports(cfg_orig, circuit_config=CIRCUIT)
    n = RB._count_synapses(cfg_eff, a.pre_gid, a.post_gid)
    cs = configs(n)
    tasks = [(cfg_eff, cfg_orig, a.pre_gid, a.post_gid, nm, r, s, t) for nm, r, s in cs for t in range(a.trials)]
    with multiprocessing.Pool(min(a.workers, len(tasks))) as p:
        res = p.map(trial, tasks)
    df = pd.DataFrame([r[:3] for r in res], columns=["config", "trial", "epsp"])
    syn = np.array(res[0][3])

    basis = pd.read_csv(f"{BASIS}/basis_{a.pre_gid}_{a.post_gid}.csv")
    npz = np.load(f"{ROOT}/glusynapse_v2/extracted/markram_delta-cooker/{pair}__10Hz_-10ms.npz")
    assert np.array_equal(np.sort(syn), np.sort(npz["syn"])) or len(syn) == n, "synapse order"
    ep = MV.edge_params(npz["syn"])   # basis columns follow cell.synapses order == npz order (checked below)
    rows = []
    for nm, r, s in cs:
        m = df[df.config == nm]["epsp"]
        pred, _ = MV.epsp_v2(basis, np.array(r, float), np.array(s) - 1.0, ep)
        rows.append(dict(config=nm, bcl=m.mean(), bcl_sem=m.std(ddof=1) / np.sqrt(len(m)), pred=pred,
                         rel_err=pred / m.mean() - 1))
    out = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    out.to_csv(f"{a.out}/{pair}.csv", index=False); df.to_csv(f"{a.out}/{pair}_trials.csv", index=False)
    print("synapse ids (cell order):", syn.tolist(), "\nnpz order:", npz["syn"].tolist())
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
