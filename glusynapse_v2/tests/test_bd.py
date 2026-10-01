"""delta-bd: runtime densities == the hoc distribute_distance expressions; og params reproduce og-delta exactly."""
import math, os, sys, re
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "emodel_bd"))
import bd


def hoc_eval(line, d):
    expr = re.search(r'"[^"]+", "([^"]+)"', line).group(1).replace("%.17g", repr(float(d)))
    return eval(expr, {"exp": math.exp})


def test_og_is_og_delta():
    g = bd.densities({}, np.array([0.0, 100.0, 400.0]))
    assert np.all(g["gNaTgbar_NaTg"] == 0.003) and np.all(g["gbar_Ka_kampa"] == 0.002)
    assert np.all(g["gCa_LVAstbar_Ca_LVAst"] == bd.OG_LVA)


def test_hoc_matches_runtime():
    for P in ({}, bd.ANTIC, dict(na=0.01, ka0=0.004, lam=80.0, m=3.0, dh=140.0)):
        lines = bd.hoc_lines(P)
        for d in (0.0, 37.5, 130.0, 260.0):
            g = bd.densities(P, d)
            for line, k in zip(lines, ("gNaTgbar_NaTg", "gbar_Ka_kampa", "gCa_LVAstbar_Ca_LVAst")):
                assert abs(hoc_eval(line, d) - float(g[k])) <= 1e-15 * max(1.0, abs(float(g[k]))), (P, d, k)


if __name__ == "__main__":
    test_og_is_og_delta(); test_hoc_matches_runtime(); print("PASS")
