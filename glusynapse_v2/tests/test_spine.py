"""Spine-VDCC variant plumbing: --glusyn-globals parsing (pairrunner_edges_fit) and the variant file."""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)


def test_parse_globals():
    from plastyfire.pairrunner_edges_fit import _parse_glusyn_globals
    g = _parse_glusyn_globals('{"ljp_VDCC_GluSynapse": 10, "gca_bar_VDCC_GluSynapse": 0.2232}')
    assert g == {"ljp_VDCC_GluSynapse": 10.0, "gca_bar_VDCC_GluSynapse": 0.2232}
    try:
        _parse_glusyn_globals('{"ljp_VDCC": 10}')
    except ValueError:
        pass
    else:
        raise AssertionError("non-GluSynapse name accepted")


def test_sv_variant_file():
    v = json.load(open(os.path.join(ROOT, "glusynapse_v2/spine/delta_sv.json")))
    assert set(v["globals"]) == {"ljp_VDCC_GluSynapse", "gca_bar_VDCC_GluSynapse"}
    assert abs(v["globals"]["gca_bar_VDCC_GluSynapse"] / 0.0744 - 3.0) < 1e-9     # 3x the mod default


if __name__ == "__main__":
    test_parse_globals(); test_sv_variant_file(); print("PASS")
