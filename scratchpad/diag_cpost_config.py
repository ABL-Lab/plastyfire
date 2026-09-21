"""Isolate why the c_post single-AP search fails in _find_cpre_cpost but
succeeded in simwriter's check_electrical_constraint for the same gid.

Runs the identical search under three configs:
  A: the calibration config simwriter used   (no `conditions` block -> v_init -65)
  B: the production prefire config           (v_init -80)
  C: production config with v_init stripped  (isolates v_init as the variable)
"""
import json
import logging
import os
import sys
import tempfile

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                    format="%(levelname)s %(message)s")
logging.getLogger("plastyfire").setLevel(logging.DEBUG)

from plastyfire.simulator import spike_threshold_finder

GID = int(sys.argv[1]) if len(sys.argv) > 1 else 205995
PAIR = sys.argv[2] if len(sys.argv) > 2 else "13255-205995"

BASE = ("/lustre06/project/6077694/dhuruva/plastyfire/refitting_results/"
        "fitting/n100/seed19091997/L23PC_L5TTPC_STDP/simulations")
CALIB = os.path.join(BASE, "single_cells", "simulation_config.json")
PROD = os.path.join(BASE, PAIR, "10Hz_5ms", "prefire_simulation_config.json")

WIDTHS = [1.5, 3, 5]


def strip(cfg_path, drop_v_init):
    cfg = json.load(open(cfg_path))
    cfg.get("conditions", {}).get("mechanisms", {}).pop("GluSynapse", None)
    cfg.pop("reports", None)
    if drop_v_init:
        cfg.get("conditions", {}).pop("v_init", None)
    out_dir = tempfile.mkdtemp()
    cfg.setdefault("output", {})["output_dir"] = out_dir
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    json.dump(cfg, tmp)
    tmp.close()
    return tmp.name


def search(tag, cfg_path):
    print(f"\n{'=' * 70}\n{tag}\n  config: {cfg_path}")
    cfg = json.load(open(cfg_path))
    print(f"  v_init: {cfg.get('conditions', {}).get('v_init', 'ABSENT -> BCL default -65')}")
    for w in WIDTHS:
        res = spike_threshold_finder(cfg_path, GID, 1, 0.1, w, 1000.,
                                     0.05, 5., 100, "S1nonbarrel_neurons", True)
        if res is not None:
            print(f"  -> width {w} ms: FOUND amp={res['amp']:.4f} nA, "
                  f"{len(res['t_spikes'])} spike(s) at {res['t_spikes']}")
            return
        print(f"  -> width {w} ms: no single AP up to 5.0 nA")
    print("  -> FAILED at every width")


print(f"gid {GID}, pair {PAIR}")
search("A: calibration config (what qualified this cell)", CALIB)
search("B: production prefire config (what _find_cpre_cpost uses)",
       strip(PROD, drop_v_init=False))
search("C: production config, v_init removed", strip(PROD, drop_v_init=True))
