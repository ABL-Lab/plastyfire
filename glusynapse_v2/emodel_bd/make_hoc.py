"""Write the delta-bd emodel package: SSCx-AAD-delta-emodels copied as is, cADpyr_L5TPC.hoc gets the delta-bd basal
distributions appended to biophys() (after the existing distribute_distance calls), plus a circuit config copy
pointing to it. The delta package and its circuit config are never modified.

    python glusynapse_v2/emodel_bd/make_hoc.py --params '{"na":..,"ka0":..,"lam":..,"m":..,"dh":..}' --name delta-bd
"""
import argparse, json, os, shutil, sys
sys.path.insert(0, os.path.dirname(__file__))
import bd

PKG = "/project/rrg-emuller/dhuruva/DEES_cell_packages"
SRC = os.path.join(PKG, "SSCx-AAD-delta-emodels")
CIRC = "/project/rrg-emuller/dhuruva/plastyfire/data/dhuruva_delta_circuit_config.json"
ANCHOR = '  distribute_distance(CellRef.somatic, "gIhbar_Ih"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True); ap.add_argument("--name", default="delta-bd")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    P = json.loads(a.params)
    dst = os.path.join(PKG, f"SSCx-AAD-{a.name}-emodels")
    if os.path.exists(dst) and not a.force:
        sys.exit(f"{dst} exists (use --force to overwrite this delta-bd package)")
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(SRC):
        if f.endswith(".hoc"):
            shutil.copy2(os.path.join(SRC, f), dst)
    hoc = os.path.join(dst, "cADpyr_L5TPC.hoc")
    lines = open(hoc).read().split("\n")
    i = next(i for i, l in enumerate(lines) if l.startswith(ANCHOR))
    add = ["  /* delta-bd: distance-dependent basal Na, Ka, LVA (glusynapse_v2/emodel_bd/bd.py) */"] + bd.hoc_lines(P)
    open(hoc, "w").write("\n".join(lines[:i + 1] + add + lines[i + 1:]))
    json.dump(dict(model=a.name, derived_from=SRC + "/cADpyr_L5TPC.hoc", params=P,
                   basal_formula=bd.__doc__.strip()), open(os.path.join(dst, f"{a.name.replace('-', '_')}_values.json"), "w"),
              indent=1, default=str)
    c = json.load(open(CIRC))
    for pop in c["networks"]["nodes"]:
        for k, v in pop["populations"].items():
            if v.get("biophysical_neuron_models_dir", "").endswith("SSCx-AAD-delta-emodels"):
                v["biophysical_neuron_models_dir"] = dst
    out = CIRC.replace("dhuruva_delta_circuit_config", f"dhuruva_{a.name}_circuit_config")
    json.dump(c, open(out, "w"), indent=4)
    print(dst, out)


if __name__ == "__main__":
    main()
