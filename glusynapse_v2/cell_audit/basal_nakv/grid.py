"""Basal Na / Kv(SKv3_1) / Ka grid on the fitted L5->L5 protocols (diag_burst.py with pre trains).

    python glusynapse_v2/cell_audit/basal_nakv/grid.py --pair P --variants n1_k1_a1,... --out X.json
Protocols (T0 = first post pulse; dt>0 = pre first):
  sj{0.1,10,20,50}@{+10,-10}  Sjostrom 2001 / Markram 1997: 5 pairings at f (0.1 Hz: one pairing)
  s03_m25, s03_m120           Sjostrom 2003: one post AP, pre 25 / 120 ms later
  s03_b120                    5 post APs at 20 Hz, pre 120 ms after the last
  fi                          0.5 nA x 600 ms somatic step, no pre (spike count vs og)
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import diag_burst as D

OG = dict(gNaTgbar_NaTg=0.003, gSKv3_1bar_SKv3_1=0.0018178286052065377, gbar_Ka_kampa=0.002)   # og-delta basal
KA_ANTIC = 0.025
VAR = {f"n{n}_k{k}_a{a}": dict(gNaTgbar_NaTg=OG["gNaTgbar_NaTg"] * n, gSKv3_1bar_SKv3_1=OG["gSKv3_1bar_SKv3_1"] * k,
                               gbar_Ka_kampa=OG["gbar_Ka_kampa"] if a == "1" else KA_ANTIC)
       for n in (1, 3, 5) for k in (1, 3) for a in ("1", "A")}
PROTOS = []
for f in ("0.1", "10", "20", "50"):
    n = 1 if f == "0.1" else 5
    post = [k * 1000.0 / float(f) for k in range(n)]
    D.PULSE_CFG[f"sj{f}"] = f"sjostrom_{f}hz_dt+10ms"
    for dt in (10, -10):
        p = f"sj{f}@{dt:+d}"
        D.PROTOS[p] = post; D.PRE[p] = [x - dt for x in post]; D.TEND[p] = post[-1] + 150.0; PROTOS.append(p)
for p, post, pre in (("s03_m25", [0.0], [25.0]), ("s03_m120", [0.0], [120.0]),
                     ("s03_b120", [0.0, 50.0, 100.0, 150.0, 200.0], [320.0])):
    D.PULSE_CFG[p] = "sjostrom_burst5x20hz_dt-120ms" if p == "s03_b120" else "sjostrom_0.1hz_dt-25ms"
    D.PROTOS[p] = post; D.PRE[p] = pre; D.TEND[p] = pre[0] + 150.0; PROTOS.append(p)
D.PULSE_CFG["fi"] = "sjostrom_0.1hz_dt+10ms"; D.PROTOS["fi"] = [0.0]; D.PULSE_OVR["fi"] = (0.5, 600.0); D.TEND["fi"] = 650.0
PROTOS.append("fi")

if __name__ == "__main__":
    vf = os.path.join(HERE, "variants.json")
    json.dump(VAR, open(vf, "w"), indent=0)
    if "--protos" not in sys.argv:
        sys.argv += ["--protos", ",".join(PROTOS)]
    sys.argv += ["--var-cands", vf]
    D.main()
