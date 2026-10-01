"""Crop rendered figure pages (figs/src/*.png) to the panels being digitised (fractions of the page)."""
from PIL import Image
S = "/project/rrg-emuller/dhuruva/plastyfire/glusynapse_v2/cell_audit/ca_vs_data/figs/src/"
CROPS = {"ns06-06.png": [("ns06_fig5DEF", (0.50, 0.40, 1.0, 0.72))],
         "corn-02.png": [("corn_fig1E", (0.36, 0.50, 1.0, 0.72))],
         "corn-11.png": [("corn_fig7AB", (0.04, 0.06, 0.56, 0.40))]}
for f, cs in CROPS.items():
    im = Image.open(S + f); W, H = im.size
    for name, (a, b, c, d) in cs:
        im.crop((int(a * W), int(b * H), int(c * W), int(d * H))).save(S + name + ".png")
        print(name, im.size)
