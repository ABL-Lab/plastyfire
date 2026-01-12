import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Data setup
header = "pre_mtype,post_mtype,gsyn,gsynSD,nrrp,dtc,dtcSD,u,uSD,d,dSD,f,fSD,gsynSRSF,uHillCoefficient,synapse_type,nrrpSD,spinevol,spinevolSD,u_gsyn_r,gsyn_nrrp_r,spinevol_nrrp_r,spinevol_gsyn_r,gsynDist,nrrpDist,uDist,dDist,fDist,spinevolDist".split(',')

row1_data = "L5_TPC:B,L5_TPC:B,1.94,1.0,2.8,1.74,0.18,0.38,0.1,365.6,100.15,25.7,45.87,0.55,2.79,E2_L5TTPC,2.8,0.1628258728497028,0.08393086229366123,0.9,0.9,0.92,0.88,gamma,poisson,beta,gamma,gamma,gamma".split(',')
row2_data = "L5_TPC:B,L5_TPC:C,0.68,0.44,1.5,1.74,0.18,0.5,0.02,671.0,17.0,17.0,5.0,0.55,2.79,E2,1.5,0.05707298635968964,0.03692957940921094,0.9,0.9,0.92,0.88,gamma,poisson,beta,gamma,gamma,gamma".split(',')

# Helper to parse values
def parse_val(v):
    try:
        return float(v)
    except ValueError:
        return v

row1 = {k: parse_val(v) for k, v in zip(header, row1_data)}
row2 = {k: parse_val(v) for k, v in zip(header, row2_data)}

# Identify differences
diffs = {}
for k in header:
    v1 = row1[k]
    v2 = row2[k]
    if v1 != v2:
        diffs[k] = (v1, v2)

# Select numeric differences for plotting
numeric_diffs = {k: v for k, v in diffs.items() if isinstance(v[0], (int, float))}

if not numeric_diffs:
    print("No numeric differences to plot.")
    exit()

# Plotting
labels = list(numeric_diffs.keys())
v1_vals = [numeric_diffs[k][0] for k in labels]
v2_vals = [numeric_diffs[k][1] for k in labels]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(16, 6))
rects1 = ax.bar(x - width/2, v1_vals, width, label='Rest other pathways between L5_TPC:A and L5_TPC:B')
rects2 = ax.bar(x + width/2, v2_vals, width, label='L5_TPC:B,L5_TPC:A\nL5_TPC:B,L5_TPC:B\nL5_TPC:A,L5_TPC:A\nL5_TPC:A,L5_TPC:B')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Value', fontsize=18)
ax.set_title('Parameter differences within L5 synapses | SSCx circuit', fontsize=20)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=14)
ax.tick_params(axis='y', labelsize=18)
ax.legend(fontsize=14)

ax.bar_label(rects1, padding=3, fmt='%.2f', fontsize=12)
ax.bar_label(rects2, padding=3, fmt='%.2f', fontsize=12)

# Log scale if needed
max_val = max(max(v1_vals), max(v2_vals))
min_val = min(min(v1_vals), min(v2_vals))

if max_val / (min_val + 1e-9) > 100:
    ax.set_yscale('log')
    ax.set_ylabel('Value (log scale)', fontsize=18)

fig.tight_layout()

output_path = '/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/comparison_plot.png'
plt.savefig(output_path)
print(f"Plot saved to {output_path}")
