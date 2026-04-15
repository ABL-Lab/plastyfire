import sys
sys.path.append("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire")
from plot_stdp_curves import process_results
from submit_l5ttpc_traces import DHURUVA_PARAMS_V12

df = process_results("/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/DHURUVA_PARAMS_V12", DHURUVA_PARAMS_V12)
subset = df[df["pair"].isin(["187212-192870", "196748-180959"])]
print("Subset Data (the 2 pairs JAX used):")
print(subset)

print("\nMean EPSP ratio by dt for these 2 pairs:")
print(subset.groupby("dt")["epsp_ratio"].mean())
