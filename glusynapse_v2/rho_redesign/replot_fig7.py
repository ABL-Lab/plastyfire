"""Redraw fig7 from the saved scan CSVs (labels only; no scan rerun). Usage: replot_fig7.py SAVE_PREFIX FIG"""
import sys
import pandas as pd
import scan_vgate_amp as S0

save, fig = sys.argv[1], sys.argv[2]
M = pd.read_csv(save + ".csv"); S = pd.read_csv(save + "_syn.csv")
L5 = S.path == "L5"; L23 = S.path == "L23"; up = S.up.astype(bool)
grp = {"L5 LTP protos, up": S[L5 & S.proto.isin(S0.L5_LTP) & up],
       "L5 all, crossing": S[L5 & (S.vmax > 0)],
       "L2/3 1AP+10, up": S[L23 & (S.proto == "letzkus_1ap_dt+10ms") & up],
       "L2/3 3AP+10 distal, up": S[L23 & (S.proto == "letzkus_3ap_200hz_dt+10ms") & (S.letzkus_distal == True) & up],  # noqa: E712
       "L2/3 3AP+10 prox, up": S[L23 & (S.proto == "letzkus_3ap_200hz_dt+10ms") & (S.letzkus_distal == False) & up],  # noqa: E712
       "L2/3 SH50 distal, up": S[L23 & (S.proto == "sjostrom_50hz_dt+10ms") & (S.sh_distal == True) & up]}  # noqa: E712
print({k: len(v) for k, v in grp.items()})
S0.fig_(M, grp, fig)
print("wrote", fig)
