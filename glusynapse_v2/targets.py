"""Experimental EPSP-ratio targets keyed by simwriter protocol id (tsk T10).

    TARGETS[(protocol_id, condition)] = (mean, sem, n, source)

condition names an in-silico manipulation that the offline model applies without new simulations:
  control       nothing blocked
  mglu_block    MCPG / AM251 / U73122: the mGluR -> PLC -> eCB -> CB1 chain, so A_mglu = 0
  post_nmdar    intracellular MK-801: postsynaptic NMDAR plasticity blocked, so post rho frozen
  nmdar_block   bath APV / pooled NMDAR antagonists: pre and post NMDARs blocked, so rho frozen, A_mglu = 0
                and A_NO = 0 (NO synthesis needs NMDAR Ca). The model then predicts exactly 1.0 (T25).
  no_block      L-NAME (NOS) / cPTIO (NO scavenger): presynaptic NO-LTP blocked, so A_NO = 0, rho and eCB kept
At L5-L5 synapses the T pathway is eCB -> CB1 with presynaptic NR2B NMDARs, not mGluR (Sjostrom 2003:
LY341495 does not block it). mglu_block is the historic name of the T-pathway block, so AM251 (CB1) and
ifenprodil (pre NR2B) map to it too.
Manipulations that change the calcium itself (D-APV bath, nimodipine, Ni2+) need new simulations.
They are listed in PENDING_SIM, not scored.

Nevian & Sakmann 2006 (L2/3, used for L5 as in Ebner 2019) values come from nevian/nevian_protocols.py,
checked against the paper text (Results, Figs 6, 7, 9, 10). The paper's dt is to the CLOSEST AP;
simwriter ids use the FIRST AP (configs/Ebner2019_L5TTPC_L5TTPC.yaml), so post-before-pre bursts shift
by (n_ap - 1) * ISI.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analytical_method"))
sys.path.insert(0, os.path.join(HERE, "..", "nevian"))

# Markram et al. 1997, 10 Hz (analytical_method/constants.py INVITRO)
MARKRAM = {
    "10Hz_-10ms": (0.7922, 0.0259, None, "Markram 1997"),
    "10Hz_5ms":   (1.2038, 0.0644, None, "Markram 1997"),
    "10Hz_10ms":  (1.2013, 0.0626, None, "Markram 1997"),
}
# Zilberter et al. 2009 Cereb Cortex p.2315: paired L5-L5, 5 + 5 at 10 Hz, -10 ms, AM251 2 uM (CB1 block -> mglu_block).
# Drug arm of the Markram 10 Hz -10 protocol (Zilberter repeats/intervals differ slightly; simulated on the Markram workdirs).
ZILBERTER_L5 = {
    ("10Hz_-10ms", "mglu_block"): (1.07, 0.08, 3, "Zilberter 2009 p.2315: L5-L5 5+5 10 Hz -10 ms, AM251 2 uM"),
}

# fit.py --anchor-tails: no plasticity far from coincidence (not data; same anchor as the cooker fit)
MARKRAM_TAILS = {f"10Hz_{d}ms": (1.0, 0.15, None, "tail anchor (fit.py TAIL_*)") for d in (-50, -30, 30, 50)}

# nevian_protocols.py id -> simwriter id (first-AP dt)
NEVIAN_ID_MAP = {
    "pre_only_far":                    "nevian_3ap_50hz_dt-90ms",
    "LTD_3ap_50hz_dt-10ms":            "nevian_3ap_50hz_dt-50ms",
    "neutral_3ap_50hz_dt-30ms":        "nevian_3ap_50hz_dt-30ms",
    "mild_LTP_3ap_50hz_dt-10ms_burst": "nevian_3ap_50hz_dt-10ms",
    "LTP_3ap_50hz_dt+10ms":            "nevian_3ap_50hz_dt+10ms",
    "post_only_far":                   "nevian_3ap_50hz_dt+50ms",
    "epsp_only":                       "nevian_epsp_only",
    "LTD_1ap_50hz_dt-10ms":            "nevian_1ap_dt-10ms",
    "nochange_1ap_50hz_dt+10ms":       "nevian_1ap_dt+10ms",
    "LTD_2ap_50hz_dt-10ms":            "nevian_2ap_50hz_dt-30ms",
    "LTP_2ap_50hz_dt+10ms":            "nevian_2ap_50hz_dt+10ms",
    "nochange_3ap_20hz_dt+10ms":       "nevian_3ap_20hz_dt+10ms",
    "LTD_3ap_20hz_dt-10ms":            "nevian_3ap_20hz_dt-110ms",
    "LTP_3ap_100hz_dt+10ms":           "nevian_3ap_100hz_dt+10ms",
    "LTD_3ap_100hz_dt-10ms":           "nevian_3ap_100hz_dt-30ms",
    # pharmacology that the offline model can apply
    "LTD_3ap_50hz_dt-10ms_mcpg":       ("nevian_3ap_50hz_dt-50ms", "mglu_block"),
    "LTP_3ap_50hz_dt+10ms_mcpg":       ("nevian_3ap_50hz_dt+10ms", "mglu_block"),
    "LTP_from_LTD_3ap_100hz_dt-10ms_mcpg": ("nevian_3ap_100hz_dt-30ms", "mglu_block"),
    "LTD_3ap_50hz_dt-10ms_u73122":     ("nevian_3ap_50hz_dt-50ms", "mglu_block"),
    "LTD_3ap_50hz_dt-10ms_am251":      ("nevian_3ap_50hz_dt-50ms", "mglu_block"),
    "LTD_3ap_50hz_dt-10ms_mk801":      ("nevian_3ap_50hz_dt-50ms", "post_nmdar"),
    "LTP_3ap_50hz_dt+10ms_mk801":      ("nevian_3ap_50hz_dt+10ms", "post_nmdar"),
}
# T25: paired L5-L5 pharmacology, Sjostrom, Turrigiano & Nelson 2003 Neuron 39:641 (papers/). All unitary
# L5-L5 pairs, P12-P21 visual cortex, 0.1 Hz baseline. Protocols as in Sjostrom 2001 (the yaml sjostrom_* ids).
SJOSTROM_DRUGS = {
    ("sjostrom_0.1hz_dt-10ms", "mglu_block"):  (1.01, 0.04, 4, "Sjostrom 2003 text p.644 + Fig 8B: AM251 0.1 Hz tLTD"),
    ("sjostrom_20hz_dt-10ms", "mglu_block"):   (1.02, 0.07, 7, "Sjostrom 2003 text p.644 + Fig 8B: AM251 20 Hz tLTD"),
    ("sjostrom_0.1hz_dt-10ms", "nmdar_block"): (1.07, 0.04, 7, "Sjostrom 2003 text p.644 + Fig 8B: NMDAR block "
                                                "(APV/MK801/ifenprodil pooled; ifenprodil alone 1.06 +- 0.06 n=5)"),
    ("sjostrom_20hz_dt-10ms", "nmdar_block"):  (1.02, 0.07, 8, "Sjostrom 2003 Fig 8B squares at 20 Hz [digitised]"),
    ("sjostrom_50hz_dt+10ms", "mglu_block"):   (1.59, 0.19, 7, "Sjostrom 2003 Fig 7D: ifenprodil does not block "
                                                "50 Hz tLTP (ctrl 1.64 +- 0.16 n=16) [digitised]"),
    ("sjostrom_50hz_dt+10ms", "nmdar_block"):  (1.04, 0.05, 4, "Sjostrom 2003 Fig 7D: APV blocks 50 Hz tLTP [digitised]"),
}
# T25 step 6: Sjostrom 2003 Fig 9C (new control experiments, 0.1 Hz induction), protocols in
# configs/Sjostrom2003_L5TTPC_L5TTPC.yaml. Fig 9C pools dt -120/-200/-400 ms into one bar; we simulate -120 and
# -200, so each gets the pooled value with SEM * sqrt(2), i.e. the bar counts once. Burst: 5 post APs at 20 Hz,
# one pre spike 120 or 200 ms after the last AP. [digitised, +-0.01]
_R2 = float(np.sqrt(2.0))
SJOSTROM_2003_TIMING = {
    ("sjostrom_0.1hz_dt-25ms", "control"):         (0.65, 0.09, 6, "Sjostrom 2003 Fig 9C Ctrl dt -25 ms [digitised]"),
    ("sjostrom_0.1hz_dt-120ms", "control"):        (1.05, 0.07 * _R2, 11, "Sjostrom 2003 Fig 9C Ctrl dt -120/-200/-400 pooled "
                                                    "(1.05 +- 0.07) [digitised]; SEM x sqrt2, shared with -200"),
    ("sjostrom_0.1hz_dt-200ms", "control"):        (1.05, 0.07 * _R2, 11, "as -120 (pooled bar)"),
    ("sjostrom_burst5x20hz_r50_dt-120ms", "control"):  (0.79, 0.06 * _R2, 7, "Sjostrom 2003 Fig 9C Burst, pre 120 or 200 ms after "
                                                    "the burst (0.79 +- 0.06) [digitised]; SEM x sqrt2, shared with -200"),
    ("sjostrom_burst5x20hz_r50_dt-200ms", "control"):  (0.79, 0.06 * _R2, 7, "as -120 (pooled bar)"),
}
def _combine(rows):
    """Several experiments on one (protocol, condition), e.g. MCPG, AM251 and U73122 all block the
    mGluR -> eCB chain: inverse-variance mean and SEM, n summed."""
    if len(rows) == 1:
        return rows[0]
    w = np.array([1.0 / r[1] ** 2 for r in rows])
    m = float(np.sum(w * [r[0] for r in rows]) / w.sum())
    return (m, float(1.0 / np.sqrt(w.sum())), sum(r[2] or 0 for r in rows) or None,
            " + ".join(r[3] for r in rows))


# T30: Sjostrom, Turrigiano & Nelson 2007 Neuropharmacology 52:176 (papers/), unitary L5 TTPC -> L5 TTPC, visual
# cortex P14-21. Induction: paired 200 ms depolarising steps in pre and post, 30x at 0.1 Hz (APs ~40-60 Hz), protocols
# in configs/Sjostrom2007_L5TTPC_L5TTPC.yaml. Scored on the per-cell amp (>= 10 APs in 200 ms); the fixed 1.2 nA id
# sjostrom07_step200ms_1.2nA_pair is a sensitivity check, not a target (the same cells would count twice).
# NO block: L-NAME 1.30 +- 0.12 n=8 and cPTIO 1.40 +- 0.093 n=6 (Results 3.3), inverse-variance combined (the text's
# pooled "136%, n=14" has no SEM). Pre-only / post-only: Fig 1B open circles, "Control, n=12" pooled, ~1.00 over
# 15-60 min; the caption says the error bars are too small to show, so SEM 0.05 is ASSUMED; SEM x sqrt2 per id as
# the pooled bar counts once. AM251 wash-in alone (0.98 +- 0.016 n=4) is a drug-only control: nothing to score.
SJOSTROM_2007 = {
    ("sjostrom07_step200ms_pair", "control"):     (1.62, 0.07, 31, "Sjostrom 2007 Results 3 + Fig 1B: paired 200 ms "
                                                   "steps, LTP 162 +- 7% [text]"),
    ("sjostrom07_step200ms_pair", "mglu_block"):  (2.13, 0.22, 10, "Sjostrom 2007 Results 3.3 + Fig 5A: AM251 900 nM, "
                                                   "213 +- 22% [text]"),
    ("sjostrom07_step200ms_pair", "no_block"):    _combine([(1.30, 0.12, 8, "Sjostrom 2007 Results 3.3: L-NAME 100 uM "
                                                              "[text]"),
                                                             (1.40, 0.093, 6, "Sjostrom 2007 Results 3.3: cPTIO 50 uM "
                                                              "[text]")]),
    ("sjostrom07_step200ms_pre_only", "control"): (1.00, 0.05 * _R2, 12, "Sjostrom 2007 Fig 1B open circles, pre/post-only "
                                                   "pooled n=12 [digitised; SEM assumed]; SEM x sqrt2, shared with post_only"),
    ("sjostrom07_step200ms_post_only", "control"): (1.00, 0.05 * _R2, 12, "as pre_only (pooled bar)"),
}
# T31: Sjostrom, Turrigiano & Nelson 2004 J Neurophysiol 92:3338 (papers/), dLTD on unitary L5 -> L5: 1 pre spike (n=14; 3 spikes
# at 20-30 Hz n=3 pooled in the paper's n=17) + 250 ms SUBTHRESHOLD post step (threshold less 20 pA, peak Vm -52 mV), 50-60 pairings
# at 0.1-0.2 Hz; yaml configs/Sjostrom2004_L5TTPC_L5TTPC.yaml (id sjostrom04_dltd_step250ms; pre spike at step onset ASSUMED).
# Control: n-weighted mean of 0.70 +- 0.04 (n=14, single spikes) and 0.65 +- 0.14 (n=3, bursts) = 0.69 [text, Results]; SEM 0.04 =
# pooled SD (~0.16) / sqrt(17). Drugs (Fig 3B, text: ifenprodil "completely abolished" dLTD, P=0.16 vs pre-before-post control;
# AM251 P=0.81): values ~1.0 (AM251, n=4) and ~1.1 (ifenprodil, n=6) are read from Fig 3B; the paper gives no SEM, so SEM = control
# SD 0.16 / sqrt(n) is ASSUMED (0.08, 0.065). Both arms block the same eCB -> CB1 (pre NR2B) chain, so both map to the existing
# mglu_block condition (A_mglu = 0, as AM251 in SJOSTROM_DRUGS and Sjostrom 2007) and are inverse-variance combined as no_block is.
# NOT scored: pre-before-post control ~1.04 (n=7) and LFS (no pairing): no ids simulated yet.
SJOSTROM_2004 = {
    ("sjostrom04_dltd_step250ms", "control"):    (0.69, 0.04, 17, "Sjostrom 2004 Results: dLTD 70 +- 4% (n=14) and 65 +- 14% "
                                                  "(n=3), n-weighted [text; pooled SEM derived]"),
    ("sjostrom04_dltd_step250ms", "mglu_block"): _combine([(1.0, 0.08, 4, "Sjostrom 2004 Fig 3B AM251 ~1.0 [read from figure; "
                                                            "SEM assumed]"),
                                                           (1.1, 0.065, 6, "Sjostrom 2004 Fig 3B ifenprodil ~1.1 [read from "
                                                            "figure; SEM assumed]")]),
}
# Paired L5-L5 data not scored (T25). Not targets.
#   (description, value, sem, n, source)
PENDING_PROTOCOL = [
    ("Fig 9B control curve, dt -200/-120/-100/-50/-25/-10/+10/+25/+50 ms",
     "0.98/1.08/0.99/0.65/0.67/0.71/1.34/0.92/0.98", "0.08/0.09/0.06/0.06/0.06/0.08/0.07/0.06/0.10", None,
     "Sjostrom 2003 Fig 9B [digitised]; filled points 'reproduced from Sjostrom 2001', legend says induction at 0.1 Hz, "
     "but +10 ms is 1.34 vs 0.97 for 0.1 Hz in 2001 Fig 1D, so the frequency is unclear until the 2001 PDF is read"),
    ("pre-only / post-only 30 Hz trains", "~1.0", None, None, "Sjostrom 2003 Fig 1B (not 0.1 Hz; the yaml 0.1 Hz "
     "pre_only/post_only runs are model sanity checks only)"),
]
# L2/3 pathways (L23_PATHWAYS.md). Unitary paired recordings only, on protocols that already have workdirs; no
# extracellular-stimulation data. paired_l23l5 = the csv letzkus_* and sjostrom_50hz_dt+10ms@distal rows (L2/3->L5 TTPC
# pairs; '@pooled' -> no suffix) plus the rows below. Location suffixes (@proximal/@distal) need a per-pair loc
# selection (eval_l23l5.py): BatchV2.objective without a loc map would score @proximal on all pairs.
PAIRED_L23L5_EXTRA = {
    ("sjostrom_50hz_dt+10ms", "control"): (1.06, 0.09, 19, "Sjostrom & Hausser 2006 Results (Fig 1 text): all unitary "
                                           "L2/3->L5 pairs, 50 Hz +10 ms (106 +- 9%, n=19) [text]"),
    ("letzkus_3ap_200hz_dt+10ms@distal", "nmdar_block"): (1.03, 0.07, None, "Letzkus 2006 Fig 7A text: 50 uM APV "
                                                          "during induction, distal (n in figure only) [text]"),
    ("letzkus_3ap_200hz_dt-10ms@distal", "nmdar_block"): (0.97, 0.06, None, "Letzkus 2006 Fig 7B text: APV, "
                                                          "distal (n in figure only) [text]"),
    # Letzkus 2006 null controls (Fig 2C, Methods; unpaired 1AP pooled n=11 not fitted, it is these two groups pooled).
    # Offset sign not stated in the paper: dt -500 ms = post leads (pre 500 ms after the first AP), an assumption.
    # Burst +500 ms scored on ~60 pairs only: 55 of 115 were dropped by the spike-count guardrail (bursts losing an AP or over-firing), so the surviving set is biased towards cells whose bursts do not adapt.
    ("letzkus_nopost", "control"): (0.98, 0.06, 6, "Letzkus 2006 Methods: no postsynaptic activity [text]"),
    ("letzkus_3ap_200hz_dt-500ms", "control"): (0.99, 0.05, 5, "Letzkus 2006 Methods: AP burst offset by 500 ms; sign assumed -500 [text]"),
}
# L2/3 PC -> L2/3 PC, unitary pairs: Zilberter et al. 2009 Cereb Cortex 19:2308 (rat visual cortex P14-21, 32-34 C, 2 mM Ca /
# 1 mM Mg, 40 pairings every 5 s, mean +- SEM). Ids: configs/Zilberter2009_L23PC_L23PC.yaml. All values printed in the text
# (none read from figures). Full table, unusable rows and why: rho_redesign/L23L23_TARGETS.md. Hardingham 2007 gives no
# fit target (10 mM EGTA post pipette, 23 C, LTP/nc/LTD split instead of a mean): validation only.
# AM251 -> mglu_block (the CB1 / T-pathway block, as everywhere). Not scored: APV (nmdar_block makes the model exactly 1.0,
# data 0.70 / 0.73), CPCCOEt + EGLU (postsynaptic mGluR, no model condition), D890, BAPTA, Mg-free (new sims), ext. EPSP.
PAIRED_L23L23 = {
    ("zilberter_1ap_dt+10ms", "control"): (0.64, 0.07, 6, "Zilberter 2009 Results + Fig 2B, 3: single pre/post AP +10 ms [text]"),
    ("zilberter_1ap_dt-10ms", "control"): (0.56, 0.06, 4, "Zilberter 2009 Results: post-pre -10 ms, data not shown [text]"),
    ("zilberter_pre_only", "control"): (0.98, 0.08, 5, "Zilberter 2009 Results + Suppl Fig 1B: pre alone, same frequency [text]"),
    ("zilberter_5ap_10hz_dt+10ms", "control"): (0.76, 0.07, 19, "Zilberter 2009 Results + Fig 2G: 5 pre/5 post 10 Hz +10 [text]"),
    ("zilberter_5ap_20hz_dt+10ms", "control"): (1.07, 0.11, 6, "Zilberter 2009 Results + Fig 2H: 5 pre/5 post 20 Hz +10 [text]"),
    ("zilberter_5ap_20hz_dt-10ms", "control"): (0.93, 0.07, 5, "Zilberter 2009 Results: 5-5 post-pre 20 Hz -10, data not "
                                                 "shown [text]"),
    ("zilberter_train10_50hz_dt+4ms_last", "control"): (1.49, 0.12, 11, "Zilberter 2009 Results + Fig 4A,D: train-LTP, pre "
                                                         "3-5 ms before 10th AP [text]"),
    ("zilberter_train10_50hz_dt-4ms_last", "control"): (0.99, 0.09, 6, "Zilberter 2009 Results + Fig 4D: pre 3-5 ms after "
                                                        "10th AP [text]"),
    ("zilberter_train10_50hz_dt-10ms_last", "control"): (0.72, 0.05, 13, "Zilberter 2009 Results + Fig 4C,D: train-LTD, pre "
                                                          "5-12 ms after 10th AP (sim 10 ms, Fig 4C inset) [text]"),
    ("zilberter_train10_50hz_dt+5ms", "control"): (0.97, 0.06, 4, "Zilberter 2009 Results + Fig 4B: pre 5 ms before 1st AP "
                                                   "[text]"),
    ("zilberter_train10_50hz_post_only", "control"): (1.03, 0.04, 4, "Zilberter 2009 Results + Suppl Fig: train alone [text]"),
    ("zilberter_train4_50hz_dt+4ms_last", "control"): (0.76, 0.07, 4, "Zilberter 2009 Results + Fig 5D: 4 bAPs, dt 4 ms [text]"),
    ("zilberter_train8_50hz_dt+4ms_last", "control"): (1.15, 0.08, 8, "Zilberter 2009 Results + Fig 5D: 8 bAPs, dt 4 ms [text]"),
    ("zilberter_train10_50hz_dt+4ms_last", "mglu_block"): (1.73, 0.24, 4, "Zilberter 2009 Results + Fig 6C: AM251 2 uM, "
                                                            "train-LTP [text]"),
    ("zilberter_train10_50hz_dt-10ms_last", "mglu_block"): (0.73, 0.07, 7, "Zilberter 2009 Results + Fig 6C: AM251 2 uM, "
                                                             "train-LTD [text]"),
}
# L2/3 -> L2/3 extra target (2026-10-01; rho_redesign/L23L23_EXTRA.md): Egger, Feldmeyer, Sakmann 1999 Nat Neurosci 2:1098 (paired recordings,
# Chindemi egg99_01: L2/3 PC pair, rat barrel cortex P12-14, 5 pre + 5 post APs at 20 Hz, +10 ms, 10 sweeps every 10 s). Group
# paired_l23l23_egger, load together with paired_l23l23 (extraction l23l23extra_*). Banerjee 2014 (bnr14_01/02; Fig 3C)
# is a PAIRED recording (mouse, room temperature 22-24 C): kept VALIDATION ONLY (like Hardingham 23 C), in VALIDATION_L23L23.
PAIRED_L23L23_EGGER = {
    ("egger1999_5ap_20hz_dt+10ms", "control"): (1.299, 0.082, 12, "Egger 1999 (Chindemi egg99_01): 5 pre/5 post 20 Hz +10 ms, 10 sweeps at 0.1 Hz [csv]"),
}
VALIDATION_L23L23 = {
    # bnr14_02 (0.80 +- 0.06, n=6; MK-801 in the PRESYNAPTIC pipette only) is the same model protocol, a duplicate of this entry. dt -15 ms per Chindemi csv, unverified (text implies -10).
    ("banerjee2014_1ap_dt-15ms", "control"): (0.77, 0.07, 5, "Banerjee 2014 (bnr14_01): paired L2/3 pair, mouse RT, post-pre -15 ms (unverified), 100 sweeps at 0.2 Hz, EPSP slope [csv]; validation only"),
}
# L5 -> L5 extra paired targets (2026-10-01; rho_redesign/L5_EXTRA.md), Chindemi's biodata/paired_recordings.csv rows that the fits never
# simulated. Group paired_l5_extra (load together with paired_l5); extraction extracted/l5extra_delta-split1-prefire-vca.
# Markram 1997 (mrk97_01/02/04/05/06): 5 pre + 5 post, pre 5 ms before each post AP, 10 sweeps every 4 s (10 Hz = existing 10Hz_5ms).
# Sjostrom 2001 (sjs01_11..17): 0.1 Hz = single pair, 50 sweeps; 20-100 Hz = 5 + 5 APs, 15 sweeps every 10 s; dt 0 = pre spike with the post AP.
# n = csv sample_size (Sjostrom 2001 notes: between 3 and 7). 40/50 Hz dt 0 are the paper's pooled 40/50 Hz bar, split in the csv in two.
PAIRED_L5_EXTRA = {
    ("2Hz_5ms", "control"):  (0.9886, 0.04, 2, "Markram 1997 (mrk97_01) 2 Hz +5 ms [csv]"),
    ("5Hz_5ms", "control"):  (1.0161, 0.0727, 5, "Markram 1997 (mrk97_02) 5 Hz +5 ms [csv]"),
    ("20Hz_5ms", "control"): (1.3687, 0.0909, 11, "Markram 1997 (mrk97_04) 20 Hz +5 ms [csv]"),
    ("30Hz_5ms", "control"): (1.4239, 0.0745, 3, "Markram 1997 (mrk97_05) 30 Hz +5 ms [csv]"),
    ("40Hz_5ms", "control"): (1.5008, 0.12, 4, "Markram 1997 (mrk97_06) 40 Hz +5 ms [csv]"),
    ("sjostrom_0.1hz_dt0ms", "control"):   (0.910125, 0.139103, 5, "Sjostrom 2001 (sjs01_11) 0.1 Hz dt 0, 50 sweeps [csv]"),
    ("sjostrom_20hz_dt+25ms", "control"):  (0.6519, 0.18373, 5, "Sjostrom 2001 (sjs01_12) 20 Hz +25 ms [csv]"),
    ("sjostrom_20hz_dt0ms", "control"):    (0.760771, 0.099366, 5, "Sjostrom 2001 (sjs01_13) 20 Hz dt 0 [csv]"),
    ("sjostrom_20hz_dt-25ms", "control"):  (0.70462, 0.06994, 5, "Sjostrom 2001 (sjs01_14) 20 Hz -25 ms [csv]"),
    ("sjostrom_40hz_dt0ms", "control"):    (0.928899, 0.0384392, 5, "Sjostrom 2001 (sjs01_15) 40 Hz dt 0, pooled 40/50 Hz [csv]"),
    # sjs01_16 50 Hz dt 0 (0.92276 +- 0.03549) DROPPED: duplicate of the Fig 7C 40/50 Hz point = Fig 7D 40 Hz (sjs01_15), A10 SJ01_DT0_CHECK.md.
    # sjs01_17 100 Hz dt 0 (1.24456 +- 0.063635, n=5) DROPPED: only 10 of 24 pairs can fire 5 APs at 100 Hz (3 ms pulses), a biased subset.
}
# Sjostrom & Hausser 2006 L5->L5 50 Hz +10 ms (csv sjh06_01: 1.40 +- 0.06, n=34) is the stimulus of sjostrom_50hz_dt+10ms (paired_l5: Sjostrom 2001
# 1.57 +- 0.26, n=6). Group paired_l5_sjh06 REPLACES that paired_l5 entry by the inverse-variance combination (load after paired_l5).
PAIRED_L5_SJH06 = {
    ("sjostrom_50hz_dt+10ms", "control"): _combine([(1.57, 0.26, 6, "Sjostrom 2001 Fig 1D 50 Hz"),
                                                    (1.40, 0.06, 34, "Sjostrom & Hausser 2006 (csv sjh06_01)")]),
}
EBNER_CSV = os.path.join(HERE, "..", "ebner", "ebner_targets.csv")
PENDING_SIM = ["dapv", "nimodipine", "ni2", "nimo_plus_ni2"]


def _nevian():
    from nevian_protocols import __dict__ as nv
    out, seen = {}, set()
    for lst in [v for v in nv.values() if isinstance(v, list)]:
        for p in lst:
            if not hasattr(p, "protocol_id") or p.protocol_id not in NEVIAN_ID_MAP or p.protocol_id in seen:
                continue
            seen.add(p.protocol_id)
            m = NEVIAN_ID_MAP[p.protocol_id]
            sim_id, cond = (m, "control") if isinstance(m, str) else m
            out.setdefault((sim_id, cond), []).append(
                (p.expected_epsp_ratio, p.expected_sem, p.expected_n, f"Nevian 2006 {p.figure_ref} ({p.protocol_id})"))
    return {k: _combine(v) for k, v in out.items()}


def _ebner_csv():
    import csv
    out = {}
    for r in csv.DictReader(open(EBNER_CSV)):
        if not r.get("mean_ratio") or not r.get("sem"):
            continue      # no SEM: a model value (Ebner fit), not data
        loc = (r.get("synapse_location") or "all").split()[0]
        pid = r["protocol_id"] if loc == "all" else f'{r["protocol_id"]}@{loc}'
        out[(pid, "control")] = (float(r["mean_ratio"]), float(r["sem"]),
                                 int(float(r["n"])) if r.get("n") else None,
                                 f'{r["source_paper"]} {r.get("figure_or_table", "")} [{r.get("provenance", "")}]')
    return out


def load_targets(groups=("markram", "nevian", "ebner")):
    """groups: markram (data + tail anchors), nevian, ebner (csv), sjostrom_drugs, and paired_l5 (T25): unitary
    L5-L5 recordings only, i.e. Markram data without tail anchors, the csv sjostrom_* rows, SJOSTROM_DRUGS and
    SJOSTROM_2003_TIMING. sjostrom07: the Sjostrom 2007 step pairing (T30). paired_l23l5 / paired_l23l23: unitary L2/3->L5 TTPC / L2/3->L2/3 (L23_PATHWAYS.md);
    pathway groups cannot be combined."""
    t = {}
    if len({"paired_l5", "paired_l23l5", "paired_l23l23"} & set(groups)) > 1:
        raise ValueError("pathway groups share protocol ids (e.g. sjostrom_50hz_dt+10ms is L5-L5 1.57 in paired_l5 "
                         "but L2/3-L5 1.06 in paired_l23l5): load and score one pathway at a time")
    if "paired_l23l5" in groups:
        if os.path.isfile(EBNER_CSV):
            for (pid, c), v in _ebner_csv().items():
                if pid.startswith("letzkus_") or pid == "sjostrom_50hz_dt+10ms@distal":
                    t[(pid.replace("@pooled", ""), c)] = v
        t.update(PAIRED_L23L5_EXTRA)
    if "paired_l23l23" in groups:
        t.update(PAIRED_L23L23)
    if "paired_l23l23_egger" in groups:
        t.update(PAIRED_L23L23_EGGER)
    if "paired_l23l23_val" in groups:   # Banerjee 2014, score only (mouse, room temperature): never a fit target
        t.update(VALIDATION_L23L23)
    if "paired_l5" in groups:
        t.update({(k, "control"): v for k, v in MARKRAM.items()})
        t.update(ZILBERTER_L5)
        t.update(SJOSTROM_DRUGS)
        t.update(SJOSTROM_2003_TIMING)
        if os.path.isfile(EBNER_CSV):
            t.update({k: v for k, v in _ebner_csv().items()
                      if k[0].startswith("sjostrom_") and "@" not in k[0]})
    if "paired_l5_extra" in groups:     # 2026-10-01, see PAIRED_L5_EXTRA
        t.update(PAIRED_L5_EXTRA)
    if "paired_l5_sjh06" in groups:     # replaces the paired_l5 sjostrom_50hz_dt+10ms entry by the combined value
        t.update(PAIRED_L5_SJH06)
    if "sjostrom07" in groups:          # T30, L5-L5 paired: combine with paired_l5 once the step prefire exists
        t.update(SJOSTROM_2007)
    if "sjostrom04" in groups:          # T31, L5-L5 dLTD (subthreshold step): needs the prefire + extraction first
        t.update(SJOSTROM_2004)
    if "sjostrom_drugs" in groups:
        t.update(SJOSTROM_DRUGS)
    if "markram" in groups:
        t.update({(k, "control"): v for k, v in MARKRAM.items()})
        t.update({(k, "control"): v for k, v in MARKRAM_TAILS.items()})
    if "nevian" in groups:
        t.update(_nevian())
    if "ebner" in groups and os.path.isfile(EBNER_CSV):
        # ebner T7 csv (owned by the ebner_plastyfire_setup session): all 28 yaml ids.
        # Control values here take precedence over nevian_protocols.py, so a disagreement
        # is reported rather than silently kept.
        for key, val in _ebner_csv().items():
            if key in t and abs(t[key][0] - val[0]) > 1e-6:
                print(f"WARNING target mismatch {key}: nevian_protocols {t[key][:2]} vs ebner csv {val[:2]}")
            t[key] = val
    return t


if __name__ == "__main__":
    for (pid, cond), (m, s, n, src) in sorted(load_targets(tuple(sys.argv[1:]) or ("markram", "nevian", "ebner")).items()):
        print(f"{pid:32s} {cond:11s} {m:5.2f} ± {s:4.2f}  n={n}  {src}")
