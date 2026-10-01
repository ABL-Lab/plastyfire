# Can the L5TTPC model make the Letzkus 2006 dendritic Ca difference?

**Setup.** `letzkus_ca.py` builds the real post cell (bluecellulab, delta circuit, same config as `diag_burst.py`) and applies somatic pulses only, with no pre spikes. The cells are 6 L23PC→L5TTPC pairs (Ebner2019 workdirs, read-only; proximal 2509-189325, 20168-191194, 20173-183306 and distal 6188-182372, 153238-206803, 10024-203348) plus 2 L5→L5 posts (181015-184976, 182339-200396). The script records at two kinds of site:
- the pairs' real synapses (n=68);
- GluSynapse clones placed along the soma→farthest-tuft path at 50–900 µm.

The Ca measure is the integral of GluSynapse `ica_VDCC` over −1..120 ms. Each L23 pair uses its own calibrated pulse (4.7 nA, 2 ms). The output is `letzkus_ca_*.json` and the plot is made by `letzkus_ca_plot.py`. The run was slurm job 22064266, MaxRSS 1.0 GB.

## Ca ratio vs distance (trunk path, median over 8 cells; og-delta / antic-delta)
| dist µm | bAP 1AP mV (og) | 3AP200/1AP | 3AP200/2AP50 |
|---|---|---|---|
| 50  | 83 | 1.95 / 2.24 | 0.98 / 1.18 |
| 150 | 44 | 2.70 / 4.72 | 1.27 / 2.06 |
| 250 | 21 | 3.11 / 5.07 | 1.72 / 2.26 |
| 450 | 11 | 2.48 / 2.85 | 1.79 / 1.96 |
| 650 | 7  | 1.59 / 1.78 | 1.34 / 1.36 |
| 900 | 5  | 1.21 / 1.27 | 1.13 / 1.12 |

The same measures at the real L23→L5 synapses (og-delta; antic-delta is within 10% except proximally):

| path distance | n synapses | bAP (mV) | 3AP200/1AP | 3AP200/2AP50 |
|---|---|---|---|---|
| <200 µm | 24 | 74 | 2.05 | 1.00 |
| 200–450 µm | 9 | 13 | 1.86 | 1.53 |
| >450 µm | 35 | 2 | 1.11 | 1.07 |

Spine VDCC charge at distances beyond 250 µm is about 1e-6, 300–1000× below the proximal value and at the level of the resting leak.

## bAP amplitude
The single bAP collapses to 44 mV at 150 µm, 11 mV at 450 µm and 7 mV at 650 µm. Letzkus gives no numbers in the text: Fig. 6A/D only shows the decremental bAP and Ca spikes at 660 µm. The literature they cite (Stuart 1997; Larkum 1999) reports bAPs of a few tens of mV at 450–650 µm. The model therefore attenuates several-fold too steeply.

## Ca spike / boosting
- **7 of 8 cells:** no Ca spike. The third AP of the 200 Hz burst reaches only 12 mV at 650 µm, and the 3AP/1AP Ca ratio *falls* with distance past about 300 µm. That is the opposite of Letzkus.
- **1 of 8 cells (10024-203348):** 3AP@200 Hz triggers an all-or-none, tree-wide Ca_HVA2 spike. The dendrite reaches about 88 mV at 550–900 µm and the soma fires a 5-spike burst. The 3AP/1AP ratio rises with distance: 12 at 158 µm, 2400 at 449 µm, about 6000 at 650 µm. The same cell shows no spike at 100 Hz or 50 Hz.
- **Protocol side-effect:** the 3AP@200 Hz pulse gives only 2 somatic APs in 6 of 8 cells.

## Verdict
**No, not robustly.** The cell provides (ii), the proximal burst increase: 3AP/1AP is about 2 at the real proximal synapses. It fails (iii) proximally, where 3AP200/2AP50 is about 1.0 (og-delta). For (i), the distal/proximal difference exists in 7 of 8 cells but has the wrong sign: distal synapses see essentially no bAP-driven Ca, and bursting boosts them less than proximal ones. As a result, no plasticity-parameter fit can make distal inputs flip sign with timing in those cells. Letzkus requires the opposite: a small single-AP Ca signal at distal synapses but a large burst-evoked Ca spike.

The apical dendrite has uniform NaTg (0.023), Ka_kampa (0.023) and Ca_HVA2 = Ca_LVAst (2.4e-3, suspiciously identical), with no distal Ca hot zone, on a thin trunk (about 1.5 µm at 450 µm). Two findings point at the cause:
- Zeroing apical Ka_kampa only raises the bAP at 440 µm from 5 to 14 mV. The failure is therefore mainly too little apical Na for the load, with Ka_kampa as a secondary cause.
- A distal Ca_LVA/HVA hot zone is missing. When a Ca spike does occur, Ca_HVA2 carries it.

antic-delta edits only basal Na/Ka, so it cannot fix this.
