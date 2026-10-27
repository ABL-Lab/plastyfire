# New L5 emodel rebuild: ready to submit (prepared 2026-10-05; nothing submitted)
Use spine/submit_newemodel_s2.sh (+ worker spine/run_newemodel_s2.sh), NOT run_new_emodel_chain.sh. The chain lacks what s2g0321 did: synapse globals
(ljp25/gca0.0321, --glusyn-globals), K from sj07 post-only (kca, --k-ca), l5extra, Egger, cexp + Sabatini, T39 doublet audit; it also uses a different
sizing and og-delta pair indexes. The wrapper is s2g0321 verbatim with only the emodel swapped (circuit configs generated from the split2 ones by replacing
biophysical_neuron_models_dir; yamls derived from s2g0321/configs; pairs: subset24 csv = 24 L5 pairs, the same 120 L2/3->L5 pairs; audit keep lists give 22/102).
no_SK_E2 and v_init -80 are written by simwriter.py itself (so they are reproduced automatically); trace vars cai_CR,shaft_cai,ica_VDCC,v_seg are in run_newemodel_s2.sh pf().
## KEEP_SK=1 switch (default unchanged)
KEEP_SK=1 skips the no_SK_E2 circuit modification (simwriter.py) and every SK_E2 uninsert: simwriter amp search (fixhp), simulator_edges.py defaults (prefire pool
and c_pre/c_post), precompute_cpre_cpost.py, run_basis_pair_edges.py, cexp/measure_cexp.py. sbatch exports it to all jobs. py_compile check: job 22506582.
## What is skipped / reused
PATHWAYS default l5l5,l23l5. L2/3->L2/3 has L2/3 post cells (hocs copied unchanged, no new .mod): the s2g0321 extracts + basis_l23l23 + cexp stay valid; ~85 CPU-h saved.
EPSP bases are NOT reusable for L5 posts: run_basis_pair_edges simulates the EPSP on the post cell (basis5 29 + basis23_5 130 CPU-h).
## Dry run (hfcandX, delta-split2 pkg as stand-in)
36 jobs, expected 285 CPU-h (measured elapsed x CPUs; l5l5,l23l5), 370 with l23l23; sizes = s2g0321 sacct MaxRSS +25% / elapsed +50% (basis 12/15G 0:30, pfA 146G 0:45 (was
at its 117G limit), pfB 215G, pf07 81G, pfx 120G, pf23_5 265G (was at 212G limit) 1:00, pf23_5n 147G 1:15, ext_* 41-83G, kca 12G). Sims/cache/audit/cexp as measured.
## The one command (after ion_fitter gives NAME and DIR; mechanisms compiled into DEES_cell_packages/x86_64)
  cd /project/rrg-emuller/dhuruva/plastyfire && NE_NAME=<name> EMODEL_DIR=<dir> KEEP_SK=1 bash glusynapse_v2/spine/submit_newemodel_s2.sh
(drop KEEP_SK=1 if the note says no_SK_E2 stays; add DRYRUN=1 first; PATHWAYS=l5l5,l23l5,l23l23 for the full ~370 CPU-h). Outputs /scratch/dhuruva/<name>/: extracted/*_delta-<name>-ljp25g0321-*,
basis_*, cexp, doublets/*_keep_pairs.txt (needs audit_fin), k.txt. Check: seff on pfA/pfB/pf23_5 (at-limit risk), PASS line in sab_merge log, audit_fin done.
## DECISIONS.md line (add before submit; sed -i '3a - ...' glusynapse_v2/DECISIONS.md)
2026-10-05: NEW-EMODEL <name> rebuild (HEAVY, ~285 CPU-h l5l5+l23l5; 370 with l23l23): spine/submit_newemodel_s2.sh (copy of s2g0321 chain, KEEP_SK=<0|1>), tree /scratch/dhuruva/<name>, sizes from s2g0321 sacct 22230193-235 (MaxRSS +25%, elapsed +50%): basis5 12G/basis23_5 15G 0:30 (29+130 CPU-h; EPSP simulated on the post cell, not reusable), prefire pfA 146G pfB 215G pf07 81G pfx 120G pf23_5 265G pf23_5n 147G (113 CPU-h), ext 41-83G, cexp 2.2-5.5G, audit 6/23G; L2/3->L2/3 reused from s2g0321 (L2/3 hocs unchanged).

## After the rebuild (T59, 2026-10-05)
- LETZ_SPLIT=median in run_fit_v9.sh needs the `letzkus_distal_median` column in the new emodel's L23->L5 pair geometry csv. The rebuild (run_newemodel_s2.sh) does NOT call ebner/pair_geometry.py, so run it on the new emodel's L23->L5 pairs (small Slurm job) before the median-split refit. fit_v9 stops with an error if the column is missing. Rise times change with the emodel, so the median must be recomputed, not copied.

## Joint refit launch sequence (prepared 2026-10-05; all ports PASS on old emodel)
1. DECISIONS.md line (rebuild ~285 CPU-h), then `DRYRUN=1` and then the real submit_newemodel_s2.sh with KEEP_SK=1 (ion_fitter: SK kept at calibrated level).
2. After the rebuild: run ebner/pair_geometry.py on the new L23->L5 pairs (letzkus_distal_median), as a small Slurm job sized from its MEASURED line.
3. Refits (GPU 3g.40gb, fit 72G 3:30 + val 89G 0:15 afterok; basis 22507016 56.5 GB, ports 70.6 GB):
   `NE_NAME=<new> LETZ_SPLIT=median FAM=V3   bash glusynapse_v2/rho_redesign/run_fit_v9.sh submit A s5`
   plus one unseeded start (pen-offset auto). Seed from the best old-emodel V3 fit. (ECBX rejected by the user 2026-10-05: no Ca-source split.)
4. Judge V3 by joint chi2 (3 pathways, core hits) -> fig9 + data json -> reviewer -> companion tex push.
