# Plasticity protocols and targets: inventory for the joint fit (2026-09-30)

Goal: one GluSynapse_v2 parameter set fitted jointly on every paired pathway. Same parameters on every pathway;
pathway differences come only from each synapse's own traces (c_pre, c_post, spine/shaft Ca, VDCC, U_SE) and basis.
Fit targets: unitary paired recordings only. Extracellular data (Nevian 2006) is validation only.
Sources: targets.py, L23_PATHWAYS.md, TARGETS.md, extracted/, basis_results_edges_*.

## 1. Inventory

"VDCC traces" = extracted dir with the cev/vdcc fields the v4 kernel needs (-vca for L5, -vseg-rs for L2/3->L5).
"Joint" = in the current fit (run_fit_v4.sh JOINT=1, e.g. v4_C1jf: L5 61.50 / 29 targets + L2/3 40.10 / 9 targets).

| group (targets.py) | pathway | paper(s) | prep | n targets | data | VDCC traces | basis dir | joint |
|---|---|---|---|---|---|---|---|---|
| paired_l5: Markram | L5 TTPC -> L5 TTPC | Markram 1997 (10 Hz, -10/+5/+10) | paired | 3 | yes | markram_delta-prefire-vca | sabrina_n120_delta (batch_v2.BASIS_DIR) | yes (l5) |
| paired_l5: csv sjostrom_* | L5 -> L5 | Sjostrom 2001 (0.1-50 Hz, +-10) | paired | 10 | yes | ebner_delta-prefire-vca | same | yes (l5) |
| paired_l5: SJOSTROM_DRUGS | L5 -> L5 | Sjostrom 2003 (AM251, NMDAR block) | paired | 6 | yes | (same dirs) | same | yes (l5) |
| paired_l5: SJOSTROM_2003_TIMING | L5 -> L5 | Sjostrom 2003 Fig 9C (-25/-120/-200, burst r50) | paired | 5 | yes | sj03_ / sj03r50_delta-prefire-vca | same | yes (l5) |
| sjostrom07 | L5 -> L5 | Sjostrom 2007 (200 ms steps, AM251, L-NAME/cPTIO) | paired | 5 | yes | sj07_delta-prefire-vca | same | yes (l5) |
| paired_l23l5 | L2/3 PC -> L5 TTPC | Letzkus 2006, Sjostrom & Hausser 2006 | paired | 9 (7 with Letzkus 3AP -10 distal pair dropped) | yes | ebner_l23l5_delta-prefire-vseg-rs | ebner_l23l5_delta_rs (T29 rho0 re-split) | yes (l23) |
| paired_l23l23 | L2/3 PC -> L2/3 PC | Zilberter 2009 (Hardingham 2007 validation only) | paired | 15 (L23L23_TARGETS.md) | yes (2026-10-01) | pending: zilberter_l23l23 delta-prefire-vseg-rs (pilot 22128323) | pending: basis_results_edges_zilberter_l23l23_delta_rs (22128322) | no |
| nevian | L2/3 basal, extracellular | Nevian & Sakmann 2006 | extracellular | 22 rows (15 control, 7 MCPG/AM251/U73122/MK-801; combined per id) | yes | L5 cells only, no VDCC (ebner_delta-prefire); L2/3 = 3-post pipette pilot | none for L2/3 | no: validation only |
| ebner (csv, all 28 ids) | Ebner 2019 setup on L5 -> L5 | Nevian/Letzkus/Sjostrom rows | mixed | 28 | yes | Nevian/Letzkus ids on L5 have no VDCC | sabrina_n120 | no: wrong preparation for Nevian/Letzkus |
| markram tails | L5 -> L5 | anchor, not data | - | 4 | - | markram dirs | - | no |

Paired data not scored yet:

| data | pathway | n | what it needs |
|---|---|---|---|
| Letzkus no-post 0.98+-0.06, burst +500 ms 0.99+-0.05, 1AP unpaired 0.98+-0.04 | L2/3 -> L5 | 6 / 5 / 11 | new protocol ids in Ebner2019_L23PC_L5TTPC yaml, workdirs, prefire, -vseg-rs extraction |
| Letzkus NiCl2 distal +10 0.81+-0.01, -10 0.99+-0.10 | L2/3 -> L5 | fig | Ni (VDCC block) condition = new simulations |
| Sjostrom 2004 dLTD (pre + 250 ms subthreshold step) 0.69, AM251 ~1.0, ifenprodil ~1.1 | L5 -> L5 | 17 / 4 / 6 | simwriter step protocol, prefire, -vca extraction |
| Sjostrom 2003 Fig 9B curve (9 dt) | L5 -> L5 | - | read the 2001 PDF: induction frequency unclear |
| Sjostrom 2007 fixed 1.2 nA pair | L5 -> L5 | - | sensitivity check only (same cells), never a target |
| Zilberter 2009, Hardingham 2007 | L2/3 -> L2/3 | - | PDFs for numbers; then the full chain below |

## 2. What to add next (ranked) and what blocks it

No additional pathway is ready: only L5 -> L5 and L2/3 -> L5 have paired targets + VDCC traces + basis, and both are in
the joint fit already. fit_v4n.py --extra is ready for the next one.

1. **L2/3 -> L5 null controls (Letzkus no-post, +500 ms offset, 1AP unpaired; 3 targets).** Same pathway, basis and
   rho0 split exist. They constrain the excess LTP directly. Blocking: protocol ids in the yaml + workdirs (120 pairs),
   prefire (run_prefire_l23l5_rs.sh: measured 22089155 37:18, 189.7 GiB on 32 workers for 4 protocols; scale by
   protocols/steps), extraction (run_extract_l23l5_rs_vca.sh: 95G 4 CPU 0:15, measured 22090550). The extra targets
   join the existing l23 model (new rows in PAIRED_L23L5_EXTRA), not a new model.
2. **Sjostrom 2004 dLTD on L5 -> L5 (3 targets incl. AM251/ifenprodil arms).** Blocking: subthreshold-step protocol in
   simwriter (the Sj07 step machinery exists), prefire on the 24 subset pairs, -vca extraction; size from the sj07
   prefire/extraction (same pairs, 200k-step traces).
3. **Letzkus Ni distal (2 targets).** Blocking: a Ni condition needs new simulations with VDCC block (not an offline
   condition); same chain as 1.
4. **L2/3 -> L2/3 (new pathway model).** Blocking, in order: PDFs (Zilberter 2009, Hardingham 2007) for numbers ->
   targets in PAIRED_L23L23 -> yaml protocols + n120 workdirs -> edges with the conductance-based rho0 split (T29
   analog; no basis_results_edges_*l23l23* exists) -> cache (cf. run_cache_l23l5_rs.sh, 22085683: 5:11, 19.8 GiB on 16
   workers) -> EPSP basis (cf. run_basis_l23l5_rs.sh, 120 jobs 22086449-568: max 11.6 GiB, <= 14:49) -> prefire with v_seg
   (cf. 22089155) -> -vseg-rs extraction (cf. 22090550). Then add it with
   `--extra l23l23:paired_l23l23:<dirs>:basis_results_edges_<l23l23_rs>[:<geom csv>]`.
5. **Sjostrom 2003 Fig 9B.** Blocking: frequency from the 2001 PDF; then likely existing sj03 dirs.

Never fit: Nevian 2006 (extracellular, validation), Ebner csv Nevian/Letzkus rows on L5 -> L5 cells (wrong preparation),
Markram tail anchors (not data).

## 3. Using fit_v4n.py

`FITPY=fit_v4n.py JOINT=1 EXTRA="name:groups:dirs:basis[:geom[:pairs]]" sbatch run_fit_v4.sh` (several models: `;`
between specs). --joint = the l23 model exactly as in fit_v4.py. --drop-targets: `target|cond` for every model,
`name/target|cond` for one. Memory: scale the joint 2g measurement (24.1-25.9 GB for 502 + 334 records) by the extra
model's records x steps; GPU traces of L5 + L2/3 are ~11 GB, so a third pathway may need a larger MIG slice.
Reproduction check: job 22123123 (MAXITER=0, seed v4_C1jf; log prints `repro ...` and `REPRO OK (1e-9)`).
