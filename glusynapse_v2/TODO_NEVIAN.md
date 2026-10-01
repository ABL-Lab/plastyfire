# Side TODO: Nevian 2006 extracellular mGluR LTD (parked 2026-09-30)

Parked by the user ("do we need Nevian data? maybe not now"). Nevian is validation only, never a fit target.
Full analysis: `NEVIAN_EXTRACELLULAR.md`; gate check: `nevian_gate_check.py`.

## Why it matters
- The t_drive 2 fit's LTD gate is fully on for every Nevian post-before-pre protocol, including dt' -90 ms (data: no change).
  MCPG cannot block it in the model.
- Sjostrom 2003 tested LY341495 only on ACEA-induced depression (Fig 4F), not on the timing protocol,
  so L5 "mGluR-independent" does not rule out the pipette-driven mGluR LTD.

## Proposed design (not implemented)
- Pipette-only mGluR source: `T += A_M * g_i * S_M` at each pre arrival.
  - S_M: post-AP trace, tau 20 ms.
  - g_i: Gaussian of distance to the pipette tip, r_act 50 um; g = 0 for paired recordings.
- One free parameter, A_M. mglur_block (MCPG) sets A_M = 0; cb1_block sets A_mglu = 0.
- Known gap: 3AP50 dt' -10 (LTP 1.42) and 1AP -10 (LTD 0.80) have the same S_M (0.61).
  The source cannot separate them, so the LTP must come from the post/NO side. Add this as a unit-test case.

## Steps when resumed
1. Unit test on the login node (ordering 100 > 50 > 2AP ~ 20 > 1AP, none at -90, +dt unchanged, the case above).
2. Ask ebner_folder_setup for pipette mode (the session owns simwriter/prefire):
   - common ~1 ms delay;
   - per-synapse tip distance;
   - fix the EPSP-only rho drift (0.8 vs 0.97).
3. Pilot: 16 protocols x 3 posts, ~5.5 CPU.h (2 CPU / 19G / 1:25 per post). Offline 1-D A_M scan (<30 min).
4. Scale to 30 groups (~55 CPU.h, log first) only if the scan passes.
