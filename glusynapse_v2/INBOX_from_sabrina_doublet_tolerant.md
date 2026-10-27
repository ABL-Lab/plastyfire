# From sabrina (2026-10-02): opt-in doublet-tolerant gates in plastyfire

The user asked to bring the doublet cells back into the og-GluSynapse Markram 10 Hz ±10 refit on delta-split2-bk25/bk17 "with a doublet-tolerant protocol". I read your DOUBLETS_FITTING.md. Your option (b) concern still applies, so the doublet-tolerant fit is reported next to the strict 69-pair fit, not instead of it.

## Code change (default off, so it is a no-op unless the env var is set)
Switch: `PLASTYFIRE_DOUBLET_TOLERANT=1`. The rule is ion_fitter's: per pulse, 1 AP or a doublet with ISI <= 6 ms, and no AP outside the pulses. Strict is always tried first, and tolerant is only a fallback.
- `plastyfire/simulator.py`:
  - new `DOUBLET_MAX_ISI_MS`, `doublet_tolerant()`, `pulse_spike_counts()`, `pulses_doublet_ok()`;
  - `spike_threshold_finder(..., doublet_ok=False)` dispatches to `_doublet_threshold_finder`. That function binary-searches on the number of pulses that fired. It stores `t_spikes` as the first AP of each pulse (all APs in `t_spikes_all`) and sets `doublet`.
- `plastyfire/simwriter.py`:
  - `calibrate_protocol_stimuli` retries cached None keys and failed strict searches with `doublet_ok=True`;
  - the older per-cell path has the same fallback for the 10 Hz train and single_ap.
- `plastyfire/simulator_edges.py`:
  - c_post: the tolerant fallback runs after all strict widths fail, so c_post = the Ca of one 1-AP or doublet event, and the multi-spike warning logs it;
  - induction guardrail: it passes if every pulse holds 1 AP or a <= 6 ms doublet, and the results get `doublet_tolerant_pass`.
- Unit check: `abl-sabrina-paper/experiments/launchers/og_stdp_refit/test_doublet_tol.py`, job 22261880, 9/9 pass.

Your uncommitted edits in these files are untouched; mine are additive.

## Update (14:30): scan up from threshold
The first version rescued 0/41 c_post cells. A spike diagnostic (`og_stdp_refit/diag_spikes.py`, job 22269220, bk25) found:
- 31/41 are bursters: 3–14 APs at any pulse.
- 10/41 give a slow doublet (ISI 11–31 ms) at threshold, but a single AP at a stronger 1.5 ms pulse (2.5–5 nA).

So `_doublet_threshold_finder` now goes on from the threshold amplitude. It scans upward with a stride (`TOLERANT_SCAN_STEPS` = 25) for the lowest amplitude where every pulse gives 1 AP or a ≤ 6 ms doublet, then refines inside the last stride. The result records `threshold_amp`.

The induction-train fallback in simwriter widens the grid with `tolerant_amp_grid` (up to 5 nA in 0.1 nA steps; the c_post grid already reaches 5 nA). The bursters still fail, as intended.

Unit check: job 22272035, 12/12 pass.

## Runs
They are in `abl-sabrina-paper/experiments/launchers/og_stdp_refit/` with `DT=1`, and outputs go to `/scratch/dhuruva/og_stdp_{bk25,bk17}/`:
- `extracted_dt`;
- `fit_results/*mk10dt*`;
- strict cache backup `sabrina_n120_<EM>_strict.pkl`.

## Outcome (bk25 only; the user locked in bk25)
- c_post: 10/41 rescued.
- Induction-train calibration: 5 of those 10 still have no 3 ms train.
- Prefire guardrail: most of the rest give 51–96 spikes for 50 pulses. Those are slow doublets or bursts, not ≤ 6 ms doublets.
- Net gain: 3 pairs, giving 72 pairs.
- DT fit (`fit_results/og-delta-split2-bk25-mk10dt_s1.json`): a00 1.020, a01 1.164, a10 1.270, a11 2.402. This is within the strict seeds' spread, so the strict 69-pair fit stays the reference.

## User decision (2026-10-02, later): the cell's own firing counts as one event
User: "if that is the behavior of the cell leave it as it is. when you do cpost consider doublet as one cell firing".
- New mode `PLASTYFIRE_DOUBLET_TOLERANT=2` (`simulator.any_firing()`):
  - the strict search runs first;
  - otherwise it takes the lowest amp at which every pulse fires, and c_post = the Ca of that whole event;
  - the induction guardrail logs extra APs instead of rejecting the run.
- Mode 1 is unchanged.
- Same rule in `abl-sabrina-paper/experiments/sabrina/dees.py` CpostSolver.
- Unit check: 22291980.
- In use for the bk25 og refit (all pairs) and Sabrina's branch set (23 branches).
