# The `c_pre`-only threshold model

Fitting the Graupner–Brunel STDP model to Markram et al. (1997) 10 Hz data for
L5TTPC→L5TTPC pairs, with thresholds driven by the presynaptic calcium
contribution alone. Compares against **DE fit #2**, the 6-parameter fit
documented in [`REPRODUCE.md`](REPRODUCE.md).

**Verdict:** the reduction does not hold up. Measured in bluecellulab over the
same 700 pair/protocol runs, weighted error at the three scored Δt goes
**1.61 → 10.87** (§2.3). The offline fit made it look 3× worse; the simulator
says 6.7×, concentrated at +5 and +10 ms. The offline model itself is the
surprise — it over-predicts potentiation by 0.06–0.09 in this parameter regime,
far outside `validate.py`'s 0.0011 gate.

---

## 1. The model

### 1.1 Calcium

The mechanism tracks a low-pass-filtered calcium, `effcai_GB`, driven by the
compartment's own `cai_CR` above a resting floor:

```
d(effcai)/dt = -effcai / tau_effca  +  (cai_CR - min_ca_CR)
```

with `tau_effca = 278.3177658387 ms` and `min_ca_CR = 70e-6 mM`. `extract.py`
rebuilds this offline from the recorded `cai_CR` at the native 0.025 ms step
using the exact exponential integrator

```
a = exp(-dt / tau_effca)
effcai[k+1] = a * effcai[k] + tau_effca * (1 - a) * (cai_CR[k] - min_ca_CR)
```

so no discretisation error is introduced relative to NEURON's own solve.

### 1.2 Thresholds — the full model

Each synapse carries two calcium contributions measured once, offline, and
cached: `c_pre` (the peak `effcai` deflection from a lone presynaptic spike) and
`c_post` (from a lone postsynaptic AP). The depression and potentiation
thresholds are linear in the two:

```
theta_d = a00 * c_pre  +  a01 * c_post
theta_p = a10 * c_pre  +  a11 * c_post
```

Apical synapses use `a20/a21` and `a30/a31` in place of `a00/a01` and `a10/a11`.
Every fit here **ties apical to basal** (`a20 = a00`, `a21 = a01`, `a30 = a10`,
`a31 = a11`), leaving 4 free a-parameters.

### 1.3 Thresholds — the `c_pre`-only reduction

Pin both postsynaptic coefficients to zero:

```
a01 = a11 = 0
```

so the thresholds collapse to a pure scaling of the presynaptic contribution:

```
theta_d = a00 * c_pre
theta_p = a10 * c_pre
```

**Two consequences follow immediately.**

First, the ordering constraint `theta_p > theta_d` — which in the full model
depends on each synapse's own `(c_pre, c_post)` and cannot be read off the
coefficients — reduces exactly to the scalar test

```
a10 > a00          (since c_pre >= 0 for every synapse)
```

Second, the thresholds become **rank-preserving in `c_pre`**: every synapse's
`theta_d` and `theta_p` are the same multiple of its own `c_pre`, so the ratio
`theta / c_pre` is now constant across the population. In the full model that
ratio varies synapse-by-synapse through `c_post`.

### 1.4 Thresholds under background — `c_bg`

Everything above describes the **in-vitro** case the fit targets: a quiescent
pair, calcium returning to the resting floor between pairings. In vivo the
synapse never returns to that floor, so `plastyfire`'s thresholds are shifted by
a measured background term before any in-vivo run uses them.

**`c_bg` is measured, not derived.** `wilson/step15_measure_background_1hz_intrinsic.py`
runs the postsynaptic cell with its DEES afferents at 1 Hz and the intrinsic
local circuit firing at in-vivo reference rates, records `effcai_GB` for all
17,477 DEES synapses for 10 s, and writes per-synapse `mean`, `p95` and `peak`.
Plasticity is genuinely frozen during that measurement: `theta_d_GB =
theta_p_GB = 1e9`. (Setting them *negative* does the opposite of disabling —
the mechanism uses `WATCH (effcai_GB > theta_d_GB)`, so a negative threshold
latches `dep_GB` and `pot_GB` on immediately and `rho` contaminates the very
trace being measured. That was a real bug in the superseded `step7`.)

**The shift.** `wilson/step16_write_corrected_edges.py` adds that floor to both
calcium contributions and re-applies the *same* a-parameters:

```
c_pre_bg  = c_pre  + c_bg          c_bg = per-synapse P95 of background effcai_GB
c_post_bg = c_post + c_bg

theta_d_bg = a00 * c_pre_bg + a01 * c_post_bg
theta_p_bg = a10 * c_pre_bg + a11 * c_post_bg
```

Because the mapping is linear and `c_bg` enters both terms, the correction is
exactly an **additive per-synapse offset** on each threshold:

```
theta_d_bg = theta_d + (a00 + a01) * c_bg
theta_p_bg = theta_p + (a10 + a11) * c_bg
```

(verified against the written table to 1.8e-15). For DE fit #2 those multipliers
are `a00 + a01 = 3.9060` and `a10 + a11 = 4.4094`.

**Under the `c_pre`-only reduction** `a01 = a11 = 0`, so the shift collapses to a
rescaling of a single quantity:

```
theta_d_bg = a00 * (c_pre + c_bg)          offset  a00 * c_bg = 1.0047 * c_bg
theta_p_bg = a10 * (c_pre + c_bg)          offset  a10 * c_bg = 1.9778 * c_bg
```

3.89x smaller an offset on `theta_d` than the full model, 2.23x smaller on
`theta_p`.

**What this does to the numbers** (17,477 DEES synapses on cell 186168, `defit2`
c_pre/c_post, `c_bg` from the 1 Hz intrinsic measurement):

| | median | IQR |
|---|---|---|
| `c_pre` | 0.0362 | 0.0275–0.0506 |
| `c_post` | 0.0000 | 0.0000–0.0001 |
| `c_bg` | 0.1159 | 0.0921–0.2067 |
| `theta_d` (uncorrected) | 0.0364 | 0.0277–0.0510 |
| `theta_d_bg` | 0.4937 | 0.3968–0.8780 |

`c_bg` is a median **3.25x larger than `c_pre`**, and the background term
contributes a median **92.7%** of `theta_d_bg` (IQR 89.5–95.1%). The
corrected threshold is mostly not about the synapse's own presynaptic drive.
The depression window widens with it, median `theta_p - theta_d` 0.0232 →
0.0840, a factor 3.63.

**Three consequences for §1.3's two claims.**

The ordering test `a10 > a00` **survives unchanged**: both corrected thresholds
scale the same non-negative quantity `c_pre + c_bg`, so `theta_p_bg > theta_d_bg`
everywhere iff `a10 > a00`. The reduction keeps that convenience under
background; the full model does not, where the ordering depends on
`(c_pre, c_post, c_bg)` per synapse.

Rank-preservation **does not survive as stated**. Under `c_pre`-only the
thresholds are rank-preserving in `c_pre + c_bg`, not in `c_pre`, and
`corr(c_bg, c_pre) = 0.49` — related, but not the same ordering. The constant
ratio does survive exactly: `theta_p_bg / theta_d_bg = a10 / a00 = 1.9685` for
every synapse, against an IQR of 1.154–1.182 in the full model, where the shared
`c_bg` dominates both thresholds and compresses their ratio toward 1.

And **`c_post` is even more vestigial here than §3 suggests.** The 22.1% quoted
there is the L5TTPC pair population; on these 17,477 in-vivo synapses **89.4%**
have `c_post` below 1e-4, and the median `c_post / c_pre` is 0.00104. Whatever
`a01` and `a11` are worth, in this population they are multiplying a number
~1000x smaller than the one `a00` and `a10` multiply.

### 1.5 The rho ODE

Synaptic efficacy `rho` follows the Graupner–Brunel bistable ODE:

```
tau_ind * d(rho)/dt = -rho * (1 - rho) * (rho_star - rho)
                      + pot * gamma_p * (1 - rho)
                      - dep * (1 - pot) * gamma_d * rho
```

with the indicator functions

```
dep = 1  if effcai > theta_d  else 0
pot = 1  if effcai > theta_p  else 0
```

and `rho_star = 0.5`, `tau_ind = 70 s`.

> **The `(1 - pot)` factor is not optional.** It gates depression *off* while
> potentiation is active, making the two mutually exclusive and giving a real
> depression *window* `(theta_d, theta_p)` rather than a half-line. This factor
> is present in the compiled mechanism the simulations actually load
> (`DEES_cell_packages/other_mods/modified_mechanisms/GluSynapse.mod`) but is
> **absent from `plastyfire/GluSynapse.mod`**. Omitting it makes every synapse
> depress above `theta_p` and biases `rho` low by ~0.13.

Because `tau_ind = 70 s` runs against a 42 s protocol, `rho` moves slowly and an
explicit step on the decimated grid is stable.

### 1.6 EPSP ratio

Final `rho` is **binarised at 0.5** and the EPSP is a linear superposition over
the synapses above the cut, read from the precomputed basis (mirroring
`evaluator_edges`):

```
EPSP = e0 + sum_i  b_i * (e_i - e0),      b_i = 1 if rho_i >= 0.5 else 0
ratio = (EPSP_after / EPSP_before) * (1 + cv2),   cv2 = min((sigma_0/e_0)^2, 0.25)
```

That binarisation is the reason the objective is **piecewise-constant** in the
a-parameters — flat almost everywhere, jumping only when a synapse crosses the
cut. No gradient exists, which is why the search is DE/CMA-ES rather than a
gradient method or a surrogate.

### 1.7 Objective

```
err = sum over dt in {-10, +5, +10}  ((pred_dt - target_dt) / SEM_dt)^2
```

plus, under `--anchor-tails`, a loose term pulling `dt = ±30, ±50` toward 1.0 at
SEM 0.15, and under `--sign` a penalty on wrong-sign plasticity. Targets are
Markram et al. (1997), 10 Hz:

| dt (ms) | target | SEM |
|---|---|---|
| −10 | 0.7922 | 0.0259 |
| +5 | 1.2038 | 0.0644 |
| +10 | 1.2013 | 0.0626 |

---

## 2. Results

![parameter comparison](params_defit2_vs_cpreonly.png)

![STDP curves](stdp_defit2_vs_cpreonly.png)

*(dark-surface variants: `params_defit2_vs_cpreonly_dark.png`,
`stdp_defit2_vs_cpreonly_dark.png`)*

### 2.1 Parameters

| | DE fit #2 (`c_pre`+`c_post`) | `c_pre`-only |
|---|---|---|
| `a00`  (`c_pre` → θ_d) | 1.003498 | 1.004709 |
| `a01`  (`c_post` → θ_d) | 2.902478 | 0 (pinned) |
| `a10`  (`c_pre` → θ_p) | 1.644558 | 1.977768 |
| `a11`  (`c_post` → θ_p) | 2.764812 | 0 (pinned) |
| `gamma_d`  (bound 50–200) | 77.7558 | 84.3845 |
| `gamma_p`  (bound 150–300) | 299.9121 (railed) | 186.5408 |
| `tau_effca` (ms) | 278.3178 | 278.3178 |
| free parameters | 6 | 4 |
| weighted error — offline | 2.594 | 7.718 |
| weighted error — **bluecellulab** | **1.61** | **10.87** |
| evaluations | 3888 | 2592 |

### 2.2 STDP curve

| dt (ms) | DE fit #2 | `c_pre`-only | in vitro |
|---|---|---|---|
| −50 | 0.8435 | 0.8594 | ≈ 1.0 |
| −30 | 0.8423 | 0.8402 | ≈ 1.0 |
| **−10** | **0.8065** | **0.8384** | 0.7922 ± 0.0259 |
| **+5** | **1.2376** | **1.1612** | 1.2038 ± 0.0644 |
| **+10** | **1.1549** | **1.1242** | 1.2013 ± 0.0626 |
| +30 | 0.9570 | 0.9885 | ≈ 1.0 |
| +50 | 0.8403 | 0.8883 | ≈ 1.0 |

Bold rows are the three fitted points. Tails are loosely anchored at 1.0
(SEM 0.15) and are not measured.

> **This table is not apples-to-apples.** The DE fit #2 column is
> *bluecellulab-measured*; the `c_pre`-only column is the *offline model's own
> prediction*. Section 2.3 replaces it with both columns measured the same way.

---

### 2.3 Measured in bluecellulab

Both parameter sets run through the full simulator over the same 700
pair/protocol workdirs (100 pairs × 7 Δt), same EPSP basis, same script
(`plot_stdp_ic_best.py`). c_pre-only = hash `bdbf06f915d0`, DE fit #2 =
`b8c7ff3ecf0a`.

![BCL comparison](stdp_bcl_defit2_vs_cpreonly.png)

| dt (ms) | in vitro | DE fit #2 (BCL) | `c_pre`-only (BCL) | `c_pre`-only (offline) | BCL − offline |
|---|---|---|---|---|---|
| −50 | ≈ 1.0 | 0.8468 ± 0.0370 | 0.8294 ± 0.0680 | 0.8594 | −0.0300 |
| −30 | ≈ 1.0 | 0.8456 ± 0.0370 | 0.8300 ± 0.0673 | 0.8402 | −0.0102 |
| **−10** | 0.7922 ± 0.0259 | **0.8077 ± 0.0371** | **0.8221 ± 0.0669** | 0.8384 | −0.0163 |
| **+5** | 1.2038 ± 0.0644 | **1.2196 ± 0.0767** | **1.0986 ± 0.0778** | 1.1612 | **−0.0626** |
| **+10** | 1.2013 ± 0.0626 | **1.1329 ± 0.0722** | **1.0373 ± 0.0750** | 1.1242 | **−0.0869** |
| +30 | ≈ 1.0 | 0.9569 ± 0.0494 | 0.9478 ± 0.0722 | 0.9885 | −0.0407 |
| +50 | ≈ 1.0 | 0.8474 ± 0.0380 | 0.8874 ± 0.0693 | 0.8883 | −0.0009 |

Weighted error over the three scored Δt, `sum(((pred - invitro)/SEM)^2)`:

| | −10 | +5 | +10 | total |
|---|---|---|---|---|
| DE fit #2 (BCL) | 0.358 | 0.060 | 1.194 | **1.61** |
| `c_pre`-only (BCL) | 1.333 | 2.668 | 6.863 | **10.87** |
| `c_pre`-only (offline) | 3.182 | 0.438 | 1.517 | 5.14 |

The offline total of 5.14 is the *scored* part of the 7.718 that `fit.py`
reports; the remaining 2.574 is the four `--anchor-tails` terms
(`sum(((tail - 1.0)/0.15)^2)`). 5.136 + 2.574 = 7.710, matching 7.718 to
transcription rounding.

**The offline model does not transfer to this parameter regime.** `validate.py`
gates the analytical model at mean |pred − obs| = 0.0011, but that gate was
established near the DE fit #2 parameters. In the `c_pre`-only regime the model
over-predicts potentiation by **0.063 at +5 and 0.087 at +10** — 60–80× the
validated error, and in the same direction at every Δt (the `BCL − offline`
column is negative in all seven). So the offline fit chose `a10`/`gamma_p` to
hit a potentiation level the simulator does not actually deliver.

This changes the verdict on the reduction. Judged offline it looked like a
3× error increase (2.594 → 7.718); judged in the simulator, where both sets are
measured identically, it is **6.7× worse** (1.61 → 10.87), and the damage is
almost entirely at +5 and +10 rather than at −10:

- **+5**: 1.2196 → 1.0986. DE fit #2 sits inside the in-vitro SEM bar
  (1.139–1.268); the reduction falls well below it.
- **+10**: 1.1329 → 1.0373, against a bar of 1.139–1.264. Both miss, but the
  reduction misses by 2.6 SEM instead of 1.1.
- **−10**: 0.8077 → 0.8221. Both are above the bar (0.766–0.818); the reduction
  is marginally worse but this was already the weakest constraint.

The tail improvement predicted offline (+30, +50 moving toward 1.0) **did not
survive measurement**: +30 came out 0.9478 vs DE fit #2's 0.9569, i.e. slightly
worse, not better. Only +50 improved (0.8874 vs 0.8474). The "buy" in section 3
was largely an artifact of the offline model, and the paragraphs there should be
read as describing the offline prediction, not the simulator.

Error bars on the `c_pre`-only curve are roughly **twice** DE fit #2's at the
tails (e.g. ±0.068 vs ±0.037 at −50). Constant `theta_p / theta_d = a10 / a00`
across the population (§1.3) removes the per-synapse spread that `c_post`
supplied, so pairs no longer average toward a common value — the population
response is more heterogeneous, not less.

---

## 3. What the reduction costs and buys

**Cost — the scored points.** The weighted error triples, 2.594 → 7.718, and all
of that is concentrated in the three fitted Δt. At −10 the reduced model lands
at 0.8384 against a SEM bar spanning 0.766–0.818: it misses in a direction the
data excludes, not merely noisily. At +5 and +10 it undershoots a wide bar.

**Buy — the right-hand tail.** +30 improves 0.9570 → 0.9885 and +50 improves
0.8403 → 0.8883, both moving toward the physiologically expected 1.0. This is
the "tails stay depressed" residual flagged in `REPRODUCE.md`. The improvement
is **not symmetric** — −50 gets slightly *worse* (0.8435 → 0.8594) — so it is
not simply "less plasticity everywhere"; the `c_post` term was doing something
specific on the pre-before-post side.

**`gamma_p` comes off its ceiling.** In the full model `gamma_p` rails at
299.91 of a 300 bound. Removing `c_post` from `theta_p` drops it to 186.54,
mid-range. The optimiser reached comparable potentiation by raising `a10`
(1.644558 → 1.977768) instead. So part of the ceiling pressure in DE fit #2 was
the `c_post` term forcing `theta_p` up, rather than a genuine demand for more
potentiation drive.

**`a00` is essentially unchanged** (1.003498 → 1.004709) and sits at its lower
bound in both fits. The depression threshold wants to be as low as the box
allows regardless of parameterisation.

### Why dropping `c_post` was worth testing at all

22.1% of synapses have `c_post` at the resting floor (< 1e-4), so their
depression window `(a11 - a01) * c_post` is ≈ 0 no matter what the parameters
are — a structural ceiling on achievable LTD. For that fifth of the population
the `c_post` term already carries no information.

---

## 4. Where the correction is used: the wilson ABCD setup

`c_bg` exists for one consumer, so what that consumer does is worth stating —
it is the only place these thresholds meet a non-quiescent cell.

**The cell and its synapses.** One postsynaptic L5 TTPC, gid 186168, from the
`O1_2023a_Ecker` S1nonbarrel circuit, carrying **17,477 DEES synapses** (the
plastic GluSynapse population) plus an intrinsic local circuit of 10,839
synapses from 1,648 presynaptic cells firing at in-vivo reference rates. The
dendrite is 366 branches; the cell spans ~1400 um from soma to tuft tip.

**The memory task.** `step1` assigns each of the 366 **branches** one of four
letters A/B/C/D from a seeded RNG, and a synapse inherits its branch's letter —
so a whole section shares a tuning, which is what makes "this branch is being
driven" a meaningful condition. Roughly half the synapses are marked `is_signal`
(8,641 signal / 8,836 background here). During a letter's presentation window,
the signal synapses carrying that letter fire at the **stimulus rate, 10 Hz**;
every other synapse stays at the **background rate, 1 Hz** — the same 1 Hz that
`step15` measured `c_bg` under.

**The schedule.** 2 cycles x 4 letters = 8 windows of 800 ms, separated by
1200 ms gaps, `tstop` 16.3 s. Each synapse therefore sees three distinct
conditions: its own letter's window (driven), another letter's window (its
branch quiet, but the soma still spiking and delivering bAP calcium), and the
gaps (background only).

**What it reads from here.** `theta_d_low_ca_bg` / `theta_p_low_ca_bg` out of
the edges file — the §1.4 corrected values — set on `hsyn.theta_d_GB` /
`theta_p_GB` at instantiation. `gamma_d`, `gamma_p` and `tau_effca` come from
the same fit-params JSON, so `defit2` in this document and `defit2` in the
wilson run are the same six numbers.

**Why the correction is not optional there.** Uncorrected `theta_d` has a median
of 0.0364 while the measured background floor `c_bg` has a median of 0.1159 —
background alone would sit a factor ~3 above the threshold at half the synapses,
latching `dep_GB` on permanently and driving the whole population to a fixed
point that has nothing to do with the task. The in-vitro fit never had to
confront this because the Markram pairing protocol leaves the cell quiescent
between pairings and `effcai` returns to the floor.

**What the three conditions actually look like.** Measured on a 5-synapse
cluster on `apic[48]` (apical tuft, ~1250 um, letter D), peak `effcai_GB`:

| condition | peak range | vs `theta_d_bg` (0.40–0.55) |
|---|---|---|
| its own letter (D) | 0.481–0.732 | crosses on all five |
| another letter (C) | 0.418–0.535 | straddles it |
| gap, background only | 0.010–0.023 | ~30–60x below |

The margin the branch's own stimulus buys over *another branch's* stimulus is
only 1.1–1.8x — the bAP calcium arriving during other letters does most of the
work — while the margin over true background is 30–60x. The corrected threshold
is placed inside that narrow first gap, which is why 92.7% of it being
background is not an accounting curiosity.

**What the `c_pre`-only reduction would do to it — untested.** The reduction
lowers `theta_d_bg` by a median factor 3.1 (0.4937 → 0.1570) while leaving the
calcium untouched, so far more synapses would cross far more often; the
depression window also widens (median 0.0840 → 0.1520). It is not a drop-in
substitution for an in-vivo run, and no such run has been done. Note also that
`c_bg` is condition-specific: it was measured at one background rate, on one
emodel, with the local circuit in one state. Any of those changing invalidates
the stored `theta_*_bg`, and folding it into a mapping fitted on quiescent
in-vitro pairs is an extrapolation of that mapping into a regime the fit never
saw.

---

## 5. Reproducing

The reduction is a flag on the existing fitter — no new files, no change to the
default path:

```bash
cd $PF/analytical_method

# timing + feasibility check
python fit.py --cpre-only --fit-gamma --time-eval --workers 60

# the fit: 4 free params (a00, a10, gamma_d, gamma_p)
python fit.py --cpre-only --fit-gamma --de --workers 60 --maxiter 80 \
              --anchor-tails --sign --verbose
```

`--cpre-only` sets `a01 = a11 = 0` and drops the search to 2 a-parameters, so DE
runs popsize 32 per generation instead of 48. Everything upstream — the
extracted calcium, the `c_pre`/`c_post` cache, the EPSP basis — is unchanged and
shared with the full fit.

If the feasibility block reports fewer than 30% of synapses active at the lower
bound, lower the floor (`theta_d = a00 * c_pre` is floored at `c_pre` itself when
`--a-lo 1.0`):

```bash
python fit.py --cpre-only --fit-gamma --de --workers 60 --maxiter 80 \
              --anchor-tails --sign --verbose --a-lo 0.3
```

### 5.1 bluecellulab

Same flag pattern: the parameter set lives in `run_de_fit2_pool.py` as
`CPRE_ONLY_ARGS` / `CPRE_ONLY_HASH`, selected with `--cpre-only`. The default
(DE fit #2) path is untouched.

```bash
cd $PF
sbatch submit_cpreonly_pool.sh                      # 700 sims, ~18 min on 60 cores
sbatch submit_cpreonly_pool.sh --skip-existing --force   # recover guardrail failures
sbatch submit_plot_cpreonly.sh                      # both curves + summaries
```

The first pass gave **673/700**; the 27 failures were all the induction
spike-count guardrail (52 post spikes instead of 50), clustered at short |Δt|.
`--force` recovers them, matching DE fit #2's 672 clean + 28 forced. The
guardrail fires on post-spike count during the pairing train, which the
a-parameters barely influence — hence 27 vs 28.

bluecellulab computes θ **at runtime** from the `c_pre`/`c_post` cache plus the
a-parameters, so `edges.h5` θ is not read on this path. Confirmed in the run
logs: `cq` (c_post) is present and non-zero, but θ_d = `a00 × c_pre` exactly.

`a01`/`a11`/`a21`/`a31` must be passed **explicitly as 0.0**, never omitted —
`simulator_edges._apply_theta_from_a_params` does `fit_params.get("a01", 1.0)`,
so a missing key silently becomes a coefficient of 1.0.

### 5.2 edges.h5 (for neurodamus)

Neurodamus has no runtime a-parameter path; it reads θ from `edges.h5`. Written
to a **copy**, because `inject_thresholds.py` rewrites the 4.2 GB file in place
and is not versioned — overwriting would destroy the DE fit #2 state that
`out_defit2/` was produced from.

```bash
cd $PF
cp data/dhuruva_modified_edges.h5 data/dhuruva_cpreonly_edges.h5

python compute_thresholds_from_cache.py \
    --cache cpre_cpost_cache/ion_channels_tau278.pkl \
    --output-dir threshold_results_cpreonly/ \
    --a00 1.004709 --a01 0.0 --a10 1.977768 --a11 0.0 \
    --a20 1.004709 --a21 0.0 --a30 1.977768 --a31 0.0 \
    --edges-h5 data/dhuruva_modified_edges.h5

python inject_thresholds.py --results-dir threshold_results_cpreonly/ \
    --edges-h5 data/dhuruva_cpreonly_edges.h5 --force
```

→ 732/732 synapses written, no `-1` sentinels, `theta_p > theta_d` everywhere,
and `theta_p / theta_d` a constant 1.968498 = `a10 / a00` — the rank-preserving
signature of §1.3, and a cheap way to verify the injection.

`data/dhuruva_cpreonly_circuit_config.json` points at the new edges file.
The neurodamus run itself has **not** been done.

---

## 6. Caveats

- **Resolved — both sets are now measured in bluecellulab** (§2.3, 700 runs
  each). Cite **1.61 vs 10.87**, not the offline 2.594 vs 7.718. Section 2.2 is
  kept only to show what the offline model predicted.
- **`validate.py`'s 0.0011 agreement does not hold here.** That gate was
  established near the DE fit #2 parameters. In the `c_pre`-only regime the
  offline model over-predicts potentiation by 0.063 at +5 and 0.087 at +10 —
  60–80× the validated error, and negative at all seven Δt. Re-validate before
  trusting the offline model anywhere far from DE fit #2.
- **The error comparison assumes matching flags.** 2.594 vs 7.718 is
  apples-to-apples only if DE fit #2's number also came from a run with
  `--anchor-tails --sign`. `REPRODUCE.md` records it under that command. The
  bluecellulab comparison in §2.3 has no such dependency — it scores measured
  curves directly.
- **Neurodamus has not been run** for these parameters. `edges.h5` is prepared
  (§5.2) but unused.
- **`data/dhuruva_modified_edges.h5` is stale for DE fit #2.** 18 of 732
  synapses sit at the `-1` sentinel that DE fit #2's own npz files say should
  carry real thresholds — the `c_pre`/`c_post` cache was rebuilt (Aug 19 17:27)
  after that injection (Aug 18 14:28), gaining 18 synapses that were never
  written. This is REPRODUCE.md's "714 synapses" vs today's 732. A `-1` θ is not
  neutral: `GluSynapse.mod` arms `WATCH (effcai_GB > theta_d_GB)`, so it fires
  immediately and pins `dep_GB = 1`. Affects `out_defit2/` (neurodamus) only —
  bluecellulab ignores `edges.h5` θ, so the 700 BCL pkls are unaffected. Not
  fixed here, since re-injecting would alter the state `out_defit2/` came from.
- **Apical is tied to basal** in both fits. An untied 8-parameter fit needs
  per-synapse section labels, which `extract.py` does not currently save.
- **Nothing here is fitted under background.** The a-parameters come from
  quiescent in-vitro pairs; `c_bg` is bolted on afterwards (§1.4) and only for
  the in-vivo consumer (§4). The `c_pre`-only parameters have never been run
  through that path.
