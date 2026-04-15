# Analytical EPSP Ratio Prediction — Formula Derivation

> Companion document for `predict_epsp_ratios.py`

---

## Overview

Given the recipe parameters (`u`, `D`, `freq`) and fixed fit-params (`a00`–`a31`),
the predicted EPSP ratio after a plasticity induction protocol is:

```
EPSP_ratio ≈ 1 + p_ltp − 0.5 · p_ltd
```

where `p_ltp` and `p_ltd` are the fractions of synapses that flip rho state.
Everything flows from a single dimensionless number — the **plasticity efficiency ratio η**.

---

## Step 1 — Plasticity thresholds (from `simulator.py`)

Before induction, the simulator fires one isolated pre-spike and records the peak
`effcai` — that is **c_pre**. Same for a single post-spike → **c_post**.
The LTD and LTP thresholds are then set per synapse:

```
θ_d = a00 · c_pre + a01 · c_post    (LTD threshold)
θ_p = a10 · c_pre + a11 · c_post    (LTP threshold)
```

During induction:

| effcai region | flag | rho tendency |
|---|---|---|
| effcai > θ_p | `pot_GB = 1` | rho drifts **up** → LTP |
| θ_d < effcai < θ_p | `dep_GB = 1` | rho drifts **down** → LTD |
| effcai < θ_d | both = 0 | no change |

Default fit-params (Chindemi / Zenodo_O1):

| param | value |
|---|---|
| a00 | 1.0018 |
| a01 | 1.9536 |
| a10 | 1.1594 |
| a11 | 2.4828 |

---

## Step 2 — TM depression reduces Use during the train

The Tsodyks-Markram model depletes vesicle release probability over the pulse train.
At steady state (pulses ≳ 4):

```
Use_ss = u / (1 + u · D / isi)
```

where `isi = 1000 / freq` ms. Example at 10 Hz with u=0.5, D=672 ms:

```
Use_ss = 0.5 / (1 + 0.5 · 672 / 100) = 0.115
```

Release probability collapses from 0.50 → 0.115 by the 4th pulse and stays there.

---

## Step 3 — effcai accumulates across pulses

`effcai` decays exponentially with `τ_effca = 278 ms`. At 10 Hz (isi = 100 ms)
it does not fully decay between pulses, so a steady-state peak builds up.
The geometric series gives an **effective pulse count**:

```
N_eff = 1 / (1 − exp(−isi / τ_effca))
```

| freq (Hz) | isi (ms) | N_eff |
|---|---|---|
| 2 | 500 | 1.16 |
| 5 | 200 | 1.52 |
| 10 | 100 | 3.32 |
| 20 | 50 | 6.30 |

Peak `effcai` during induction ≈ `effcai_per_pulse × N_eff`.

---

## Step 4 — Coincidence amplification factor K

The c_pre calibration signal is measured with `Use = 1` and **no coincident post-spike**
(isolated pre-spike only). During the paired induction, the actual Ca²⁺ is much larger
because:

- The post-spike depolarises the spine → **NMDA Mg²⁺ unblock** → large Ca²⁺ through NMDA
- The back-AP re-opens **spine VDCCs**

So the actual induction effcai = **K** × (c_pre × Use_ss × N_eff), where K captures
the coincidence boost beyond what the single-spike calibration predicts.

**K ≈ 3.93** was measured empirically from 400 synapses in the `CHINDEMI_PARAMS_v2`
full-trace simulations (dt = +10 ms, 10 Hz):

```
K = (peak effcai during induction) / (c_pre × Use_ss × N_eff)   →  median = 3.93
```

---

## Step 5 — The plasticity efficiency ratio η

Dividing the induction effcai by θ_p:

```
η = K · c_pre · Use_ss · N_eff / (a10 · c_pre + a11 · c_post)
  = K · Use_ss · N_eff / (a10 + a11 · ζ)
```

where **ζ = c_post / c_pre ≈ 0.013** is the post-to-pre Ca²⁺ ratio.

> **Key simplification:** both c_pre and c_post scale similarly with spine volume
> (c_pre ∝ V⁻⁰·³⁶, c_post ∝ V⁻⁰·³⁷), so ζ is nearly constant and the volume
> dependence cancels. **η depends only on u, D, and freq** (plus the fixed a-params).

---

## Step 6 — Zone boundaries

```
η > 1                    → LTP zone
θ_d/θ_p  <  η  <  1     → LTD zone
η < θ_d/θ_p              → no plasticity
```

The LTD/LTP boundary ratio is fixed by the a-params:

```
θ_d/θ_p = (a00 + a01·ζ) / (a10 + a11·ζ)  ≈  0.862
```

So the **LTD window spans only [0.862, 1.0]** — a narrow 14% band.
η must land precisely in this range for LTD to occur.

### Operating window in terms of u

Solving η = 1 and η = θ_d/θ_p for u (at fixed D, freq):

```
u_ltp_min = 1 / (K · N_eff / (a10 + a11·ζ)  −  D / isi)
u_ltd_min = 1 / (K · N_eff / (a00 + a01·ζ)  −  D / isi)
```

LTD occurs when  `u_ltd_min ≤ u ≤ u_ltp_min`.

Example at 10 Hz, D = 672 ms:

| condition | u threshold |
|---|---|
| LTD lower bound (η = θ_d/θ_p) | u ≥ 0.169 |
| LTP lower bound (η = 1) | u ≥ 0.239 |

---

## Step 7 — EPSP ratio from the rho ODE

The rho variable evolves as (from `GluSynapse.mod`):

```
drho/dt = [ pot_GB · γ_p · (1 − rho)  −  dep_GB · γ_d · rho ] / τ_ind
```

with `τ_ind = 70 s` (very slow). The fraction of induction time (T = 140 s) with
`effcai > θ_p` is approximated by a sigmoid centred at η = 1:

```
f_ltp = σ(12 · (η − 1.02))          (sigmoid)

Δrho_ltp = γ_p · T / τ_ind · f_ltp  (clipped to [0, 1])
Δrho_ltd = γ_d · T / τ_ind · f_ltd
```

Default values: γ_p = 216.2,  γ_d = 101.5.

The flip probabilities are:

```
p_ltp = frac(rho=0) · Δrho_ltp    (depressed synapses flipping to potentiated)
p_ltd = frac(rho=1) · Δrho_ltd    (potentiated synapses flipping to depressed)
```

And the EPSP ratio:

```
EPSP_ratio ≈ 1  +  p_ltp  −  0.5 · p_ltd
```

The `+1` / `−0.5` weights come from the fact that the potentiated conductance
(gmax_p) is 2× the depressed conductance (gmax_d). A 0→1 flip doubles gmax
(+1 unit), a 1→0 flip halves it (−0.5 units).

---

## Full pipeline summary

```
recipe.csv
  u, D, gsyn, spinevol, correlations
        │
        ▼  Gaussian copula (mirrors epg.py)
  3000 sampled synapses  (u_i, D_i per synapse)
        │
        ▼
  Use_ss_i = u_i / (1 + u_i · D_i / isi)
  N_eff    = 1 / (1 − exp(−isi / τ_effca))
  η_i      = K · Use_ss_i · N_eff / (a10 + a11 · ζ)
        │
        ▼
  Zone assignment per synapse:
    η > 1            → p_ltp contribution
    θ_d/θ_p < η < 1  → p_ltd contribution
    η < θ_d/θ_p      → no contribution
        │
        ▼
  EPSP_ratio = mean(1 + p_ltp − 0.5 · p_ltd)
```

---

## Caveats and limitations

| limitation | effect |
|---|---|
| K = 3.93 calibrated at dt=+10ms only | dt=−10ms has slightly lower K (less NMDA coincidence); both dt values give same η in the script |
| ζ treated as constant | introduces ~5% error across the spinevol range 0.02–0.30 |
| sigmoid approximation for f_ltp | ignores within-train Use fluctuations; rho ODE shape is approximate |
| K calibrated at Chindemi params | if gsyn or spinevol changes drastically, K may shift |
| morphology variance ignored | actual c_pre has large scatter (p5–p95 spans 13×) not captured by recipe alone |

---

## Key insight: the operating window

The recipe parameters define whether plasticity is *possible* at all.
With fixed Chindemi a-params (a00–a31), the **only levers that move η** are:

1. **u** (Use0) — linear effect on Use_ss
2. **D** (tau_dep) — controls TM depression depth
3. **freq** — controls both Use_ss and N_eff

`spinevol`, `gsyn`, and `nrrp` affect **absolute Ca²⁺ levels** but their effects
cancel in η because the thresholds scale proportionally.
They matter only **indirectly** via the copula correlations
(high spinevol → high gsyn → tends to co-sample high u → higher η).

This means: **if you change a-params to fit a different circuit, you must also
ensure the recipe's (u, D) distribution places η near the LTD–LTP boundary**.
The two cannot be fitted independently.
