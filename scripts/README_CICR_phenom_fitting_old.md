# Phenomenological CICR Model for Plasticity Fitting

## Motivation

With the Chindemi model alone (depression time constant changed from 670 to 365 ms), there is no parameter space that can simultaneously fit LTP and LTD across STDP protocols. The 10 Hz pre-post and post-pre calcium traces are too similar for threshold-based separation.

Biological evidence shows:

- **mGluR activation is necessary for LTD but not LTP.** Blocking mGluRs (with MCPG) converts a LTD protocol into LTP.
- During LTD protocols, mGluRs trigger IP3 production, which initiates CICR (Calcium-Induced Calcium Release) from the ER.
- CICR provides **moderate but prolonged** calcium elevation sufficient for depression.
- At steady state, SERCA balances IP3R release, producing a constant calcium "pedestal."

## Model Overview

The model inserts a **frequency-dependent calcium booster** between raw spine calcium (`cai_CR`) and the Graupner-Brunel plasticity rule. It has two stages:

```
cai_CR ──► mGluR Bucket ──► CICR Gate ──► ER Dynamics ──► cai_total ──► effcai ──► rho ODE ──► EPSP ratio
           (freq detect)     (on/off)     (SERCA+release)   (cai+pedestal)
```

### Stage 1: Frequency Detector ("Leaky Bucket")

A proxy for mGluR activation. Since we only have postsynaptic calcium traces, we assume that repeated synaptic activity (= repeated calcium spikes) implies mGluR activation by glutamate spillover.

**Event detection:** Rising-edge crossing of `cai_CR > ca_event_thresh`.

**Leaky integrator:**
```
mGluR(t + dt) = mGluR(t) * exp(-dt / tau_mGluR) + event(t)
```

Each detected spike adds +1 to the bucket. The bucket leaks with time constant `tau_mGluR`.

**Frequency discrimination:**
- **2 Hz:** Events are spaced 500 ms apart. The bucket drains between events. Level stays below threshold. No CICR.
- **10 Hz:** Events every 100 ms. Bucket fills faster than it drains. Level crosses threshold. CICR activates.
- **50 Hz:** Bucket fills instantly. CICR active, but the massive direct calcium influx dominates anyway.

**CICR gate:**
```
cicr_on = 1  if  mGluR > mGluR_thresh
           0  otherwise
```

### Stage 2: ER Dynamics ("Calcium Pedestal")

When the CICR gate is ON, the ER releases calcium into the cytoplasm. SERCA continuously pumps calcium back into the ER (when cytoplasmic calcium exceeds a per-synapse threshold). At balance, this produces a steady calcium pedestal.

**ER loading threshold (per synapse):**
```
gamma_er = a_er * Cpre + b_er * Cpost
```
SERCA is only active when `cai_total > gamma_er`. This makes ER loading depend on synapse-specific properties.

**SERCA pump (Hill kinetics, from cawave.py):**
```
J_serca = g_serca * cai_um^2 / (K_serca^2 + cai_um^2)    if cai_total > gamma_er
        = 0                                                 otherwise

where cai_um = 1000 * cai_total   (mM -> uM conversion)
      K_serca = 0.1 uM            (fixed, from cawave.py)
```

**CICR release:**
```
J_release = g_release * (ca_er - cai_total)    if cicr_on AND ca_er > cai_total
          = 0                                   otherwise
```

**ER conservation:**
```
d(ca_er)/dt = -(fc/fe) * (J_release - J_serca)
```
where `fc = 0.83` (cytoplasmic volume fraction), `fe = 0.17` (ER volume fraction).

**Cytoplasmic CICR contribution:**
```
d(ca_cicr)/dt = J_release - J_serca
ca_cicr >= 0  (clamped, cannot go negative)
```

**Total calcium:**
```
cai_total(t) = cai_CR(t) + ca_cicr(t)
```

This `cai_total` replaces `cai_CR` in the downstream effcai computation.

### Stage 3: Existing Pipeline (unchanged)

```
effcai(t+dt) = effcai(t) * exp(-dt/tau) + (cai_total(t) - min_ca) * tau * (1 - exp(-dt/tau))
```
followed by threshold crossing (theta_d, theta_p) and rho ODE integration.

## Parameters

### Existing (18 params, unchanged)

| # | Name | Description |
|---|------|-------------|
| 0-1 | `gamma_d`, `gamma_p` | Depression/potentiation rates |
| 2-5 | `a00, a01, a02, d0` | Basal theta_d coefficients |
| 6-9 | `a10, a11, a12, p0` | Basal theta_p coefficients |
| 10-13 | `a20, a21, a22, d0_ap` | Apical theta_d coefficients |
| 14-17 | `a30, a31, a32, p0_ap` | Apical theta_p coefficients |

### New CICR (7 params)

| # | Name | Default | Bounds | Source | Description |
|---|------|---------|--------|--------|-------------|
| 18 | `tau_mGluR` | 500 ms | [100, 5000] | ~ip3rtau | Leaky bucket time constant. Controls frequency selectivity. |
| 19 | `mGluR_thresh` | 3.0 | [0.5, 20] | fit | Bucket level that activates CICR. |
| 20 | `ca_event_thresh` | 0.0003 mM | [0.0001, 0.002] | CICR_phenom.mod | Calcium spike detection threshold. |
| 21 | `a_er` | 1.0 | [0, 5] | fit | Cpre coefficient for ER loading threshold. |
| 22 | `b_er` | 1.0 | [0, 5] | fit | Cpost coefficient for ER loading threshold. |
| 23 | `g_serca_cicr` | 1.9565 mM/ms | [0.01, 50] | cawave.cfg | SERCA maximal pump rate. |
| 24 | `g_release_cicr` | 0.5 /ms | [0.001, 10] | fit | CICR release rate constant. |

### Fixed Constants (from cawave.py / cawave.cfg)

| Name | Value | Description |
|------|-------|-------------|
| `K_serca` | 0.1 uM | SERCA Michaelis constant |
| `fc` | 0.83 | Cytoplasmic volume fraction |
| `fe` | 0.17 | ER volume fraction |
| `caAvg` | 0.0017 mM | Conserved total calcium (fc*cai + fe*ca_er) |

## Numerical Implementation

- **mGluR bucket:** Exact exponential integrator (stable for any dt).
- **ER dynamics:** Forward Euler with subcycling. Maximum substep `dt_sub = 0.1 ms` to ensure stability when g_release is large.
- **Event detection:** Rising-edge only (prevents counting one sustained calcium transient as multiple events).
- **Floors:** `ca_er >= 0`, `ca_cicr >= 0` enforced at each step.

## Why This Solves the LTD/LTP Problem

| Protocol | mGluR Bucket | CICR | Effect |
|----------|-------------|------|--------|
| 2 Hz, +5 ms | Low (drains between events) | OFF | No change. Calcium stays low. |
| 5 Hz, +5 ms | Borderline | OFF/weak | Minimal effect. |
| 10 Hz, +10 ms (LTP) | High | ON | Pedestal present, but large NMDA spikes push past theta_p anyway. LTP preserved. |
| 10 Hz, -10 ms (LTD) | High | ON | Pedestal lifts valleys between spikes above theta_d. Prolonged depression drive. LTD rescued. |
| 50 Hz, +10 ms | Very high | ON | Pedestal negligible compared to massive burst calcium. LTP preserved. |

## File Changes

All changes are in `fit_params.py`:

1. **Constants section:** Added CICR fixed constants (K_serca, fc, fe, caAvg).
2. **FIT_PARAMS / DEFAULT_PARAMS:** Added 7 CICR parameters (indices 18-24).
3. **`_apply_cicr_batch()`:** New numba JIT function implementing the full CICR model.
4. **`_load_pkl()`:** No longer pre-computes effcai. Stores raw `cai` only.
5. **`objective()` / `evaluate_params()`:** Now call `_apply_cicr_batch()` then `_compute_effcai_batch()` per evaluation, using the candidate CICR parameters.
6. **JIT warmup:** Warms up `_apply_cicr_batch` alongside existing functions.

## Usage

No change to CLI. The 7 new parameters are optimized alongside the existing 18:

```bash
python fit_params.py --method de --max-iter 2000 --workers 8
python fit_params.py --eval-only   # evaluate default params (including CICR defaults)
```
