# Nevian & Sakmann 2006 — Simulation Setup

Simulation framework for replicating the STDP experiments from:

> Nevian T, Sakmann B (2006) **Spine Ca²⁺ signaling in spike-timing-dependent plasticity.**  
> *J Neurosci* 26(43):11001–11013.

---

## Files

| File | Purpose |
|---|---|
| `nevian_protocols.py` | All 28 protocols as structured dataclasses (7 experiment groups) |
| `setup_nevian_simulations.py` | Creates workdir folders (spike trains, JSON configs, H5 files) |
| `submit_nevian.py` | Submits jobs to SLURM or runs locally |
| `Nevian and Sakmann - 2006 - ....md` | Paper transcription |

---

## Workflow

### Step 1 — Create simulation folders

```bash
# Dry run first (see what would be created)
python setup_nevian_simulations.py --protocols control --dry-run

# Create control protocol folders using existing L5TTPC pairs
python setup_nevian_simulations.py --protocols control

# All protocols (including blocker experiments)
python setup_nevian_simulations.py --protocols all

# Specific protocols only
python setup_nevian_simulations.py --protocols LTP_3ap_50hz_dt+10ms,LTD_3ap_50hz_dt-10ms

# Limit pairs for testing
python setup_nevian_simulations.py --protocols control --max-pairs 3 --dry-run
```

Folders are created under `nevian/simulations/{pre_gid}-{post_gid}/{protocol_id}/`.

### Step 2 — Submit / run simulations

```bash
# Submit to SLURM (default: v7 params)
python submit_nevian.py --protocols control

# Run locally with 8 workers
python submit_nevian.py --protocols control --execution-mode cpu --workers 8

# Different parameter set
python submit_nevian.py --params no_cicr --protocols control

# Test with 2 pairs
python submit_nevian.py --max-pairs 2 --protocols LTP_3ap_50hz_dt+10ms,LTD_3ap_50hz_dt-10ms
```

Results are saved to `nevian/results/{param_label}/{pair}/{protocol_id}/simulation_traces.pkl`.

---

## Protocol Groups

| Group | # Protocols | Key finding |
|---|---|---|
| Timing sweep (3AP 50Hz) | 6 | STDP curve shape |
| N-APs sweep (50Hz ±10ms) | 5 | ≥2 APs needed for LTP |
| Frequency sweep (3AP ±10ms) | 4 | >20Hz needed for LTP |
| NMDAR blockers | 4 | MK-801: LTP killed, LTD spared |
| VDCC blockers | 5 | Combined L+T block kills LTD |
| mGluR blockers | 3 | MCPG flips LTD→LTP at 100Hz |
| LTD cascade blockers | 3 | PLC, CB1 required; IP3R not required |

### Canonical protocols
- **LTP**: `LTP_3ap_50hz_dt+10ms` → expected ratio ≈ 2.01
- **LTD**: `LTD_3ap_50hz_dt-10ms` → expected ratio ≈ 0.68

---

## Notes

- Spike timing uses **Δt = time from EPSP onset to closest AP** (positive = AP after EPSP).
- `setup_nevian_simulations.py` reads current amplitude from the existing `single_cells/` cache.
- `pairrunner.py` and `simulator.py` are **not modified** — these scripts are drop-in compatible.
