# Evaluation Script Usage Guide

The `evaluate_best_solution.py` script now supports evaluating multiple parameter sets and comparing them.

## Quick Start

### 1. Evaluate the Best Solution (Default)
```bash
python evaluate_best_solution.py --sample-size 100 --workers 60
```

### 2. Evaluate Custom Parameter Set (custom_v1)
```bash
python evaluate_best_solution.py --param-set custom_v1 --sample-size 100 --workers 60
```

### 3. Compare All Parameter Sets (Best + All Custom Sets)
```bash
python evaluate_best_solution.py --compare-all --sample-size 100 --workers 60
```

## Available Options

### Parameter Set Selection
- `--param-set {best,custom_v1}`: Choose which parameter set to evaluate
  - `best`: Use the optimized solution from bestsol.pkl or checkpoint.pkl
  - `custom_v1`: Use the custom parameters defined in the script

- `--compare-all`: Evaluate all parameter sets and create comparison plots/tables

### Performance Options
- `--sample-size N`: Number of pairs to evaluate (default: 100)
- `--workers N`: Number of parallel workers for EPSP computation (default: 60)
- `--sim-workers N`: Number of workers for simulations (default: same as --workers)
- `--ephys-workers N`: Number of workers for ephys generation (default: same as --workers)

### Other Options
- `--partial`: Skip ephys generation, use only existing ephys files
- `--seed N`: Random seed for reproducibility (default: 1234)

## Adding More Custom Parameter Sets

Edit the `CUSTOM_PARAM_SETS` dictionary in [evaluate_best_solution.py:65-80](evaluate_best_solution.py#L65-L80):

```python
CUSTOM_PARAM_SETS = {
    "custom_v1": {
        "gamma_d_GB_GluSynapse": 101.5,
        "gamma_p_GB_GluSynapse": 216.2,
        "a00": 1.002,
        "a01": 1.954,
        "a10": 1.159,
        "a11": 2.483,
        "a20": 1.127,
        "a21": 2.456,
        "a30": 5.236,
        "a31": 1.782,
    },
    "custom_v2": {  # Add your new parameter set here
        "gamma_d_GB_GluSynapse": ...,
        "gamma_p_GB_GluSynapse": ...,
        # ... rest of parameters
    },
}
```

## Output Structure

Results are saved to `/lustre06/project/6077694/dhuruva/plastyfire/evaluation_results/`:

### Single Parameter Set Evaluation
```
evaluation_results/
└── {param_set_name}/
    ├── detailed_results.csv          # Individual simulation results
    ├── comparison_summary.csv         # Protocol-level comparison with biodata
    ├── params.pkl                     # Parameters used
    └── evaluation_comparison.png      # Visualization
```

### Multi-Parameter Set Comparison
```
evaluation_results/
├── best/                              # Results for best solution
│   ├── detailed_results.csv
│   ├── comparison_summary.csv
│   ├── params.pkl
│   └── evaluation_comparison.png
├── custom_v1/                         # Results for custom_v1
│   ├── detailed_results.csv
│   ├── comparison_summary.csv
│   ├── params.pkl
│   └── evaluation_comparison.png
├── combined_comparison.csv            # Combined comparison table
└── combined_comparison.png            # Side-by-side comparison plot
```

## Example Commands

### Quick test with custom parameters (10 samples)
```bash
python evaluate_best_solution.py --param-set custom_v1 --sample-size 10 --workers 10
```

### Full evaluation comparing best vs custom
```bash
python evaluate_best_solution.py --compare-all --sample-size 100 --workers 60
```

### Resume evaluation with existing ephys files
```bash
python evaluate_best_solution.py --param-set custom_v1 --partial
```

## Understanding the Results

### comparison_summary.csv columns:
- `protocol_id`: Evaluation protocol (mrk97_03, mrk97_07, mrk97_08)
- `mean_epsp_ratio`: In vitro measured EPSP ratio
- `sem_epsp_ratio`: Standard error of in vitro measurement
- `mean_epsp_ratio_sim`: In silico simulated EPSP ratio (mean)
- `sem_epsp_ratio_sim`: Standard error of simulation
- `error`: Normalized error (difference / SEM)
- `param_set`: Which parameter set was used

Lower `error` values indicate better match with experimental data.
