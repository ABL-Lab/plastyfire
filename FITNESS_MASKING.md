# Dynamic Fitness Masking for Protocol Optimization

## Overview

This feature allows you to dynamically enable/disable specific protocols during different phases of the optimization process. This is useful when you want to:

1. Focus optimization on specific protocols during early generations
2. Expand to all protocols mid-optimization
3. Return to focused optimization in later generations

## How It Works

Fitness masking works by multiplying the error values for each protocol by a mask weight:
- Weight `0.0` = protocol does not contribute to fitness (effectively disabled)
- Weight `1.0` = protocol fully contributes to fitness (normal behavior)
- Weight `0.0 < w < 1.0` = protocol partially contributes to fitness

The mask is applied per generation based on a YAML configuration file.

## Usage

### 1. Create a Fitness Schedule YAML File

See `fitness_schedule_example.yaml` for a complete example:

```yaml
phases:
  # Phase 1: Generations 1-10 - optimize only mrk97_07 and mrk97_08
  - generation_start: 1
    generation_end: 11
    mask:
      mrk97_01: 0.0
      mrk97_02: 0.0
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 0.0

  # Phase 2: Generations 11-20 - optimize all protocols
  - generation_start: 11
    generation_end: 21
    mask:
      mrk97_01: 1.0
      mrk97_02: 1.0
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 1.0

  # Phase 3: Generations 21-30 - back to mrk97_07 and mrk97_08 only
  - generation_start: 21
    generation_end: 31
    mask:
      mrk97_01: 0.0
      mrk97_02: 0.0
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 0.0
```

### 2. Run Optimization with Fitness Schedule

Pass the fitness schedule file using the `--fitness-schedule` argument:

```bash
python -m plastyfire.modelfitter \
    --gen 30 \
    --pop_size 128 \
    --fitness-schedule fitness_schedule_example.yaml
```

### 3. Monitor the Logs

The logs will show when fitness masks are applied:

```
INFO - Starting Generation 1 ...
INFO - Generation 1: Applied fitness mask [0. 0. 1. 1. 0.] for protocols ['mrk97_01', 'mrk97_02', 'mrk97_07', 'mrk97_08', 'sjh06_02']
...
INFO - Starting Generation 11 ...
INFO - Generation 11: Applied fitness mask [1. 1. 1. 1. 1.] for protocols ['mrk97_01', 'mrk97_02', 'mrk97_07', 'mrk97_08', 'sjh06_02']
...
INFO - Starting Generation 21 ...
INFO - Generation 21: Applied fitness mask [0. 0. 1. 1. 0.] for protocols ['mrk97_01', 'mrk97_02', 'mrk97_07', 'mrk97_08', 'sjh06_02']
```

## Implementation Details

### Key Changes

1. **modelfitter.py**:
   - Added `--fitness-schedule` command-line argument
   - Enabled generation tracking via `generation_aware_map` function
   - Passes fitness schedule file to Evaluator

2. **evaluator.py**:
   - Added `fitness_schedule_file` parameter to `__init__`
   - Added `get_fitness_mask()` method to compute mask based on current generation
   - Modified `_process_evaluation_results()` to apply fitness mask to errors

### Mask Application Formula

```python
error = (np.array(error) * self.weight_reduce * fitness_mask).tolist()
```

Where:
- `error`: Raw error values per protocol
- `weight_reduce`: Standard protocol weights (unchanged from original)
- `fitness_mask`: Generation-specific mask from YAML schedule

## Advanced Usage

### Partial Masking

You can use weights between 0.0 and 1.0 for gradual transitions:

```yaml
phases:
  - generation_start: 1
    generation_end: 10
    mask:
      mrk97_01: 0.0
      mrk97_02: 0.0
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 0.0

  # Gradual transition phase
  - generation_start: 10
    generation_end: 12
    mask:
      mrk97_01: 0.5  # Half weight
      mrk97_02: 0.5  # Half weight
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 0.5  # Half weight
```

### No Fitness Schedule

If you don't provide a `--fitness-schedule` argument, all protocols will be optimized with equal weights throughout (default behavior, backward compatible).

## Troubleshooting

### Generation Counter Not Starting at 1

If resuming from a checkpoint, the generation counter continues from where it left off. Make sure your fitness schedule accounts for this.

### Mask Not Applied

Verify that:
1. The YAML file path is correct
2. Generation ranges in YAML match your `--gen` argument
3. Protocol IDs in YAML exactly match those in `PROTOCOL_IDX` in modelfitter.py

### All Protocols Masked to 0.0

This will result in zero fitness signal. Make sure at least one protocol has weight > 0.0 in each phase.
