# Summary of Changes - Dynamic Fitness Masking & Seed Individual

## Overview
Two new features have been added to the optimization system with minimal, surgical changes to the codebase:

1. **Dynamic Fitness Masking**: Ability to turn protocols on/off during different phases of optimization
2. **Seed Individual**: Ability to initialize the population with known good parameter values

## Changes Made

### 1. Fitness Masking Feature

#### Files Modified:
- `plastyfire/modelfitter.py`
- `plastyfire/evaluator.py`

#### Changes in `modelfitter.py`:
- **Line 86-91**: Added `--fitness-schedule` command-line argument
- **Line 157**: Pass `fitness_schedule_file` to Evaluator
- **Line 176-205**: Enabled generation tracking (uncommented existing code)

#### Changes in `evaluator.py`:
- **Line 307**: Added `fitness_schedule_file` parameter to `__init__`
- **Line 331-338**: Load and store fitness schedule from YAML file
- **Line 407-433**: New `get_fitness_mask()` method
- **Line 597-599**: Apply fitness mask to error calculation

### 2. Seed Individual Feature

#### Files Modified:
- `plastyfire/modelfitter.py`

#### Changes in `modelfitter.py`:
- **Line 92-97**: Added `--seed-individual` command-line argument
- **Line 223-235**: Parse seed values and create parent population

## Usage

### Fitness Masking

Create a YAML file (e.g., `fitness_schedule.yaml`):

```yaml
phases:
  - generation_start: 1
    generation_end: 11
    mask:
      mrk97_01: 0.0
      mrk97_02: 0.0
      mrk97_07: 1.0
      mrk97_08: 1.0
      sjh06_02: 0.0
```

Run with:
```bash
python -m plastyfire.modelfitter --fitness-schedule fitness_schedule.yaml
```

### Seed Individual

Run with:
```bash
python -m plastyfire.modelfitter \
    --seed-individual "199.887594498312,194.23843405523158,1.0005007769557286,2.706871658650942,2.5132952169332246,3.9826946063778337,9.594554277367934,3.0337323982995814,1.6205328327559538,2.692343865436747"
```

### Combined Usage

```bash
python -m plastyfire.modelfitter \
    --fitness-schedule fitness_schedule.yaml \
    --seed-individual "199.887594498312,194.23843405523158,1.0005007769557286,2.706871658650942,2.5132952169332246,3.9826946063778337,9.594554277367934,3.0337323982995814,1.6205328327559538,2.692343865436747"
```

## Files Created

1. `fitness_schedule_example.yaml` - Example fitness schedule configuration
2. `FITNESS_MASKING.md` - Detailed documentation for fitness masking feature
3. `test_fitness_schedule.py` - Test script to verify YAML parsing
4. `CHANGES_SUMMARY.md` - This file

## Key Design Decisions

### Surgical Precision
All changes were made with surgical precision to minimize impact:
- Fitness masking integrated into existing error calculation
- Seed individual uses existing `parent_population` parameter
- No changes to optimizer or core algorithm logic

### Backward Compatibility
- Without `--fitness-schedule`: All protocols optimized equally (default behavior)
- Without `--seed-individual`: Random initial population (default behavior)
- Existing optimizations continue to work unchanged

### Generation Tracking
- Enabled existing (commented-out) generation tracking code
- Required for fitness masking to know current generation
- Also provides useful timing information

## How It Works

### Fitness Masking
1. YAML file defines phases with generation ranges and protocol weights
2. `get_fitness_mask()` returns appropriate mask for current generation
3. Mask is multiplied with errors: `error = error * weight_reduce * fitness_mask`
4. Protocols with mask=0.0 contribute 0 to fitness (effectively disabled)
5. Protocols with mask=1.0 contribute fully to fitness

### Seed Individual
1. Parse comma-separated parameter values from command line
2. Call `opt.setup_deap()` to initialize DEAP framework
3. Create individual using `opt.toolbox.Individual()`
4. Replace random values with seed values
5. Pass as `parent_population=[seed_ind]` to optimizer
6. Optimizer includes seed in initial population alongside random individuals

## Testing

Run the fitness schedule test:
```bash
cd /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire
python test_fitness_schedule.py
```

Expected output shows masks changing across generations as defined in YAML.

## Log Messages

Look for these in logs to verify features are working:

**Fitness Schedule Loaded:**
```
INFO - Loaded fitness schedule: {'phases': [{'generation_start': 1, ...
```

**Seed Individual:**
```
INFO - Seeding initial population with: [199.887594498312, ...
```

**Fitness Mask Applied (per generation):**
```
INFO - Generation 1: Applied fitness mask [0. 0. 1. 1. 0.] for protocols ['mrk97_01', 'mrk97_02', 'mrk97_07', 'mrk97_08', 'sjh06_02']
```

## Notes

- The SLURM submission error in your logs is unrelated to these changes
- Both features work independently and can be used separately or together
- All changes preserve the existing optimization behavior when features are not used
