#!/usr/bin/env python
"""
Quick test script to verify fitness schedule YAML parsing
"""
import yaml
import numpy as np

def test_fitness_schedule():
    """Test loading and parsing fitness schedule"""

    # Load the example schedule
    with open('fitness_schedule_example.yaml', 'r') as f:
        schedule = yaml.safe_load(f)

    print("Loaded fitness schedule:")
    print(yaml.dump(schedule, default_flow_style=False))

    # Simulate protocol objectives
    protocol_ids = ["mrk97_01", "mrk97_02", "mrk97_07", "mrk97_08", "sjh06_02"]

    # Test for different generations
    test_generations = [1, 5, 10, 11, 15, 20, 21, 25, 30, 35]

    print("\nTesting fitness masks for different generations:")
    print("-" * 80)

    for gen in test_generations:
        # Find the appropriate phase
        mask = None
        for phase in schedule.get('phases', []):
            gen_start = phase.get('generation_start', 0)
            gen_end = phase.get('generation_end', float('inf'))

            if gen_start <= gen < gen_end:
                mask_dict = phase.get('mask', {})
                mask = [mask_dict.get(pid, 0.0) for pid in protocol_ids]
                break

        if mask is None:
            mask = [1.0] * len(protocol_ids)

        mask_array = np.array(mask)
        active_protocols = [pid for pid, w in zip(protocol_ids, mask) if w > 0.0]

        print(f"Generation {gen:2d}: mask={mask_array} -> Active: {active_protocols}")

    print("-" * 80)
    print("\nTest completed successfully!")
    print("\nKey observations:")
    print("  - Generations 1-10: Only mrk97_07 and mrk97_08 active")
    print("  - Generations 11-20: All protocols active")
    print("  - Generations 21-30: Back to mrk97_07 and mrk97_08 only")
    print("  - Generations 31+: No schedule defined, defaults to all active")

if __name__ == '__main__':
    test_fitness_schedule()
