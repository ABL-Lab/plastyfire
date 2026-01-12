"""
Main run script for parameter optimization that parses command line arguments, runs optimization, and saves results
authors: Giuseppe Chindemi (12.2020) + minor modifications by András Ecker (06.2024)
"""

import argparse
import logging
import multiprocessing
import os
import pickle
import time

import numpy as np
import pandas as pd
from bluepyopt.deapext.optimisations import IBEADEAPOptimisation
from deap.tools import ParetoFront
import ipyparallel as ipp

import plastyfire.evaluator as eval

logger = logging.getLogger("modelfitter")
CSVF_NAME = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/paired_recordings.csv"
PROTOCOL_IDX = [
    "mrk97_01",
    "mrk97_02",
    "mrk97_07",
    "mrk97_08",
    "sjh06_02",
]  # protocols to use for optimization
# PROTOCOL_IDX = ["mrk97_07", "mrk97_08"]  # protocols to use for optimization
# model parameters to be optimized (and their boundaries)
FIT_PARAMS = [  # ("tau_effca_GB_GluSynapse", 150., 350.),
    ("gamma_d_GB_GluSynapse", 50.0, 200.0),
    ("gamma_p_GB_GluSynapse", 150.0, 300.0),
    ("a00", 1.0, 5.0),
    ("a01", 1.0, 5.0),
    ("a10", 1.0, 5.0),
    ("a11", 1.0, 5.0),
    ("a20", 1.0, 10.0),
    ("a21", 1.0, 5.0),
    ("a30", 1.0, 10.0),
    ("a31", 1.0, 5.0),
]


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_size",
        type=int,
        default=100,
        help="Number of in silico connections per protocol",
    )
    parser.add_argument("--seed", type=int, default=1234, help="RNG master seed")
    parser.add_argument(
        "-s", "--pop_size", type=int, default=128, help="Population size"
    )
    parser.add_argument(
        "-g", "--gen", type=int, default=30, help="Number of generations"
    )
    parser.add_argument("-e", "--eta", type=float, default=20.0, help="Eta parameter")
    parser.add_argument(
        "-m", "--mutpb", type=float, default=0.7, help="Mutation probability"
    )
    parser.add_argument(
        "-c", "--cxpb", type=float, default=0.3, help="Crossover probability"
    )

    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true", help="Verbose messages"
    )
    parser.add_argument(
        "--debug", default=False, action="store_true", help="Enable debug mode"
    )
    parser.add_argument(
        "--log-file", type=str, help="Custom log filename (without extension)"
    )
    parser.add_argument(
        "--max-jobs", type=int, default=900, help="Maximum concurrent SLURM jobs"
    )
    parser.add_argument(
        "--use-multiprocessing",
        action="store_true",
        default=True,
        help="Use multiprocessing instead of SLURM",
    )
    parser.add_argument(
        "--use-slurm",
        action="store_true",
        help="Use SLURM batch jobs instead of multiprocessing",
    )
    parser.add_argument(
        "--fitness-schedule",
        type=str,
        default=None,
        help="Path to fitness masking schedule YAML file",
    )
    parser.add_argument(
        "--seed-individual",
        type=str,
        default=None,
        help="Comma-separated parameter values to seed initial population",
    )
    parser.add_argument(
        "--use-ipp",
        action="store_true",
        help="Use IPyParallel for distributed execution",
    )
    parser.add_argument(
        "--profile",
        type=str,
        default="default",
        help="IPyParallel profile name",
    )
    parser.add_argument(
        "--recipe-path",
        type=str,
        default=None,
        help="Path to custom recipe file",
    )
    args = parser.parse_args()
    # Configure logger
    # Create logs directory if it doesn't exist
    if not os.path.exists("logs"):
        os.makedirs("logs")

    # Configure logging to both file and console
    log_level = logging.DEBUG if args.debug else logging.INFO
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Create formatters
    formatter = logging.Formatter(log_format)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Clear any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File handler with custom name or timestamp
    if args.log_file:
        log_filename = f"logs/{args.log_file}.log"
    else:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_filename = f"logs/optimization_{timestamp}.log"

    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    logger.info(f"Logging to file: {log_filename}")
    logger.info(f"Log level: {'DEBUG' if args.debug else 'INFO'}")
    # Enable debugging mode
    if args.debug:
        eval.DEBUG = True
        # Also enable debug mode in simulator
        import plastyfire.simulator as sim

        sim.DEBUG = True
        args.pop_size = 4
        args.gen = 3
    # Set support directories
    if not os.path.exists(".cache"):
        os.makedirs(".cache")

    # Load in vitro results
    invitro_db = pd.read_csv(CSVF_NAME)
    invitro_db = invitro_db.loc[invitro_db["protocol_id"].isin(PROTOCOL_IDX)]
    # Create `bluepyopt` evaluator
    np.random.seed(args.seed)
    work_dir = os.path.abspath(os.getcwd())  # Get absolute current working directory
    # Determine execution method
    use_mp = args.use_multiprocessing and not args.use_slurm
    ev = eval.Evaluator(
        FIT_PARAMS,
        invitro_db,
        args.seed,
        args.sample_size,
        None,
        work_dir,
        args.max_jobs,
        use_multiprocessing=use_mp,
        fitness_schedule_file=args.fitness_schedule,
        recipe_path=args.recipe_path,
    )
    # Set map function
    if args.use_ipp:
        logger.info(f"Initializing IPyParallel client with profile: {args.profile}")
        rc = ipp.Client(profile=args.profile)
        
        # Wait for at least one engine to be ready
        logger.info("Waiting for engines to register...")
        try:
            rc.wait_for_engines(1, timeout=60)
        except ipp.TimeoutError:
            logger.error("Timeout waiting for engines! Check if ipengine is running.")
            raise

        dview = rc[:]
        # dview.use_dill()  # Dill not installed, relying on default pickle
        
        # Ensure NEURON is initialized on all engines
        logger.info("Initializing NEURON on all engines...")
        # We need to make sure the evaluator is available on engines if we were to use it directly,
        # but bluepyopt handles the mapping. However, we need to ensure the worker function
        # initializes the simulator.
        
        # Use load_balanced_view for better distribution of varying task lengths
        lview = rc.load_balanced_view()
        map_function = lview.map_sync
        logger.info(f"Connected to {len(rc)} engines")
    else:
        pool = multiprocessing.Pool(args.pop_size)
        map_function = pool.map

    # Create `bluepyopt` optimization
    logger.info(
        "Optimization parameters\nEta = %f Mut = %f Cx = %f"
        % (args.eta, args.mutpb, args.cxpb)
    )
    opt = IBEADEAPOptimisation(
        ev,
        offspring_size=args.pop_size,
        eta=args.eta,
        mutpb=args.mutpb,
        cxpb=args.cxpb,
        map_function=map_function,
        hof=ParetoFront(),
        seed=args.seed + 1,
    )

    # Custom map function to track generations with timing
    generation_counter = [0]
    generation_times = []

    def generation_aware_map(func, population):
        generation_counter[0] += 1
        generation_start_time = time.time()

        logger.info(f"Starting generation {generation_counter[0]} processing...")
        ev.set_generation(generation_counter[0])

        if args.use_ipp:
            # For IPP, we need to use the method that initializes the simulator
            # BluePyOpt expects the map function to take a function and a list of inputs.
            # The 'func' passed here is usually ev.evaluate_with_lists (or similar).
            # We want to force it to be ev.init_simulator_and_evaluate_with_lists
            # BUT, 'func' is bound to the evaluator instance method by BluePyOpt.
            
            # Actually, BluePyOpt calls map_function(self.evaluator.evaluate_with_lists, population)
            # We can override the function being mapped if we want, but it's cleaner to 
            # let BluePyOpt do its thing and just ensure the evaluator uses the right method.
            
            # However, since we can't easily change what BluePyOpt passes, we rely on 
            # the fact that we passed 'ev' to IBEADEAPOptimisation.
            # If we want to ensure initialization, we might need to wrap the function.
            
            # Let's try using the standard map first. If we need explicit initialization,
            # we might need to change how 'ev' behaves or use a wrapper.
            # Given the 'init_simulator_and_evaluate_with_lists' method exists in evaluator.py,
            # we should probably use it.
            
            # BluePyOpt's DEAP optimisation uses:
            # fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            # toolbox.evaluate is registered as evaluator.evaluate_with_lists
            
            # We can hack this by changing the registered function in the toolbox if we had access,
            # but here we are passing 'ev' to IBEADEAPOptimisation.
            
            # The cleanest way is to rely on 'ev' having the right method called.
            # IBEADEAPOptimisation registers 'evaluate' as 'evaluator.evaluate_with_lists'.
            
            # Wait, if we use IPP, the engines are persistent. They might need initialization ONCE.
            # But 'init_simulator_and_evaluate_with_lists' does it every time.
            # That's fine, it's safer.
            
            # To force this, we can monkey-patch the evaluator instance before passing it?
            # Or just update the map call here:
            
            results = map_function(ev.init_simulator_and_evaluate_with_lists, population)
        else:
            results = pool.map(func, population)

        generation_end_time = time.time()
        generation_duration = generation_end_time - generation_start_time
        generation_times.append(generation_duration)

        logger.info(
            f"Generation {generation_counter[0]} completed in {time.strftime('%H:%M:%S', time.gmtime(generation_duration))}"
        )
        if len(population) > 0:
            logger.info(
                f"Average time per individual: {generation_duration / len(population):.2f} seconds"
            )

        return results

    opt.map_function = generation_aware_map

    # Parse seed individual if provided and create parent population
    parent_population = None
    if args.seed_individual:
        seed_values = [float(x.strip()) for x in args.seed_individual.split(',')]
        logger.info(f"Seeding initial population with: {seed_values}")
        # Trigger DEAP setup to create Individual class
        opt.setup_deap()
        # Now create individual using toolbox
        seed_ind = opt.toolbox.Individual()
        # Replace its values with our seed values
        for i, val in enumerate(seed_values):
            seed_ind[i] = val
        parent_population = [seed_ind]

    # Run optimization
    cpf_name = "checkpoint.pkl"
    continue_cp = os.path.isfile(cpf_name)
    logger.info(
        "Resuming optimization" if continue_cp else "Starting a new optimization"
    )
    pop, hof, log, history = opt.run(
        max_ngen=args.gen, continue_cp=continue_cp, cp_filename=cpf_name, cp_frequency=1,
        parent_population=parent_population
    )
    # Gather and store best solution
    best = hof[np.argmin([np.linalg.norm(np.array(ind.fitness.values)) for ind in hof])]
    with open("bestsol.pkl", "wb") as f:
        pickle.dump(ev.get_param_dict(best), f, -1)

    logger.info("Optimization concluded")

    # Properly close the multiprocessing pool
    # Properly close the multiprocessing pool
    if not args.use_ipp:
        pool.close()
        pool.join()
