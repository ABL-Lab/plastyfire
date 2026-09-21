"""
Main run script for parameter optimisation using dhuruva_modified_edges.h5.

Replicates modelfitter.py but drives simulations via the prefire+edges pipeline
(EvaluatorEdges → pairrunner_edges_fit.py → simulator_edges.py) instead of the
recipe-based EPG pipeline.

Parameters optimised (tau_effca_GB_GluSynapse is fixed at FITTED_TAU):
  gamma_d_GB_GluSynapse, gamma_p_GB_GluSynapse, a00..a31

Usage:
    python modelfitter_edges.py \\
        --edges-h5 /path/to/dhuruva_modified_edges.h5 \\
        [--cpre-cpost-cache /path/to/cpre_cpost.pkl] \\
        [--sample_size 100] [--pop_size 128] [--gen 30]
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

import plastyfire.evaluator_edges as eval_edges
import plastyfire.simulator_edges as sim_mod

logger = logging.getLogger("modelfitter_edges")

CSVF_NAME = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/biodata/paired_recordings.csv"
# Only protocols with frequency_train=10 Hz are covered by the STDP simulation index
# (index_L5TTPC_L5TTPC_STDP.csv contains only frequency=10.0).
PROTOCOL_IDX = [
    "mrk97_07",  # freq=10 Hz, dt=+10 ms
    "mrk97_08",  # freq=10 Hz, dt=-10 ms
]
# Parameters to optimise (tau_effca is fixed; see evaluator_edges.FITTED_TAU)
FIT_PARAMS = [
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

# --same-apical-basal: drop the apical-specific coefficients from the search and
# mirror the basal ones into them (a20=a00, a21=a01, a30=a10, a31=a11). Halves
# the a-param search space from 8 dimensions to 4.
FIT_PARAMS_TIED = [
    ("gamma_d_GB_GluSynapse", 50.0, 200.0),
    ("gamma_p_GB_GluSynapse", 150.0, 300.0),
    ("a00", 1.0, 5.0),
    ("a01", 1.0, 5.0),
    ("a10", 1.0, 5.0),
    ("a11", 1.0, 5.0),
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit Graupner-Brunel plasticity model using dhuruva_modified_edges.h5"
    )
    parser.add_argument("--sample_size", type=int, default=100,
                        help="Number of in silico connections per protocol")
    parser.add_argument("--seed", type=int, default=1234, help="RNG master seed")
    parser.add_argument("-s", "--pop_size", type=int, default=128, help="Population size")
    parser.add_argument("-g", "--gen", type=int, default=30, help="Number of generations")
    parser.add_argument("-e", "--eta", type=float, default=20.0, help="Eta parameter")
    parser.add_argument("-m", "--mutpb", type=float, default=0.7, help="Mutation probability")
    parser.add_argument("-c", "--cxpb", type=float, default=0.3, help="Crossover probability")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose messages")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    parser.add_argument("--log-file", type=str, help="Custom log filename (without extension)")
    parser.add_argument("--max-jobs", type=int, default=900,
                        help="Maximum concurrent SLURM jobs (unused in multiprocessing mode)")
    parser.add_argument("--use-multiprocessing", action="store_true", default=True,
                        help="Use multiprocessing instead of SLURM (default: True)")
    parser.add_argument("--use-slurm", action="store_true",
                        help="Use SLURM batch jobs instead of multiprocessing")
    parser.add_argument("--fitness-schedule", type=str, default=None,
                        help="Path to fitness masking schedule YAML file")
    parser.add_argument("--same-apical-basal", action="store_true",
                        help="Force apical and basal a-params to be identical "
                             "(a20=a00, a21=a01, a30=a10, a31=a11). The GA then "
                             "searches 6 parameters instead of 10. --seed-individual "
                             "accepts either 6 values or the full 10 (apical entries "
                             "are dropped).")
    parser.add_argument("--seed-individual", type=str, default=None,
                        help="Comma-separated parameter values to seed initial population")
    parser.add_argument("--use-ipp", action="store_true",
                        help="Use IPyParallel for distributed execution")
    parser.add_argument("--profile", type=str, default="default",
                        help="IPyParallel profile name")
    # Edges-specific options
    parser.add_argument("--edges-h5", type=str, default=sim_mod.EDGES_H5_DEFAULT,
                        help="Path to dhuruva_modified_edges.h5")
    parser.add_argument("--cpre-cpost-cache", type=str, default=None,
                        help="Path to precomputed c_pre/c_post cache pkl")
    parser.add_argument("--basis-dir", type=str, required=True,
                        help="Directory containing basis_{pre}_{post}.csv files produced "
                             "by run_basis_pair_edges.py (e.g. basis_results_edges_mini/)")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Max parallel simulation workers (default: all CPUs)")
    parser.add_argument("--circuit-config", type=str, default=None,
                        help="Override the 'network' field in each prefire_simulation_config.json "
                             "(e.g. for new ion channels). Default: use path baked into each config.")
    args = parser.parse_args()

    # Configure logging
    os.makedirs("logs", exist_ok=True)
    log_level = logging.DEBUG if args.debug else logging.INFO
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(log_format)
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    for h in root_logger.handlers[:]:
        root_logger.removeHandler(h)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    log_filename = (
        f"logs/{args.log_file}.log" if args.log_file
        else f"logs/optimization_edges_{time.strftime('%Y%m%d_%H%M%S')}.log"
    )
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
    logger.info("Logging to: %s", log_filename)

    if args.debug:
        args.pop_size = 4
        args.gen = 3

    os.makedirs(".cache", exist_ok=True)

    # Load in vitro data
    invitro_db = pd.read_csv(CSVF_NAME)
    invitro_db = invitro_db.loc[invitro_db["protocol_id"].isin(PROTOCOL_IDX)]

    np.random.seed(args.seed)
    work_dir = os.path.abspath(os.getcwd())
    use_mp = args.use_multiprocessing and not args.use_slurm

    # Tie apical a-params to basal before the evaluator is built: the flag lives
    # on the evaluator module because that is where the GA vector is turned into
    # the 10 named parameters pairrunner_edges_fit.py expects.
    eval_edges.SAME_APICAL_BASAL = args.same_apical_basal
    fit_params = FIT_PARAMS_TIED if args.same_apical_basal else FIT_PARAMS
    if args.same_apical_basal:
        logger.info("--same-apical-basal: optimising %d params (a20=a00, a21=a01, "
                    "a30=a10, a31=a11)", len(fit_params))

    ev = eval_edges.EvaluatorEdges(
        fit_params,
        invitro_db,
        args.seed,
        args.sample_size,
        None,          # ipp_id (unused here)
        work_dir,
        args.max_jobs,
        use_multiprocessing=use_mp,
        max_workers=args.max_workers,
        fitness_schedule_file=args.fitness_schedule,
        edges_h5=args.edges_h5,
        cpre_cpost_cache=args.cpre_cpost_cache,
        basis_dir=args.basis_dir,
        circuit_config=args.circuit_config,
    )

    # Set map function
    # Architecture: sequential outer map (one individual at a time) + parallel inner
    # via ThreadPoolExecutor(--max-workers) in EvaluatorEdges.evaluate_with_multiprocessing.
    #
    # Why sequential outer?
    #   With a multiprocessing Pool outer, eval workers are daemon processes which
    #   historically caused issues. More importantly, with N CPUs it's more efficient
    #   to use all N CPUs for one individual's 200 sims (ThreadPoolExecutor) rather
    #   than splitting CPUs across multiple individuals simultaneously.
    #
    # Expected per-generation wall time:
    #   pop_size × ceil(200 / max_workers) × ~95 s/sim
    #   e.g. pop_size=16, max_workers=32 → 16 × 7 × 95 ≈ 2.8 h/gen
    if args.use_ipp:
        logger.info("Initialising IPyParallel client, profile=%s", args.profile)
        rc = ipp.Client(profile=args.profile)
        try:
            rc.wait_for_engines(1, timeout=60)
        except ipp.TimeoutError:
            logger.error("Timeout waiting for engines!")
            raise
        lview = rc.load_balanced_view()
        map_function = lview.map_sync
        logger.info("Connected to %d engines", len(rc))
    else:
        # Sequential outer — parallelism is handled inside evaluate_with_multiprocessing
        map_function = map

    logger.info("Optimisation parameters: eta=%f mut=%f cx=%f", args.eta, args.mutpb, args.cxpb)
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

    generation_counter = [0]
    generation_times = []

    def generation_aware_map(func, population):
        generation_counter[0] += 1
        ev.set_generation(generation_counter[0])
        t0 = time.time()
        logger.info("Starting generation %d ...", generation_counter[0])
        if args.use_ipp:
            results = list(map_function(ev.init_simulator_and_evaluate_with_lists, population))
        else:
            results = list(map(func, population))
        elapsed = time.time() - t0
        generation_times.append(elapsed)
        logger.info(
            "Generation %d done in %s (avg %.1f s/ind)",
            generation_counter[0],
            time.strftime("%H:%M:%S", time.gmtime(elapsed)),
            elapsed / len(population) if population else 0,
        )
        return results

    opt.map_function = generation_aware_map

    # Optionally seed initial population
    parent_population = None
    if args.seed_individual:
        seed_values = [float(x.strip()) for x in args.seed_individual.split(",")]
        if args.same_apical_basal and len(seed_values) == len(FIT_PARAMS):
            # Full 10-value seed given with tying on: keep the basal half and
            # drop the apical entries, which are now mirrored rather than free.
            # (Averaging the two halves instead would invent a value that was
            # never fitted, so prefer the basal set verbatim.)
            logger.info("--same-apical-basal: seed has %d values, using the basal "
                        "half %s and dropping apical %s",
                        len(seed_values), seed_values[:6], seed_values[6:])
            seed_values = seed_values[:6]
        if len(seed_values) != len(fit_params):
            raise ValueError(
                f"--seed-individual has {len(seed_values)} values but the search "
                f"space has {len(fit_params)} parameters "
                f"({[n for n, _, _ in fit_params]})")
        logger.info("Seeding initial population with: %s", seed_values)
        opt.setup_deap()
        # Generate a full population, then replace index 0 with the seed.
        # bluepyopt.run() overwrites offspring_size to len(parent_population),
        # so parent_population MUST have exactly pop_size individuals.
        parent_population = opt.toolbox.population(n=args.pop_size)
        for i, val in enumerate(seed_values):
            parent_population[0][i] = val
        logger.info("Seeded individual 0 with chindemi params; %d random individuals fill the rest",
                    args.pop_size - 1)

    # Run optimisation
    cpf_name = "checkpoint_edges.pkl"
    continue_cp = os.path.isfile(cpf_name)
    logger.info("%s optimisation", "Resuming" if continue_cp else "Starting new edges-based")
    pop, hof, log_obj, history = opt.run(
        max_ngen=args.gen,
        continue_cp=continue_cp,
        cp_filename=cpf_name,
        cp_frequency=1,
        parent_population=parent_population,
    )

    best = hof[np.argmin([np.linalg.norm(np.array(ind.fitness.values)) for ind in hof])]
    with open("bestsol_edges.pkl", "wb") as f:
        pickle.dump(ev.get_param_dict(best), f, -1)

    logger.info("Optimisation concluded — best solution written to bestsol_edges.pkl")

    if args.use_ipp:
        rc.close()
