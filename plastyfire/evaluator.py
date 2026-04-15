"""
Custom BluePyOpt evaluator for the Graupner & Brunel model
authors: Giuseppe Chindemi (12.2020) + minor modifications by András Ecker (06.2024)
"""

import glob
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import subprocess
import sys
import time
import traceback

import numpy as np
import pandas as pd
from bluepyopt.ephys.simulators import NrnSimulator
from bluepyopt.evaluators import Evaluator
from bluepyopt.objectives import Objective
from bluepyopt.parameters import Parameter

from plastyfire.config import OptConfig
from plastyfire.ephysutils import get_epsp_vector
from plastyfire.pyslurm import wait_for_job_slots

MIN2MS = 60 * 1000.0
FITTED_TAU = 278.3177658387  # previously optimized time constant of Ca*
# could use `SingletonWeightObjective`s, but it's easier to just multiply the (ordered) errors with the values below...
# WEIGHT_REDUCE will be set dynamically based on number of protocols in __init__
CONFIGS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/configs"
FIT_PARAM_NAMES = [
    "tau_effca_GB_GluSynapse",
    "gamma_d_GB_GluSynapse",
    "gamma_p_GB_GluSynapse",
    "a00",
    "a01",
    "a10",
    "a11",
    "a20",
    "a21",
    "a30",
    "a31",
]
logger = logging.getLogger(__name__)
DEBUG = False


def compute_epsp_prefire(pkl_file, window=100):
    """Safely compute EPSP values with bounds checking"""
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    logger.info(f"Data keys: {data.keys()}")
    t, v, spikes = data["t"], data["v"], data["pre_spikes"]

    # Extract rho values
    if len(data["rho_GB"]) > 100:
        data["rho_GB"] = np.transpose(data["rho_GB"])

    initial_rho = [0 if k[0] < 0.5 else 1 for k in data["rho_GB"]]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in data["rho_GB"]]

    valid_spikes = spikes[spikes < (t[-1] - window)][:60]

    epsp_values = get_epsp_vector(t, v, valid_spikes, window)
    avg_epsp = np.mean(epsp_values)

    return avg_epsp


def compute_epsp_ratio_prefire_batch(param_values, sim_dict, workdir, log_details=True, recipe_path=None):
    """Reads results from batch job output files and computes EPSP ratio"""
    import os
    import pickle

    # Reduced logging: only log if there are issues
    logger.info(
        f"Starting result processing for protocol {sim_dict['protocol_id']} in {workdir}"
    )

    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[
        :12
    ]  # Use first 12 chars

    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")

    if not os.path.exists(pkl_file):
        logger.error(f"Specific pickle file not found: {pkl_file}")
        raise FileNotFoundError(
            f"Simulation result file not found: simulation_{param_hash}.pkl"
        )

    logger.info(f"Reading results from: {pkl_file}")

    try:
        with open(pkl_file, "rb") as f:
            raw_results = pickle.load(f)
    except (pickle.UnpicklingError, EOFError, KeyError, ValueError) as e:
        logger.warning(
            f"Corrupted pickle file detected: {pkl_file} - {e}. Deleting and will re-run simulation."
        )
        os.remove(pkl_file)
        return None, None, []

    # get initial rho and final rho
    initial_rho = [0 if k[0] < 0.5 else 1 for k in raw_results["rho_GB"].T]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in raw_results["rho_GB"].T]
    initial_rho_str = "_".join(map(str, initial_rho))
    final_rho_str = "_".join(map(str, final_rho))

    pre_gid, post_gid = workdir.split("/")[-2].split("-")
    sim_config_path = os.path.join(workdir, "simulation_config.json")

    ephys_before_filename = f"ephys_data_{pre_gid}_{post_gid}_{initial_rho_str}.pkl"
    ephys_before_path = os.path.join(
        os.path.dirname(sim_config_path),
        "..",
        "..",
        "..",
        "..",
        "ephys_data",
        ephys_before_filename,
    )

    ephys_after_filename = f"ephys_data_{pre_gid}_{post_gid}_{final_rho_str}.pkl"
    ephys_after_path = os.path.join(
        os.path.dirname(sim_config_path),
        "..",
        "..",
        "..",
        "..",
        "ephys_data",
        ephys_after_filename,
    )

    # Return metadata about missing ephys files instead of generating immediately
    missing_ephys = []
    for rho_config, file_path in [
        (initial_rho_str, ephys_before_path),
        (final_rho_str, ephys_after_path),
    ]:
        if not os.path.exists(file_path):
            missing_ephys.append(
                {
                    "file_path": file_path,
                    "sim_config_path": sim_config_path,
                    "pre_gid": pre_gid,
                    "post_gid": post_gid,
                    "rho_config": rho_config,
                }
            )

    # If there are missing files, generate them immediately
    if missing_ephys:
        logger.info(
            f"Protocol {sim_dict['protocol_id']} needs {len(missing_ephys)} ephys files. Generating them now..."
        )
        from plastyfire.compute_test_pulse_epsp_ratio import get_epsp_value
        
        # Reconstruct fit_params
        fit_params = dict(zip(FIT_PARAM_NAMES, param_values))
        
        for ephys_info in missing_ephys:
            logger.info(f"Generating missing ephys file: {ephys_info['file_path']}")
            try:
                get_epsp_value(
                    ephys_info['sim_config_path'],
                    int(ephys_info['pre_gid']),
                    int(ephys_info['post_gid']),
                    ephys_info['rho_config'],
                    node_pop='S1nonbarrel_neurons', # Default
                    trial=0, # Default
                    output_dir=None,
                    no_cache=False,
                    recipe_path=recipe_path,
                    synapse_ids_str=None,
                    fit_params=fit_params
                )
            except Exception as e:
                logger.error(f"Failed to generate ephys file {ephys_info['file_path']}: {e}")
                return None, None, [] # Fail this simulation

    # If files exist (or were just generated), compute expected EPSP ratio
    if os.path.exists(ephys_before_path) and os.path.exists(ephys_after_path):
        epsp_before = compute_epsp_prefire(ephys_before_path)
        epsp_after = compute_epsp_prefire(ephys_after_path)
        epsp_ratio = epsp_after / epsp_before
    else:
        logger.error(
            f"Ephys files not found: {ephys_before_path} or {ephys_after_path}"
        )
        return None, None

    data = {
        "pkl_file": [pkl_file],
        "protocol_id": [sim_dict["protocol_id"]],
        "param_hash": [param_hash],
        "epsp_ratio": [epsp_ratio],
    }
    df = pd.DataFrame(data)

    results_file = os.path.join(os.path.dirname(__file__), "simulation_epsp_df.csv")
    if os.path.exists(results_file):
        df.to_csv(results_file, mode="a", header=False, index=False)
    else:
        df.to_csv(results_file, mode="w", header=True, index=False)

    logger.info(
        f"Successfully computed EPSP ratio: {epsp_ratio} for protocol {sim_dict['protocol_id']}"
    )
    return sim_dict["protocol_id"], epsp_ratio, []


def run_simulation_worker(args):
    """Worker function for multiprocessing that runs pairrunner.py directly"""
    if len(args) == 4:
        param_values, sim_dict, param_hash, recipe_path = args
    else:
        param_values, sim_dict, param_hash = args
        recipe_path = None

    try:
        # Build parameter arguments string
        param_dict = dict(zip(FIT_PARAM_NAMES, param_values))
        param_args = [f"--{name}={value}" for name, value in param_dict.items()]
        param_args.append(f"--fastforward={sim_dict['fastforward']}")
        param_args.append(f"--param_hash={param_hash}")

        # Run pairrunner.py directly
        script_path = os.path.join(os.path.dirname(__file__), "pairrunner.py")
        cmd = [sys.executable, script_path] + param_args

        logger.info(
            f"Running simulation for {sim_dict['protocol_id']} in {sim_dict['workdir']}"
        )
        result = subprocess.run(
            cmd, cwd=sim_dict["workdir"], capture_output=True, text=True, timeout=1800
        )

        if result.returncode != 0:
            logger.error(
                f"Simulation failed for {sim_dict['protocol_id']}: {result.stderr}"
            )
            return None

        # Check if PKL file was created
        pkl_file = os.path.join(sim_dict["workdir"], f"simulation_{param_hash}.pkl")
        if not os.path.exists(pkl_file):
            logger.error(f"PKL file not created for {sim_dict['protocol_id']}")
            return None

        # Process results
        return compute_epsp_ratio_prefire_batch(
            param_values, sim_dict, sim_dict["workdir"], log_details=True, recipe_path=recipe_path
        )

    except Exception as e:
        logger.error(f"Worker error for {sim_dict['protocol_id']}: {e}")
        return None


def compute_epsp_ratio_batch(param_values, sim_dict, workdir, log_details=True):
    """Reads results from batch job output files and computes EPSP ratio"""
    import os
    import pickle

    from plastyfire.ephysutils import Experiment

    # Reduced logging: only log if there are issues
    logger.info(
        f"Starting result processing for protocol {sim_dict['protocol_id']} in {workdir}"
    )

    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[
        :12
    ]  # Use first 12 chars

    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")

    if not os.path.exists(pkl_file):
        logger.error(f"Specific pickle file not found: {pkl_file}")
        raise FileNotFoundError(
            f"Simulation result file not found: simulation_{param_hash}.pkl"
        )

    logger.info(f"Reading results from: {pkl_file}")

    try:
        with open(pkl_file, "rb") as f:
            raw_results = pickle.load(f)
    except (pickle.UnpicklingError, EOFError, KeyError, ValueError) as e:
        logger.warning(
            f"Corrupted pickle file detected: {pkl_file} - {e}. Deleting and will re-run simulation."
        )
        os.remove(pkl_file)
        return None, None

    exp_handler = Experiment(
        raw_results,
        c01duration=sim_dict["c01duration"],
        c02duration=sim_dict["c02duration"],
        period=sim_dict["period"],
    )

    # Reduced logging: removed verbose debug
    logger.info(f"Computing EPSP ratio with nepsp={sim_dict['nepsp']}")
    epsp_ratio = exp_handler.compute_epsp_ratio(sim_dict["nepsp"])

    # Save data as pandas DataFrame (surgical addition - no return modification)
    data = {
        "pkl_file": [pkl_file],
        "protocol_id": [sim_dict["protocol_id"]],
        "param_hash": [param_hash],
        "epsp_ratio": [epsp_ratio],
    }
    df = pd.DataFrame(data)

    # Append to a single global results file (thread-safe approach)
    # Use absolute path to ensure consistent location
    results_file = os.path.join(os.path.dirname(__file__), "simulation_epsp_df.csv")
    if os.path.exists(results_file):
        df.to_csv(results_file, mode="a", header=False, index=False)
    else:
        df.to_csv(results_file, mode="w", header=True, index=False)

    # Reduced logging: only log essential info
    logger.info(
        f"Successfully computed EPSP ratio: {epsp_ratio} for protocol {sim_dict['protocol_id']}"
    )
    return sim_dict["protocol_id"], epsp_ratio


class Evaluator(Evaluator):
    """Graupner & Brunel plasticity model evaluator"""

    def __init__(
        self,
        fit_params,
        invitro_db,
        seed,
        sample_size,
        ipp_id,
        work_dir=None,
        max_jobs=900,
        use_multiprocessing=True,
        max_workers=None,
        fitness_schedule_file=None,
        recipe_path=None,
    ):
        """
        :param fit_params: list of tuples (name, min, max)
        :param invitro_db: pandas DataFrame with in vitro data
        :param seed: RNG seed
        :param sample_size: number of in silico connections per protocol
        :param ipp_id: IPyParallel client ID (optional)
        :param work_dir: working directory for batch scripts
        :param max_jobs: maximum concurrent SLURM jobs
        :param use_multiprocessing: use multiprocessing instead of SLURM
        :param max_workers: maximum number of workers for multiprocessing
        :param fitness_schedule_file: path to fitness schedule YAML file
        :param recipe_path: path to custom recipe file
        """
        super(Evaluator, self).__init__()
        self.params = [
            Parameter(param_name, bounds=(min_bound, max_bound))
            for param_name, min_bound, max_bound in fit_params
        ]
        self.param_names = [param.name for param in self.params]
        self.invitro_db = invitro_db
        self.seed = seed
        self.sample_size = sample_size
        self.ipp_id = ipp_id
        self.sim = NrnSimulator()
        self.current_generation = 0  # Track current generation
        self.work_dir = work_dir or os.getcwd()  # Working directory for batch scripts
        self.glost_dir = os.path.join(
            self.work_dir, "glost_jobs"
        )  # Directory for GLOST files
        os.makedirs(self.glost_dir, exist_ok=True)  # Create glost_jobs directory
        self._processed_files = set()  # Track processed files to prevent infinite loops
        self.max_jobs = max_jobs  # Maximum concurrent SLURM jobs
        self.use_multiprocessing = use_multiprocessing
        self.max_workers = max_workers or min(mp.cpu_count(), max_jobs)
        self.max_workers -= 4  # save some for other processes
        # Find all simulations
        self.all_sims, self.objectives = [], []
        # Store all protocol IDs to support dynamic masking
        self.all_protocol_ids = invitro_db["protocol_id"].unique().tolist()
        # Load fitness schedule if provided
        self.fitness_schedule = None
        if fitness_schedule_file:
            import yaml

            with open(fitness_schedule_file, "r") as f:
                self.fitness_schedule = yaml.safe_load(f)
            logger.info(f"Loaded fitness schedule: {self.fitness_schedule}")
        
        self.recipe_path = recipe_path
        if self.recipe_path:
            logger.info(f"Using custom recipe file: {self.recipe_path}")

        # Set weight reduce based on number of protocols
        num_protocols = len(invitro_db["protocol_id"].unique())
        if num_protocols <= 2:
            self.weight_reduce = np.array([1 / 8] * num_protocols)
        else:
            self.weight_reduce = np.array([1 / 8] * 2 + [1 / 4] * (num_protocols - 2))
        for elem in self.invitro_db.itertuples():
            # Add objective
            self.objectives.append(Objective(elem.protocol_id))
            # Load simulation config and extract simulation global parameters
            config = OptConfig(
                os.path.join(
                    CONFIGS_DIR, "%s_%s.yaml" % (elem.pre_mtype, elem.post_mtype)
                )
            )
            np.testing.assert_almost_equal(config.T / 1000.0, elem.period_sweep)
            fastforward = config.fastforward
            if fastforward is None:
                fastforward = config.C01_duration * MIN2MS + config.nreps * config.T
            # Ensure fastforward is a clean float to avoid precision issues
            fastforward = float(fastforward)
            nepsp = int(config.C01_duration * MIN2MS / config.T)
            # Load simulation index (witten by `simwriter.py`)
            sim_idx = pd.read_csv(
                os.path.join(
                    os.path.split(os.path.split(config.out_dir)[0])[0],
                    "index_%s_%s.csv" % (elem.pre_mtype, elem.post_mtype),
                )
            )
            sim_idx.set_index(["frequency", "dt"], inplace=True)
            sim_idx.sort_index(inplace=True)

            paths = sim_idx.loc[elem.frequency_train, elem.dt_train]["path"]
            # if paths is one element convert it inot pandas series
            if isinstance(paths, str):
                paths = pd.Series([paths])

            logger.info(f"paths: {paths}")
            # Sample paths based on sample_size
            if DEBUG:
                np.random.seed(self.seed)
                paths = paths.sample(min(3, len(paths)), random_state=self.seed)
            else:
                np.random.seed(self.seed)
                sample_size = min(self.sample_size, len(paths))
                paths = paths.sample(sample_size, random_state=self.seed)

            self.all_sims.extend(
                [
                    {
                        "protocol_id": elem.protocol_id,
                        "period": elem.period_sweep,
                        "c01duration": config.C01_duration,
                        "c02duration": config.C02_duration,
                        "fastforward": fastforward,
                        "nepsp": nepsp,
                        "simpath": path,
                        "workdir": "/".join(path.split("/")[:-1]),
                    }
                    for path in paths
                ]
            )

        logger.info(f"Found {len(self.all_sims)} simulations")
        # Reduced logging: removed verbose debug output
        # logger.debug("Available sims:")
        # for sim in self.all_sims:
        #     logger.debug(sim)

    def set_generation(self, generation):
        """Set current generation for batch job naming"""
        if generation != self.current_generation:
            logger.info(
                f"Starting Generation {generation} (Population size: {len(self.all_sims) if hasattr(self, 'all_sims') else 'unknown'} simulations per individual) ==="
            )
        self.current_generation = generation

    def get_fitness_mask(self):
        """Get fitness mask for current generation based on schedule"""
        if not self.fitness_schedule:
            # No schedule, return all 1.0 (all protocols active)
            return np.ones(len(self.objectives))

        # Find which phase we're in
        for phase in self.fitness_schedule.get("phases", []):
            gen_start = phase.get("generation_start", 0)
            gen_end = phase.get("generation_end", float("inf"))

            if gen_start <= self.current_generation < gen_end:
                # Get mask for this phase
                mask = phase.get("mask", {})
                fitness_mask = []

                for obj in self.objectives:
                    protocol_id = obj.name
                    weight = mask.get(protocol_id, 0.0)
                    fitness_mask.append(weight)

                fitness_mask = np.array(fitness_mask)
                logger.info(
                    f"Generation {self.current_generation}: Applied fitness mask {fitness_mask} for protocols {[obj.name for obj in self.objectives]}"
                )
                return fitness_mask

        # Default: all protocols active
        return np.ones(len(self.objectives))

    def get_param_dict(self, param_values):
        """Build dictionary of parameters for the Graupner & Brunel model
        from an ordered list of values (i.e. an individual)"""
        return dict(zip(self.param_names, param_values))

    def create_glost_task_list(self, param_values, all_sims):
        """Create GLOST task list for all simulations"""
        param_dict = self.get_param_dict(param_values)
        param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]

        param_args = " ".join(
            [f"--{name}={value}" for name, value in param_dict.items()]
        )

        task_list_path = os.path.join(self.glost_dir, f"glost_tasks_{param_hash}.txt")

        with open(task_list_path, "w") as f:
            for sim_dict in all_sims:
                cmd = f"cd {sim_dict['workdir']} && python /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plastyfire/pairrunner.py {param_args} --fastforward={float(sim_dict['fastforward'])} --param_hash={param_hash}"
                f.write(cmd + "\n")

        return task_list_path, param_hash

    def submit_glost_job(self, task_list_path, param_hash, num_tasks):
        """Submit GLOST job and return job ID"""
        wait_for_job_slots(max_jobs=self.max_jobs)

        template_path = os.path.join(
            os.path.dirname(__file__), "templates", "glost.batch.tmpl"
        )
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template not found: {template_path}")

        with open(template_path, "r") as f:
            template = f.read()

        # Determine resources based on task count
        # Narval nodes have 128 cores, using 63 to leave some headroom/system overhead or use half-nodes
        max_tasks_per_node = 63
        nodes = max(
            1, min(10, (num_tasks + max_tasks_per_node - 1) // max_tasks_per_node)
        )

        # Distribute tasks evenly across nodes if possible
        if num_tasks <= nodes * max_tasks_per_node:
            ntasks_per_node = (num_tasks + nodes - 1) // nodes
        else:
            ntasks_per_node = max_tasks_per_node

        # Calculate memory per node (approx 1.3GB per task based on seff analysis of 0.96GB avg usage)
        mem_per_node = f"{ntasks_per_node * 1300}M"

        template_vars = {
            "name": f"glost_{param_hash}",
            "nodes": nodes,
            "ntasks_per_node": ntasks_per_node,
            "mem_per_node": mem_per_node,
            "cpu_time": "01:50:00",
            "log_file": os.path.join(self.glost_dir, f"glost_{param_hash}.log"),
            "env": "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh",
            "task_list": task_list_path,
        }

        batch_script = template.format(**template_vars)
        batch_path = os.path.join(self.glost_dir, f"glost_{param_hash}.batch")

        with open(batch_path, "w") as f:
            f.write(batch_script)

        result = subprocess.run(
            ["sbatch", batch_path], capture_output=True, text=True, check=False
        )
        if result.returncode != 0:
            logger.error(f"sbatch failed with return code {result.returncode}")
            logger.error(f"sbatch stdout: {result.stdout}")
            logger.error(f"sbatch stderr: {result.stderr}")
            logger.error("Batch file contents:")
            with open(batch_path, "r") as f:
                logger.error(f.read())
            raise subprocess.CalledProcessError(
                result.returncode, ["sbatch", batch_path], result.stdout, result.stderr
            )

        output = result.stdout.strip()
        job_id = int(output.split()[-1])

        logger.info(
            f"Submitted GLOST job {job_id} with {num_tasks} tasks on {nodes} node(s)"
        )
        return job_id, batch_path

    def submit_ephys_job_array(self, missing_ephys_list):
        """Submit ephys generation as SLURM job array"""

        if not missing_ephys_list:
            return None

        # Deduplicate based on file_path
        unique_ephys = {}
        for ephys_info in missing_ephys_list:
            unique_ephys[ephys_info["file_path"]] = ephys_info

        missing_ephys_list = list(unique_ephys.values())
        logger.info(
            f"Submitting {len(missing_ephys_list)} unique ephys files as job array"
        )

        # Create ephys data directory
        ephys_dir = os.path.dirname(missing_ephys_list[0]["file_path"])
        os.makedirs(ephys_dir, exist_ok=True)

        # Create task file with one ephys job per line
        task_file = os.path.join(ephys_dir, f"ephys_tasks_{int(time.time())}.txt")
        with open(task_file, "w") as f:
            for ephys_info in missing_ephys_list:
                f.write(
                    f"{ephys_info['sim_config_path']} {ephys_info['pre_gid']} {ephys_info['post_gid']} {ephys_info['rho_config'].replace('_', ',')}\n"
                )

        # Create job array script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        setupenv_path = os.path.expanduser(
            "~/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh"
        )

        batch_script = f"""#!/bin/bash
#SBATCH --job-name=ephys_array
#SBATCH --account=ctb-emuller
#SBATCH --array=1-{len(missing_ephys_list)}
#SBATCH --cpus-per-task=4
#SBATCH --mem=4g
#SBATCH --time=01:00:00
#SBATCH --output={ephys_dir}/ephys_array_%A_%a.log

source {setupenv_path}

# Read the task line for this array task
TASK_LINE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {task_file})

# Parse arguments
read -r sim_config pre_gid post_gid rho_config <<< "$TASK_LINE"

# Run ephys generation
python {script_dir}/compute_test_pulse_epsp_ratio.py \\
    "$sim_config" \\
    "$pre_gid" \\
    "$post_gid" \\
    "$rho_config"
"""

        batch_path = os.path.join(ephys_dir, f"batch_ephys_array_{int(time.time())}.sh")
        with open(batch_path, "w") as f:
            f.write(batch_script)

        # Submit job array
        try:
            wait_for_job_slots(max_jobs=self.max_jobs)
            result = subprocess.run(
                ["sbatch", batch_path], capture_output=True, text=True, check=True
            )

            output = result.stdout.strip()
            job_id_match = output.split()[-1]
            job_id = int(job_id_match)

            logger.info(
                f"Submitted ephys job array {job_id} with {len(missing_ephys_list)} tasks"
            )
            return job_id

        except Exception as e:
            logger.error(f"Failed to submit ephys job array: {e}")
            return None

    def _process_evaluation_results(self, results, param_values):
        """Process evaluation results and compute errors - shared logic"""
        # Check if we have enough results
        if len(results) < len(self.objectives):
            logger.error(
                f"Insufficient results: got {len(results)}, need {len(self.objectives)}"
            )
            high_error = [1000.0] * len(self.objectives)
            logger.warning(
                f"Returning high error values due to insufficient results: {high_error}"
            )
            return high_error

        # Process results
        res_db = pd.DataFrame(results, columns=["protocol_id", "epsp_ratio"])

        logger.info("Simulations completed, results:")
        logger.info(res_db)
        insilico_db = (
            res_db.groupby("protocol_id")["epsp_ratio"]
            .agg(["mean", "sem"])
            .add_suffix("_epsp_ratio")
            .reset_index()
        )

        logger.info("Aggregating results")
        logger.info(insilico_db)

        # Compute error
        merged_db = pd.merge(
            self.invitro_db,
            insilico_db,
            on="protocol_id",
            suffixes=("_invitro", "_insilico"),
        )

        logger.info("Joining in vitro and in silico results")
        logger.info(merged_db)
        merged_db["error"] = np.abs(
            (
                merged_db["mean_epsp_ratio_invitro"]
                - merged_db["mean_epsp_ratio_insilico"]
            )
            / merged_db["sem_epsp_ratio_invitro"]
        )
        error = [
            float(merged_db.loc[merged_db["protocol_id"] == obj.name, "error"].iloc[0])
            for obj in self.objectives
        ]
        # Apply fitness mask based on generation
        fitness_mask = self.get_fitness_mask()
        error = (np.array(error) * self.weight_reduce * fitness_mask).tolist()

        logger.info("Sorting errors")
        logger.info(error)
        outcome = [
            float(
                merged_db.loc[
                    merged_db["protocol_id"] == obj.name, "mean_epsp_ratio_insilico"
                ].iloc[0]
            )
            for obj in self.objectives
        ]

        # Store results in cache
        cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
        pklf_name = os.path.join(".cache", "%s.pkl" % cachekey)
        with open(pklf_name, "wb") as f:
            pickle.dump(
                {
                    "error": error,
                    "outcome": outcome,
                    "individual": list(param_values),
                    "resdb": res_db,
                },
                f,
                -1,
            )

        return error

    def evaluate_with_multiprocessing(self, param_values):
        """Evaluate individual using multiprocessing instead of SLURM"""
        try:
            logger.info("Evaluating individual with multiprocessing: %s", param_values)

            # Check cache first
            cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
            pklf_name = os.path.join(".cache", "%s.pkl" % cachekey)
            if os.path.isfile(pklf_name):
                with open(pklf_name, "rb") as f:
                    cache_data = pickle.load(f)
                np.testing.assert_array_equal(param_values, cache_data["individual"])
                logger.info(f"Returning cached results for individual {cachekey[:12]}")

                if "resdb" in cache_data:
                    res_db = cache_data["resdb"]
                    current_protocol_ids = [obj.name for obj in self.objectives]
                    res_db_filtered = res_db[
                        res_db["protocol_id"].isin(current_protocol_ids)
                    ]

                    if len(res_db_filtered) == len(self.objectives):
                        results = [
                            (row["protocol_id"], row["epsp_ratio"])
                            for _, row in res_db_filtered.iterrows()
                        ]
                        error = self._process_evaluation_results(results, param_values)
                        return error
                    else:
                        logger.warning(
                            f"Cache missing protocols - need {len(self.objectives)}, have {len(res_db_filtered)}. Recomputing."
                        )
                else:
                    if len(cache_data["error"]) == len(self.objectives):
                        return cache_data["error"]
                    else:
                        logger.warning(
                            f"Cached error length mismatch - need {len(self.objectives)}, have {len(cache_data['error'])}. Recomputing."
                        )

            param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
            existing_results = []
            jobs_to_run = []

            # Check for existing results first
            for sim_dict in self.all_sims:
                pkl_file = os.path.join(
                    sim_dict["workdir"], f"simulation_{param_hash}.pkl"
                )
                file_key = f"{pkl_file}_{param_hash}"

                if file_key in self._processed_files:
                    logger.warning(
                        f"Already processed file {pkl_file} in this evaluation, skipping"
                    )
                    continue

                if os.path.exists(pkl_file):
                    logger.info(
                        f"Found existing PKL file for {sim_dict['protocol_id']}: {pkl_file}"
                    )
                    self._processed_files.add(file_key)
                    try:
                        result = compute_epsp_ratio_prefire_batch(
                            param_values,
                            sim_dict,
                            sim_dict["workdir"],
                            log_details=False,
                        )
                        if result is not None and result[1] is not None:
                            if len(result) == 3:
                                protocol_id, epsp_ratio, _ = result
                                if protocol_id is not None and epsp_ratio is not None:
                                    existing_results.append((protocol_id, epsp_ratio))
                            else:
                                existing_results.append(result)
                            continue
                        else:
                            logger.warning(
                                f"Invalid result from existing PKL file {pkl_file}, will run new simulation"
                            )
                    except Exception as e:
                        logger.info(f"Failed to load existing PKL file {pkl_file}: {e}")

                jobs_to_run.append((param_values, sim_dict, param_hash, self.recipe_path))

            # Run simulations - check if we're in a daemon process
            results = existing_results.copy()
            if jobs_to_run:
                # Check if current process is daemon (nested multiprocessing not allowed)
                current_process = mp.current_process()
                is_daemon = (
                    current_process.daemon
                    if hasattr(current_process, "daemon")
                    else False
                )

                if is_daemon:
                    logger.info(
                        f"Running {len(jobs_to_run)} simulations sequentially (daemon process detected)"
                    )
                    # Run sequentially in daemon processes
                    for job in jobs_to_run:
                        result = run_simulation_worker(job)
                        if result is not None:
                            results.append(result)
                else:
                    logger.info(
                        f"Running {len(jobs_to_run)} simulations using multiprocessing with {self.max_workers} workers"
                    )

                    with mp.Pool(processes=self.max_workers) as pool:
                        mp_results = pool.map(run_simulation_worker, jobs_to_run)

                    # Filter out None results and add to results
                    for result in mp_results:
                        if result is not None:
                            if len(result) == 3:
                                protocol_id, epsp_ratio, _ = result
                                if protocol_id is not None and epsp_ratio is not None:
                                    results.append((protocol_id, epsp_ratio))
                            else:
                                results.append(result)

            # Process results using shared logic
            error = self._process_evaluation_results(results, param_values)

            logger.info("Multiprocessing evaluation completed successfully")
            return error

        except Exception:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))

    def check_pkl_completion(self, workdir, param_hash):
        """Check if simulation completed by looking for the specific PKL file"""

        try:
            # Look for the specific PKL file pattern: simulation_{param_hash}.pkl
            pkl_pattern = os.path.join(workdir, f"simulation_{param_hash}.pkl")
            pkl_files = glob.glob(pkl_pattern)

            if pkl_files:
                # Check if the specific PKL file is complete (not being written)
                for pkl_file in pkl_files:
                    try:
                        # Check file size and stability
                        file_size = os.path.getsize(pkl_file)
                        if file_size > 0:
                            batch_file = os.path.join(
                                workdir, f"simulation_{param_hash}.batch"
                            )
                            if os.path.exists(batch_file):
                                try:
                                    os.remove(batch_file)
                                except FileNotFoundError:
                                    pass
                            return True
                    except (OSError, IOError):
                        # File might be being written, continue checking
                        continue

            return False

        except Exception:
            return False

    def _log_job_failure_details(self, job_id, job_info):
        """Extract and log detailed failure information from job status"""
        try:
            # Parse job info to extract failure details
            lines = job_info.split("|")
            for line in lines:
                if "=" in line:
                    key, value = line.split("=", 1)
                    if key.strip() in [
                        "JobState",
                        "ExitCode",
                        "DerivedExitCode",
                        "Reason",
                        "Comment",
                    ]:
                        logger.error(f"Job {job_id} {key.strip()}: {value.strip()}")

            # Try to get job history for more details
            hist_result = subprocess.run(
                [
                    "sacct",
                    "-j",
                    str(job_id),
                    "-o",
                    "JobID,JobName,State,ExitCode,DerivedExitCode,Reason,MaxRSS,MaxVMSize,Elapsed",
                    "--noheader",
                ],
                capture_output=True,
                text=True,
            )
            if hist_result.returncode == 0 and hist_result.stdout.strip():
                logger.error(f"Job {job_id} detailed history:")
                for line in hist_result.stdout.strip().split("\n"):
                    if line.strip():
                        logger.error(f"  {line.strip()}")

            # Check if there are any output/error files
            self._check_job_output_files(job_id)

        except Exception as e:
            logger.error(f"Failed to parse failure details for job {job_id}: {e}")

    def _check_job_output_files(self, job_id):
        """Check for job output files that might contain error information"""
        try:
            # Look for SLURM output files in current directory
            for filename in os.listdir("."):
                if filename.startswith(f"slurm-{job_id}"):
                    logger.info(f"Found SLURM output file: {filename}")
                    try:
                        with open(filename, "r") as f:
                            content = f.read()
                            # Log last few lines that might contain error info
                            lines = content.strip().split("\n")
                            if len(lines) > 10:
                                logger.error(f"Last 10 lines of {filename}:")
                                for line in lines[-10:]:
                                    logger.error(f"  {line}")
                            else:
                                logger.error(f"Content of {filename}:")
                                for line in lines:
                                    logger.error(f"  {line}")
                    except Exception as e:
                        logger.error(f"Could not read {filename}: {e}")
        except Exception as e:
            logger.error(f"Error checking job output files: {e}")

    def evaluate_with_lists(self, param_values):
        """Evaluate individual using batch jobs or multiprocessing"""
        if self.use_multiprocessing:
            return self.evaluate_with_multiprocessing(param_values)

        try:
            # Reduced logging: only log essential evaluations
            logger.info("Evaluating individual: %s", param_values)

            # Add a simple check to prevent infinite loops
            if hasattr(self, "_evaluation_count"):
                self._evaluation_count += 1
                if self._evaluation_count > 1000:  # Prevent infinite loops
                    logger.error(
                        "Too many evaluations detected, possible infinite loop. Exiting."
                    )
                    return [1000.0] * len(self.objectives)
            else:
                self._evaluation_count = 1

            cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
            pklf_name = os.path.join(".cache", "%s.pkl" % cachekey)
            if os.path.isfile(pklf_name):
                with open(pklf_name, "rb") as f:
                    cache_data = pickle.load(f)
                np.testing.assert_array_equal(
                    param_values, cache_data["individual"]
                )  # verify no collision (OMG)
                logger.info(f"Returning cached results for individual {cachekey[:12]}")

                param_hash = cachekey[:12]
                if "resdb" in cache_data:
                    res_db = cache_data["resdb"]
                    # Filter resdb to only include current protocols
                    current_protocol_ids = [obj.name for obj in self.objectives]
                    res_db_filtered = res_db[
                        res_db["protocol_id"].isin(current_protocol_ids)
                    ]

                    # Check if we have results for all current protocols
                    if len(res_db_filtered) == len(self.objectives):
                        # Recompute error with current weight_reduce
                        results = [
                            (row["protocol_id"], row["epsp_ratio"])
                            for _, row in res_db_filtered.iterrows()
                        ]
                        error = self._process_evaluation_results(results, param_values)

                        # Save filtered results to CSV
                        cached_data = []
                        for _, row in res_db_filtered.iterrows():
                            cached_data.append(
                                {
                                    "pkl_file": f"cached_{param_hash}_{row['protocol_id']}.pkl",
                                    "protocol_id": row["protocol_id"],
                                    "param_hash": param_hash,
                                    "epsp_ratio": row["epsp_ratio"],
                                }
                            )
                        if cached_data:
                            cached_df = pd.DataFrame(cached_data)
                            results_file = os.path.join(
                                os.path.dirname(__file__), "simulation_epsp_df.csv"
                            )
                            if os.path.exists(results_file):
                                cached_df.to_csv(
                                    results_file, mode="a", header=False, index=False
                                )
                            else:
                                cached_df.to_csv(
                                    results_file, mode="w", header=True, index=False
                                )
                        return error
                    else:
                        logger.warning(
                            f"Cache missing protocols - need {len(self.objectives)}, have {len(res_db_filtered)}. Recomputing."
                        )
                else:
                    # Old cache format, just return it (assuming it matches)
                    if len(cache_data["error"]) == len(self.objectives):
                        return cache_data["error"]
                    else:
                        logger.warning(
                            f"Cached error length mismatch - need {len(self.objectives)}, have {len(cache_data['error'])}. Recomputing."
                        )

            param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]
            existing_results = []
            sims_to_run = []

            for sim_dict in self.all_sims:
                pkl_file = os.path.join(
                    sim_dict["workdir"], f"simulation_{param_hash}.pkl"
                )
                file_key = f"{pkl_file}_{param_hash}"

                if file_key in self._processed_files:
                    logger.warning(
                        f"Already processed file {pkl_file} in this evaluation, skipping to prevent infinite loop"
                    )
                    continue

                if os.path.exists(pkl_file):
                    logger.info(
                        f"Found existing PKL file for {sim_dict['protocol_id']}: {pkl_file}"
                    )
                    self._processed_files.add(file_key)
                    try:
                        result = compute_epsp_ratio_prefire_batch(
                            param_values,
                            sim_dict,
                            sim_dict["workdir"],
                            log_details=False,
                        )
                        if result is not None and result[1] is not None:
                            existing_results.append(result)
                            continue
                        else:
                            logger.warning(
                                f"Invalid result from existing PKL file {pkl_file}, will run via GLOST"
                            )
                    except Exception as e:
                        logger.info(f"Failed to load existing PKL file {pkl_file}: {e}")
                else:
                    logger.info(
                        f"No existing PKL file found for {sim_dict['protocol_id']}: {pkl_file}"
                    )

                sims_to_run.append(sim_dict)

            glost_job_id = None
            batch_path = None
            successful_jobs = []

            if sims_to_run:
                task_list_path, param_hash = self.create_glost_task_list(
                    param_values, sims_to_run
                )
                glost_job_id, batch_path = self.submit_glost_job(
                    task_list_path, param_hash, len(sims_to_run)
                )

                logger.info(f"Waiting for GLOST job {glost_job_id} to complete...")

                max_wait_time = 7200
                poll_interval = 10
                start_time = time.time()

                while (time.time() - start_time) < max_wait_time:
                    all_complete = True
                    for sim_dict in sims_to_run:
                        if not self.check_pkl_completion(
                            sim_dict["workdir"], param_hash
                        ):
                            all_complete = False
                            break

                    if all_complete:
                        logger.info(f"GLOST job {glost_job_id} completed successfully")
                        successful_jobs = [
                            (glost_job_id, batch_path, sim_dict)
                            for sim_dict in sims_to_run
                        ]
                        break

                    try:
                        result = subprocess.run(
                            ["squeue", "-j", str(glost_job_id), "-h"],
                            capture_output=True,
                            text=True,
                            timeout=10,
                        )
                        if not result.stdout.strip():
                            logger.error(
                                f"GLOST job {glost_job_id} not in queue but not all results found"
                            )
                            break
                    except (
                        subprocess.TimeoutExpired,
                        subprocess.CalledProcessError,
                    ) as e:
                        logger.warning(
                            f"Failed to check GLOST job {glost_job_id} status: {e}"
                        )

                    logger.info(
                        f"Still waiting for GLOST job {glost_job_id}... (elapsed: {int(time.time() - start_time)}s)"
                    )
                    time.sleep(poll_interval)

                if (time.time() - start_time) >= max_wait_time:
                    logger.error(
                        f"Timeout waiting for GLOST job {glost_job_id} after {max_wait_time}s"
                    )
            else:
                logger.info(
                    "All results found in existing PKL files, no GLOST job needed"
                )

            logger.info(
                f"Job completion summary: {len(existing_results)} from existing files, {len(successful_jobs)} new successful jobs"
            )

            # Start with existing results and collect missing ephys files
            results = []
            all_missing_ephys = []
            failed_processing = []

            # Process existing results first
            for result in existing_results:
                if len(result) == 3:  # New format with missing_ephys
                    protocol_id, epsp_ratio, missing_ephys = result
                    if missing_ephys:
                        all_missing_ephys.extend(missing_ephys)
                    if protocol_id is not None and epsp_ratio is not None:
                        results.append((protocol_id, epsp_ratio))
                else:  # Old format
                    results.append(result)

            # Process newly completed jobs and collect missing ephys
            for _, batch_path, sim_dict in successful_jobs:
                try:
                    result = compute_epsp_ratio_prefire_batch(
                        param_values, sim_dict, sim_dict["workdir"], log_details=False
                    )
                    if len(result) == 3:
                        protocol_id, epsp_ratio, missing_ephys = result
                        if missing_ephys:
                            all_missing_ephys.extend(missing_ephys)
                        if protocol_id is not None and epsp_ratio is not None:
                            results.append((protocol_id, epsp_ratio))
                    else:
                        results.append(result)
                except Exception as e:
                    logger.error(
                        f"Failed to process results for {sim_dict['protocol_id']}: {str(e)}"
                    )
                    failed_processing.append((sim_dict, str(e)))

            if batch_path and os.path.exists(batch_path):
                try:
                    os.remove(batch_path)
                except FileNotFoundError:
                    pass

            # If there are missing ephys files, submit them as a job array
            if all_missing_ephys:
                logger.info(
                    f"Found {len(all_missing_ephys)} missing ephys files, submitting as job array"
                )
                ephys_job_id = self.submit_ephys_job_array(all_missing_ephys)

                if ephys_job_id:
                    # Wait for ephys job array to complete
                    logger.info(
                        f"Waiting for ephys job array {ephys_job_id} to complete..."
                    )
                    max_ephys_wait = 7200  # 2 hours
                    ephys_start_time = time.time()

                    while (time.time() - ephys_start_time) < max_ephys_wait:
                        # Check if all ephys files now exist
                        all_exist = True
                        for ephys_info in all_missing_ephys:
                            if not os.path.exists(ephys_info["file_path"]):
                                all_exist = False
                                break

                        if all_exist:
                            logger.info("All ephys files generated successfully")
                            break

                        # Check if job array is still running
                        result = subprocess.run(
                            ["squeue", "-j", str(ephys_job_id), "-h"],
                            capture_output=True,
                            text=True,
                            timeout=10,
                        )
                        if not result.stdout.strip():
                            # Job not in queue anymore
                            logger.warning(
                                f"Ephys job array {ephys_job_id} completed but some files missing"
                            )
                            break

                        time.sleep(30)  # Check every 30 seconds

                    logger.info(
                        "Retrying result processing with generated ephys files..."
                    )
                    for _, _, sim_dict in successful_jobs:
                        needs_retry = False
                        for result in results:
                            if (
                                len(result) == 2
                                and result[0] == sim_dict["protocol_id"]
                            ):
                                needs_retry = False
                                break
                        else:
                            needs_retry = True

                        if needs_retry:
                            try:
                                result = compute_epsp_ratio_prefire_batch(
                                    param_values,
                                    sim_dict,
                                    sim_dict["workdir"],
                                    log_details=False,
                                )
                                if len(result) == 3:
                                    protocol_id, epsp_ratio, _ = result
                                    if (
                                        protocol_id is not None
                                        and epsp_ratio is not None
                                    ):
                                        results.append((protocol_id, epsp_ratio))
                                else:
                                    results.append(result)
                            except Exception as e:
                                logger.error(
                                    f"Retry failed for {sim_dict['protocol_id']}: {e}"
                                )

            # Process results using shared logic
            error = self._process_evaluation_results(results, param_values)

            # Log failed processing details if any
            if failed_processing:
                logger.error("Failed processing details:")
                for sim_dict, error_msg in failed_processing:
                    logger.error(f"  - {sim_dict['protocol_id']}: {error_msg}")

            logger.info("Evaluation completed successfully")
            return error
        except Exception:
            # Make sure exception and backtrace are thrown back to parent process
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))

    def init_simulator_and_evaluate_with_lists(self, param_values):
        """
        Set NEURON variables and run evaluation with lists.
        Setting the NEURON variables is necessary when using `ipyparallel`,
        since the new subprocesses have pristine NEURON.
        """
        self.sim.initialize()
        return self.evaluate_with_lists(param_values)
