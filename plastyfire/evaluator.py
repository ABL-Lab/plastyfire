"""
Custom BluePyOpt evaluator for the Graupner & Brunel model
authors: Giuseppe Chindemi (12.2020) + minor modifications by András Ecker (06.2024)
"""

import os
import sys
import time
import pickle
import glob
import logging
import traceback
import hashlib
import numpy as np
import pandas as pd
from bluepyopt.evaluators import Evaluator
from bluepyopt.ephys.simulators import NrnSimulator
from bluepyopt.objectives import Objective
from bluepyopt.parameters import Parameter
from itertools import product
import subprocess
from plastyfire.ephysutils import get_epsp_vector
    


from plastyfire.config import OptConfig
from plastyfire.pyslurm import submitjob, canceljob, wait_for_job_slots

MIN2MS = 60 * 1000.
FITTED_TAU = 278.3177658387  # previously optimized time constant of Ca*
# could use `SingletonWeightObjective`s, but it's easier to just multiply the (ordered) errors with the values below...
WEIGHT_REDUCE = np.array([1 / 8] * 2 + [1 / 4] * 3)  # weights of each protocol (lower for the first two)
CONFIGS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/configs"
logger = logging.getLogger(__name__)
DEBUG = False

def compute_epsp_prefire(pkl_file, window=100):
    """Safely compute EPSP values with bounds checking"""
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    
    logger.info(f"Data keys: {data.keys()}")
    t, v, spikes = data['t'], data['v'], data['pre_spikes']
    
    # Extract rho values
    if len(data['rho_GB']) > 100:
        data['rho_GB'] = np.transpose(data['rho_GB'])
    
    initial_rho = [0 if k[0] < 0.5 else 1 for k in data['rho_GB']]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in data['rho_GB']]
    

    valid_spikes = spikes[spikes < (t[-1] - window)][:60]
    
    epsp_values = get_epsp_vector(t, v, valid_spikes, window)  
    avg_epsp = np.mean(epsp_values)

    return avg_epsp


def compute_epsp_ratio_prefire_batch(param_values, sim_dict, workdir, log_details=True):
    """Reads results from batch job output files and computes EPSP ratio"""
    from plastyfire.ephysutils import Experiment
    import pickle
    import os
    import glob
    
    # Reduced logging: only log if there are issues
    logger.info(f"Starting result processing for protocol {sim_dict['protocol_id']} in {workdir}")
    
    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]  # Use first 12 chars
    
    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")
    
    if not os.path.exists(pkl_file):
        logger.error(f"Specific pickle file not found: {pkl_file}")
        raise FileNotFoundError(f"Simulation result file not found: simulation_{param_hash}.pkl")
            
    logger.info(f"Reading results from: {pkl_file}")

    with open(pkl_file, 'rb') as f:
        raw_results = pickle.load(f)
    
    # get initial rho and final rho
    initial_rho = [0 if k[0] < 0.5 else 1 for k in raw_results['rho_GB'].T]
    final_rho = [0 if k[-1] < 0.5 else 1 for k in raw_results['rho_GB'].T]
    initial_rho_str = '_'.join(map(str, initial_rho))
    final_rho_str = '_'.join(map(str, final_rho))
    
    os.remove(pkl_file)

    pre_gid, post_gid = workdir.split('/')[-2].split('-')
    sim_config_path = os.path.join(workdir, "simulation_config.json")    

    ephys_before_filename = f"ephys_data_{pre_gid}_{post_gid}_{initial_rho_str}.pkl"
    ephys_before_path = os.path.join(os.path.dirname(sim_config_path), '..', '..', '..', '..', 'ephys_data', ephys_before_filename)
    
    ephys_after_filename = f"ephys_data_{pre_gid}_{post_gid}_{final_rho_str}.pkl"
    ephys_after_path = os.path.join(os.path.dirname(sim_config_path), '..', '..', '..', '..', 'ephys_data', ephys_after_filename)
    
    for rho_config, file_path in [(initial_rho_str, ephys_before_path), (final_rho_str, ephys_after_path)]:
        if not os.path.exists(file_path):
            logger.info(f"Ephys file not found, generating: {file_path}")
            try:
                # Run compute_test_pulse_epsp_ratio.py
                cmd = [
                    sys.executable,
                    "compute_test_pulse_epsp_ratio.py",
                    sim_config_path,
                    pre_gid,
                    post_gid,
                    rho_config.replace('_', ',')
                ]
                logger.info(f"Running command: {cmd}")
                # Change to the directory containing the script
                script_dir = os.path.dirname(__file__)
                result = subprocess.run(cmd, cwd=script_dir, capture_output=True, text=True, timeout=1800)  # Reduced timeout to 5 minutes
                
                if result.returncode != 0:
                    logger.error(f"Failed to generate ephys file: {result.stderr}")
                    # Return early if ephys file generation fails
                    logger.error(f"Ephys file generation failed, returning None for protocol {sim_dict['protocol_id']}")
                    return None, None
                else:
                    logger.info(f"Successfully completed subprocess")
                    
            except subprocess.TimeoutExpired:
                logger.error(f"Timeout generating ephys file: {file_path}")
                return None, None
            except Exception as e:
                logger.error(f"Error generating ephys file: {e}")
                return None, None
 
    # If files exist, compute expected EPSP ratio
    if os.path.exists(ephys_before_path) and os.path.exists(ephys_after_path):
        epsp_before = compute_epsp_prefire(ephys_before_path)
        epsp_after = compute_epsp_prefire(ephys_after_path)
        epsp_ratio = epsp_after / epsp_before
    else:
        logger.error(f"Ephys files not found: {ephys_before_path} or {ephys_after_path}")
        return None, None

    data = {
        'pkl_file': [pkl_file],
        'protocol_id': [sim_dict['protocol_id']],
        'param_hash': [param_hash],
        'epsp_ratio': [epsp_ratio]
    }
    df = pd.DataFrame(data)
    
    results_file = os.path.join(os.path.dirname(__file__), "simulation_epsp_df.csv")
    if os.path.exists(results_file):
        df.to_csv(results_file, mode='a', header=False, index=False)
    else:
        df.to_csv(results_file, mode='w', header=True, index=False)
    
    logger.info(f"Successfully computed EPSP ratio: {epsp_ratio} for protocol {sim_dict['protocol_id']}")
    return sim_dict["protocol_id"], epsp_ratio

def compute_epsp_ratio_batch(param_values, sim_dict, workdir, log_details=True):
    """Reads results from batch job output files and computes EPSP ratio"""
    from plastyfire.ephysutils import Experiment
    import pickle
    import os
    import glob
    
    # Reduced logging: only log if there are issues
    logger.info(f"Starting result processing for protocol {sim_dict['protocol_id']} in {workdir}")
    
    param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]  # Use first 12 chars
    
    pkl_file = os.path.join(workdir, f"simulation_{param_hash}.pkl")
    
    if not os.path.exists(pkl_file):
        logger.error(f"Specific pickle file not found: {pkl_file}")
        raise FileNotFoundError(f"Simulation result file not found: simulation_{param_hash}.pkl")
            
    logger.info(f"Reading results from: {pkl_file}")

    with open(pkl_file, 'rb') as f:
        raw_results = pickle.load(f)
    
    
    exp_handler = Experiment(raw_results, c01duration=sim_dict["c01duration"],
                            c02duration=sim_dict["c02duration"], period=sim_dict["period"])
    
    # Reduced logging: removed verbose debug
    logger.info(f"Computing EPSP ratio with nepsp={sim_dict['nepsp']}")
    epsp_ratio = exp_handler.compute_epsp_ratio(sim_dict["nepsp"])
    
    # Save data as pandas DataFrame (surgical addition - no return modification)
    data = {
        'pkl_file': [pkl_file],
        'protocol_id': [sim_dict['protocol_id']],
        'param_hash': [param_hash],
        'epsp_ratio': [epsp_ratio]
    }
    df = pd.DataFrame(data)
    
    # Append to a single global results file (thread-safe approach)
    # Use absolute path to ensure consistent location
    results_file = os.path.join(os.path.dirname(__file__), "simulation_epsp_df.csv")
    if os.path.exists(results_file):
        df.to_csv(results_file, mode='a', header=False, index=False)
    else:
        df.to_csv(results_file, mode='w', header=True, index=False)
    
    # Reduced logging: only log essential info
    logger.info(f"Successfully computed EPSP ratio: {epsp_ratio} for protocol {sim_dict['protocol_id']}")
    return sim_dict["protocol_id"], epsp_ratio

class Evaluator(Evaluator):
    """Graupner & Brunel plasticity model evaluator"""
    def __init__(self, fit_params, invitro_db, seed, sample_size, ipp_id, work_dir=None, max_jobs=900):
        super(Evaluator, self).__init__()
        self.params = [Parameter(param_name, bounds=(min_bound, max_bound))
                       for param_name, min_bound, max_bound in fit_params]
        self.param_names = [param.name for param in self.params]
        self.invitro_db = invitro_db
        self.seed = seed
        self.sample_size = sample_size
        self.ipp_id = ipp_id
        self.sim = NrnSimulator()
        self.current_generation = 0  # Track current generation
        self.work_dir = work_dir or os.getcwd()  # Working directory for batch scripts
        self._processed_files = set()  # Track processed files to prevent infinite loops
        self.max_jobs = max_jobs  # Maximum concurrent SLURM jobs
        # Find all simulations
        self.all_sims, self.objectives = [], []
        for elem in self.invitro_db.itertuples():
            # Add objective
            self.objectives.append(Objective(elem.protocol_id))
            # Load simulation config and extract simulation global parameters
            config = OptConfig(os.path.join(CONFIGS_DIR, "%s_%s.yaml" % (elem.pre_mtype, elem.post_mtype)))
            np.testing.assert_almost_equal(config.T / 1000., elem.period_sweep)
            fastforward = config.fastforward
            if fastforward is None:
                fastforward = config.C01_duration * MIN2MS + config.nreps * config.T
            # Ensure fastforward is a clean float to avoid precision issues
            fastforward = float(fastforward)
            nepsp = int(config.C01_duration * MIN2MS / config.T)
            # Load simulation index (witten by `simwriter.py`)
            sim_idx = pd.read_csv(os.path.join(os.path.split(os.path.split(config.out_dir)[0])[0],
                                               "index_%s_%s.csv" % (elem.pre_mtype, elem.post_mtype)))
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

            self.all_sims.extend([{"protocol_id": elem.protocol_id, "period": elem.period_sweep,
                                   "c01duration": config.C01_duration, "c02duration": config.C02_duration,
                                   "fastforward": fastforward, "nepsp": nepsp, "simpath": path, "workdir": "/".join(path.split('/')[:-1])} for path in paths])
            
        logger.info(f"Found {len(self.all_sims)} simulations")
        # Reduced logging: removed verbose debug output
        # logger.debug("Available sims:")
        # for sim in self.all_sims:
        #     logger.debug(sim)

    def set_generation(self, generation):
        """Set current generation for batch job naming"""
        if generation != self.current_generation:
            logger.info(f"Starting Generation {generation} (Population size: {len(self.all_sims) if hasattr(self, 'all_sims') else 'unknown'} simulations per individual) ===")
        self.current_generation = generation

    def get_param_dict(self, param_values):
        """Build dictionary of parameters for the Graupner & Brunel model
        from an ordered list of values (i.e. an individual)"""
        return dict(zip(self.param_names, param_values))

    def create_batch_script(self, param_values, sim_dict):
        """Create batch script for individual evaluation in the protocol directory"""
        param_dict = self.get_param_dict(param_values)
        
        # Create hash from parameter values
        param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]  # Use first 12 chars
        
        # Build parameter arguments string
        param_args = " ".join([f"--{name}={value}" for name, value in param_dict.items()])
        
        # Read template
        template_path = os.path.join(os.path.dirname(__file__), "templates", "simulation.batch.tmpl")
        # Reduced logging: removed verbose debug
        # logger.debug(f"Template path: {template_path}")
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template not found: {template_path}")
        with open(template_path, 'r') as f:
            template = f.read()
        
        # Template variables
        template_vars = {
            "name": f"param_{param_hash}_{sim_dict['protocol_id']}",
            "cpu_time": "01:30:00",
            "log": f"{param_hash}_{sim_dict['protocol_id']}",
            "env": "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/setupenv.sh",
            "run": "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/plastyfire/pairrunner.py",
            "param_args": param_args,
            "fastforward": float(sim_dict["fastforward"]),  
            "workdir": sim_dict["workdir"],  
            "param_hash_arg": f" --param_hash={param_hash}",  
            "individual_id_arg": "", 
            "generation_arg": "" 
        }
        
        # Fill template
        batch_script = template.format(**template_vars)
        
        # Write batch script in the protocol directory with unique name to prevent race conditions
        batch_path = os.path.join(sim_dict["workdir"], f"simulation_{param_hash}.batch")
        # Reduced logging: removed verbose debug
        # logger.debug(f"Creating batch script at: {batch_path}")
        with open(batch_path, 'w') as f:
            f.write(batch_script)
        
        return batch_path

    def submit_batch_job(self, batch_path):
        """Submit batch job and return job ID with job limit enforcement"""
        # Wait for available job slots before submitting
        wait_for_job_slots(max_jobs=self.max_jobs)
        
        # Check if batch script exists
        if not os.path.exists(batch_path):
            logger.error(f"Batch script not found: {batch_path}")
            raise FileNotFoundError(f"Batch script not found: {batch_path}")
        
        # Submit without --parsable flag to avoid potential issues
        result = subprocess.run(["sbatch", batch_path], 
                               capture_output=True, text=True, check=True)
        
        # Parse job ID from output (format: "Submitted batch job 12345")
        output = result.stdout.strip()
        job_id_match = output.split()[-1]  # Get last word which should be the job ID
        try:
            job_id = int(job_id_match)
        except ValueError:
            logger.error(f"Could not parse job ID from output: {output}")
            raise ValueError(f"Could not parse job ID from output: {output}")
            
        logger.info(f"Submitted job {job_id}")
        return job_id
    
    def check_pkl_completion(self, workdir, param_hash):
        """Check if simulation completed by looking for the specific PKL file"""
        import glob
        
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
                            batch_file = os.path.join(workdir, f"simulation_{param_hash}.batch")
                            if os.path.exists(batch_file):
                                os.remove(batch_file)
                            return True
                    except (OSError, IOError):
                        # File might be being written, continue checking
                        continue
            
            return False
            
        except Exception as e:
            return False
    
    def _log_job_failure_details(self, job_id, job_info):
        """Extract and log detailed failure information from job status"""
        try:
            # Parse job info to extract failure details
            lines = job_info.split('|')
            for line in lines:
                if '=' in line:
                    key, value = line.split('=', 1)
                    if key.strip() in ['JobState', 'ExitCode', 'DerivedExitCode', 'Reason', 'Comment']:
                        logger.error(f"Job {job_id} {key.strip()}: {value.strip()}")
            
            # Try to get job history for more details
            hist_result = subprocess.run(["sacct", "-j", str(job_id), "-o", "JobID,JobName,State,ExitCode,DerivedExitCode,Reason,MaxRSS,MaxVMSize,Elapsed", "--noheader"], 
                                        capture_output=True, text=True)
            if hist_result.returncode == 0 and hist_result.stdout.strip():
                logger.error(f"Job {job_id} detailed history:")
                for line in hist_result.stdout.strip().split('\n'):
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
            for filename in os.listdir('.'):
                if filename.startswith(f'slurm-{job_id}'):
                    logger.info(f"Found SLURM output file: {filename}")
                    try:
                        with open(filename, 'r') as f:
                            content = f.read()
                            # Log last few lines that might contain error info
                            lines = content.strip().split('\n')
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
        """Evaluate individual using batch jobs instead of ipyparallel"""
        try:
            # Reduced logging: only log essential evaluations
            logger.info("Evaluating individual: %s", param_values)
            
            # Add a simple check to prevent infinite loops
            if hasattr(self, '_evaluation_count'):
                self._evaluation_count += 1
                if self._evaluation_count > 1000:  # Prevent infinite loops
                    logger.error("Too many evaluations detected, possible infinite loop. Exiting.")
                    return [1000.0] * len(self.objectives)
            else:
                self._evaluation_count = 1
                
            cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
            pklf_name = os.path.join(".cache", "%s.pkl" % cachekey)
            if os.path.isfile(pklf_name):
                with open(pklf_name, "rb") as f:
                    cache_data = pickle.load(f)  
                np.testing.assert_array_equal(param_values, cache_data["individual"])  # verify no collision (OMG)
                logger.info(f"Returning cached results for individual {cachekey[:12]}")
                
                param_hash = cachekey[:12]
                if "resdb" in cache_data:
                    res_db = cache_data["resdb"]
                    cached_data = []
                    for _, row in res_db.iterrows():
                        cached_data.append({
                            'pkl_file': f"cached_{param_hash}_{row['protocol_id']}.pkl",
                            'protocol_id': row['protocol_id'],
                            'param_hash': param_hash,
                            'epsp_ratio': row['epsp_ratio']
                        })
                    
                    if cached_data:
                        cached_df = pd.DataFrame(cached_data)
                        results_file = os.path.join(os.path.dirname(__file__), "simulation_epsp_df.csv")
                        if os.path.exists(results_file):
                            cached_df.to_csv(results_file, mode='a', header=False, index=False)
                        else:
                            cached_df.to_csv(results_file, mode='w', header=True, index=False)
                
                return cache_data["error"]  
            
            param_dict = self.get_param_dict(param_values)
            param_hash = hashlib.md5(str(param_values).encode()).hexdigest()[:12]  # Use first 12 chars
            batch_jobs = []
            existing_results = []
            
            for i, sim_dict in enumerate(self.all_sims):
                pkl_file = os.path.join(sim_dict["workdir"], f"simulation_{param_hash}.pkl")
                
                file_key = f"{pkl_file}_{param_hash}"
                if file_key in self._processed_files:
                    logger.warning(f"Already processed file {pkl_file} in this evaluation, skipping to prevent infinite loop")
                    continue
                
                if os.path.exists(pkl_file):
                    logger.info(f"Found existing PKL file for {sim_dict['protocol_id']}: {pkl_file}")
                    self._processed_files.add(file_key)  # Mark as processed
                    try:
                        result = compute_epsp_ratio_prefire_batch(param_values, sim_dict, sim_dict["workdir"], log_details=False)
                        if result is not None and result[1] is not None:  # Check for valid result
                            existing_results.append(result)
                            continue
                        else:
                            logger.warning(f"Invalid result from existing PKL file {pkl_file}, will submit new batch job")
                    except Exception as e:
                        logger.info(f"Failed to load existing PKL file {pkl_file}: {e}")
                else:
                    logger.info(f"No existing PKL file found for {sim_dict['protocol_id']}: {pkl_file}, submitting batch job")
                
                batch_path = self.create_batch_script(param_values, sim_dict)
                job_id = self.submit_batch_job(batch_path)
                batch_jobs.append((job_id, batch_path, sim_dict, param_hash))
            
            failed_jobs = []
            successful_jobs = []
            
            if batch_jobs:
                logger.info(f"Waiting for {len(batch_jobs)} batch jobs to complete...")
                for job_id, batch_path, sim_dict, param_hash in batch_jobs:
                    
                    if self.check_pkl_completion(sim_dict["workdir"], param_hash):
                        logger.info(f"Job {job_id} completed successfully")
                        successful_jobs.append((job_id, batch_path, sim_dict))
                    else:
                        logger.error(f"Job {job_id} failed")
                        failed_jobs.append((job_id, batch_path, sim_dict))
            else:
                logger.info("All results found in existing PKL files, no batch jobs needed")
            
            # Log summary of job results
            logger.info(f"Job completion summary: {len(existing_results)} from existing files, {len(successful_jobs)} new successful jobs, {len(failed_jobs)} failed jobs")
                        
            # Start with existing results
            results = existing_results.copy()
            failed_processing = []
            
            # Add results from newly completed jobs
            for job_id, batch_path, sim_dict in successful_jobs:
                try:
                    # Read results from job output files
                    result = compute_epsp_ratio_prefire_batch(param_dict, sim_dict, os.path.dirname(sim_dict["simpath"]))
                    results.append(result)
                    # Reduced logging: removed verbose debug
                    # logger.debug(f"Successfully processed results for {sim_dict['protocol_id']}: {result}")
                except Exception as e:
                    logger.error(f"Failed to process results for {sim_dict['protocol_id']}: {str(e)}")
                    failed_processing.append((sim_dict, str(e)))
                
                # Clean up batch script (now with unique naming)
                if os.path.exists(batch_path):
                    os.remove(batch_path)
                [os.remove(f) for f in glob.glob(os.path.join(os.path.dirname(batch_path), "*.log"))]
            
            # Check if we have enough results to proceed
            if len(results) < len(self.objectives):
                logger.error(f"Insufficient results: got {len(results)}, need {len(self.objectives)}")
                logger.error("Failed processing details:")
                for sim_dict, error in failed_processing:
                    logger.error(f"  - {sim_dict['protocol_id']}: {error}")
                
                # Return high error values for insufficient results
                high_error = [1000.0] * len(self.objectives)
                logger.warning(f"Returning high error values due to insufficient results: {high_error}")
                return high_error
            
            # Assemble results
            res_db = pd.DataFrame(results, columns=["protocol_id", "epsp_ratio"])

            logger.info("Simulations completed, results:")
            logger.info(res_db)
            insilico_db = res_db.groupby("protocol_id")["epsp_ratio"].agg(["mean", "sem"]).add_suffix("_epsp_ratio").reset_index()

            logger.info("Aggregating results")
            logger.info(insilico_db)
            # Compute error, ensuring order
            merged_db = pd.merge(self.invitro_db, insilico_db, on="protocol_id", suffixes=("_invitro", "_insilico"))

            logger.info("Joining in vitro and in silico results")
            logger.info(merged_db)
            merged_db["error"] = np.abs((merged_db["mean_epsp_ratio_invitro"] - merged_db["mean_epsp_ratio_insilico"])
                                        / merged_db["sem_epsp_ratio_invitro"])
            error = [float(merged_db.loc[merged_db["protocol_id"] == obj.name, "error"].iloc[0])
                     for obj in self.objectives]
            error = (error * WEIGHT_REDUCE).tolist()  # weight protocols
            # Reduced logging: removed verbose debug output
            logger.info("Sorting errors")
            logger.info(error)
            outcome = [float(merged_db.loc[merged_db["protocol_id"] == obj.name, "mean_epsp_ratio_insilico"].iloc[0])
                       for obj in self.objectives]
            # Store results in cache
            cachekey = hashlib.md5(str(param_values).encode()).hexdigest()
            pklf_name = os.path.join(".cache", "%s.pkl" % cachekey)
            with open(pklf_name, "wb") as f:
                pickle.dump({"error": error, "outcome": outcome, "individual": list(param_values), "resdb": res_db},
                            f, -1)
            # Reduced logging: removed verbose debug
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
