"""
Small helper functions to submit jobs to Slurm from Python
authors: Giuseppe Chindemi (12.2020) + minor modifications and docs by András Ecker (05.2024)
"""


import re
import time
import random
import logging
import subprocess

logger = logging.getLogger(__name__)


def get_active_job_count():
    """Get the current number of active SLURM jobs for the user"""
    import os
    try:
        username = os.environ.get('USER', os.environ.get('USERNAME', 'unknown'))
        result = subprocess.run(["squeue", "-u", username, "-h"], 
                               capture_output=True, text=True, check=True)
        lines = [line for line in result.stdout.strip().split('\n') if line.strip()]
        job_count = len(lines) if result.stdout.strip() else 0
        return job_count
    except subprocess.CalledProcessError as e:
        logger.warning(f"Failed to get job count: {e}")
        return 0


def wait_for_job_slots(max_jobs=900, check_interval=30):
    """Wait until there are available job slots"""
    while True:
        current_jobs = get_active_job_count()
        logger.info(f"Current active jobs: {current_jobs}/{max_jobs}")
        
        if current_jobs < max_jobs:
            break
            
        logger.info(f"Job limit reached ({current_jobs}/{max_jobs}), waiting {check_interval}s...")
        time.sleep(check_interval)


def submitjob(block=True):
    """Submit job and get its ID"""
    time.sleep(random.randint(0, 60))
    job_out = subprocess.run(["sbatch", "--parsable", "ipp.sh"], check=True, text=True, stdout=subprocess.PIPE).stdout
    jobid = int(re.findall(r'\d+', job_out)[0])
    logger.debug("Submitted job %i", jobid)
    if block:
        # Wait for job to start running
        while True:
            logger.debug("Waiting for job to start")
            job_info = subprocess.run(["scontrol", "show", "job", "%i" % jobid, "-o"],
                                      check=True, text=True, stdout=subprocess.PIPE).stdout
            status = re.findall(r"JobState=(\w*)", job_info)[0]
            logger.debug("Job status = %s", status)
            if status == "RUNNING":
                break
            else:  # try again in a bit
                time.sleep(10)
    return jobid


def canceljob(jobid):
    subprocess.run(["scancel", "%i" % jobid], check=True, text=True, stdout=subprocess.PIPE).stdout