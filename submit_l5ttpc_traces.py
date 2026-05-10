#!/usr/bin/env python3
"""
Submit L5TTPC_L5TTPC plasticity simulations with Chindemi parameters to SLURM,
recording synaptic traces (effcai_GB, cai_CR, vsyn, rho_GB) during induction.

Usage:
    # Submit all 100 pairs as SLURM jobs (overnight run)
    python submit_l5ttpc_traces.py
    
    # Test with 2 pairs first
    python submit_l5ttpc_traces.py --max-pairs 2
    
    # Use multiprocessing instead of SLURM
    python submit_l5ttpc_traces.py --execution-mode cpu --workers 30

Author: Generated for trace recording
"""

import argparse
import logging
import os
import pickle
import subprocess
import sys
from pathlib import Path

from plastyfire.epg_dhuruva_custom_ratio import RHO_RATIO_ENV, parse_rho_ratio

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Base directories
SIMULATIONS_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations"
SIMULATIONS_DIR_STDP = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC_STDP/simulations"
OUTPUT_BASE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/full_trace_results"
PLASTYFIRE_DIR = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire"
DEFAULT_EPG_FULL_PAIRS_FILE = os.path.join(PLASTYFIRE_DIR, "data", "pairs_n100.txt")

# Chindemi parameters (from evaluate_best_solution.py)
CHINDEMI_PARAMS = {
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
    "tau_effca_GB_GluSynapse": 278.318,
}

DHURUVA_PARAMS = {
    "tau_effca_GB_GluSynapse": 111.4966,
    "gamma_d_GB_GluSynapse": 88.70172982129458,
    "gamma_p_GB_GluSynapse": 268.32068636660046,
    "a00": 0.4777981934076907,
    "a01": 0.6934394338463974,
    "a10": 1.0586973073436239,
    "a11": 3.465063729118567,
    "a20": 0.7322653809967383,
    "a21": 0.14746008318800907,
    "a30": 1.88136100637871,
    "a31": 4.447189652961518,
}

DHURUVA_PARAMS_V15 = {
    "tau_effca_GB_GluSynapse": 292.0600,
    "gamma_d_GB_GluSynapse": 418.2236,
    "gamma_p_GB_GluSynapse": 150.1040,
    "a00": 1.0021,   # basal depression pre
    "a01": 2.2192,   # basal depression post
    "a10": 3.2803,   # basal potentiation pre
    "a11": 1.7763,   # basal potentiation post
    "a20": 1.0021,   # apical depression pre  (= a00, shared)
    "a21": 2.2192,   # apical depression post (= a01, shared)
    "a30": 3.2803,   # apical potentiation pre  (= a10, shared)
    "a31": 1.7763,   # apical potentiation post (= a11, shared)
}

DHURUVA_PARAMS_V2 = {
    "tau_effca_GB_GluSynapse": 73.21360951908201,
    "gamma_d_GB_GluSynapse": 109.4093089925484,
    "gamma_p_GB_GluSynapse": 596.8444303473575,
    "a00": 0.2933119882453923,
    "a01": 1.5969307382598152,
    "a10": 0.8408881959013523,
    "a11": 4.8170123961535785,
    "a20": 0.6607260497664416,
    "a21": 0.010296203791785974,
    "a30": 3.161785481914619,
    "a31": 3.6158549749054814,
    "tau_prime_CICR_GluSynapse": 1109.1151367901289,
    "k_prime_CICR_GluSynapse": 17516.914653528733,
    "K_prime_limit_CICR_GluSynapse": 0.00046795676802132444,
    "K_RyR_CICR_GluSynapse": 0.0020996724134375957,
    "n_RyR_CICR_GluSynapse": 4.457570893183441,
    "Vmax_CICR_GluSynapse": 0.0005044554223042098,
    "tau_rel_CICR_GluSynapse": 223.5880115766854,
}

DHURUVA_PARAMS_V3 = {
    "gamma_d_GB_GluSynapse": 50.4690,
    "gamma_p_GB_GluSynapse": 161.2494,
    "a00": 0.5169,
    "a01": 4.6079,
    "a10": 1.8423,
    "a11": 4.9230,
    "a20": 0.1282,
    "a21": 0.4070,
    "a30": 0.3265,
    "a31": 4.9275,
    "tau_prime_CICR_GluSynapse": 436.0976,
    "k_prime_CICR_GluSynapse": 44310.3178,
    "K_prime_limit_CICR_GluSynapse": 0.0003,
    "K_RyR_CICR_GluSynapse": 0.0013,
    "n_RyR_CICR_GluSynapse": 5.0744,
    "Vmax_CICR_GluSynapse": 0.0022,
    "tau_rel_CICR_GluSynapse": 313.2375,
    "tau_effca_GB_GluSynapse": 100.7488,
}

DHURUVA_PARAMS_V4 = {
    "gamma_d_GB_GluSynapse": 171.5880,
    "gamma_p_GB_GluSynapse": 598.7721,
    "a00": 2.3404,
    "a01": 2.2347,
    "a10": 2.7113,
    "a11": 4.0093,
    "a20": 2.6263,
    "a21": 2.6543,
    "a30": 4.4278,
    "a31": 3.9129,
    "tau_prime_CICR_GluSynapse": 1989.2330,
    "k_prime_CICR_GluSynapse": 41683.2555,
    "K_prime_limit_CICR_GluSynapse": 0.0010,
    "K_RyR_CICR_GluSynapse": 0.0049,
    "n_RyR_CICR_GluSynapse": 4.5044,
    "Vmax_CICR_GluSynapse": 0.0019,
    "tau_rel_CICR_GluSynapse": 450.8465,
    "tau_effca_GB_GluSynapse": 344.6220,
}

DHURUVA_PARAMS_V5 = {
    "gamma_d_GB_GluSynapse": 69.35582333947616,
    "gamma_p_GB_GluSynapse": 172.8162926754338,
    "a00": 1.0338281805504534,
    "a01": 2.0824027912480663,
    "a10": 1.50981772882111,
    "a11": 4.937595300732932,
    "a20": 1.024923299559302,
    "a21": 1.0040979194963504,
    "a30": 4.869465155545409,
    "a31": 4.37081939321435,
    "tau_prime_CICR_GluSynapse": 870.7718136117252,
    "k_prime_CICR_GluSynapse": 17653.30973320602,
    "K_prime_limit_CICR_GluSynapse": 0.0001738667305935812,
    "K_RyR_CICR_GluSynapse": 0.002333447308131663,
    "n_RyR_CICR_GluSynapse": 4.637471212566512,
    "Vmax_CICR_GluSynapse": 0.0009454229634058555,
    "tau_rel_CICR_GluSynapse": 71.00509272517166,
    "tau_effca_GB_GluSynapse": 279.72967330670116,
}

DHURUVA_PARAMS_V6 = {
    "gamma_d_GB_GluSynapse": 88.37748728897184,
    "gamma_p_GB_GluSynapse": 75.21696897014064,
    "a00": 1.0016365565804994,
    "a01": 1.4995212257656416,
    "a10": 1.1327719903735114,
    "a11": 3.077770462888927,
    "a20": 1.370297324550421,
    "a21": 2.4040852180572143,
    "a30": 1.9262572961622118,
    "a31": 3.508641458270441,
    # CICR params — confirmed HOC global names (all end in _GluSynapse in compiled .so)
    "k_fill_CICR_GluSynapse": 1.6465937576950695,
    "tau_leak_CICR_GluSynapse": 96157.34999520841,
    "K_clip_CICR_GluSynapse": 0.0008326648885177493,
    "K_trig_CICR_GluSynapse": 0.0019013479472819725,
    "n_trig_CICR_GluSynapse": 5.737730388740266,
    "Vmax_CICR_GluSynapse": 1.6609943109456529,
    "tau_CICR_rel_GluSynapse": 807.4956732356754,
    "g_cicr_GB_GluSynapse": 47234.51951479494,
    "tau_effca_GB_GluSynapse": 484.9361541639738,
}

DHURUVA_PARAMS_V7 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 69.61242611028226,
    "gamma_p_GB_GluSynapse": 166.96578019613747,
    "a00": 1.0003490121729528,
    "a01": 1.0877438002708641,
    "a10": 1.3917207001380245,
    "a11": 4.129505808304705,
    "a20": 2.3079274077684175,
    "a21": 1.9619384872398133,
    "a30": 4.456560780131431,
    "a31": 4.80485148600644,
    "delta_IP3_CICR_GluSynapse": 2.3564846028724626,
    "tau_IP3_CICR_GluSynapse": 4198.8339980611045,
    "phi_serca_CICR_GluSynapse": 0.00010926586342370828,
    "k_leak_CICR_GluSynapse": 0.00038983484712417705,
    "V_RyR_CICR_GluSynapse": 0.0008831061403601626,
    "K_T_CICR_GluSynapse": 0.0016997377463884833,
    "n_T_CICR_GluSynapse": 9.677747356489755,
    "V_IP3R_CICR_GluSynapse": 0.009883023779226525,
    "g_RyR_CICR_GluSynapse": 14060.428640146583,
    "tau_effca_GB_GluSynapse": 306.6936809019515,
}

DHURUVA_PARAMS_V8 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 50.68462528598515,
    "gamma_p_GB_GluSynapse": 269.16205800226743,
    "a00": 1.0113045779827383,
    "a01": 1.2882604966456566,
    "a10": 2.3576812422737237,
    "a11": 4.952090470632914,
    "a20": 1.0812138620096716,
    "a21": 1.2187438504221766,
    "a30": 3.7963239035846366,
    "a31": 4.702457929298914,
    "delta_IP3_CICR_GluSynapse": 6.591298429018293,
    "tau_IP3_CICR_GluSynapse": 4052.7399183766393,
    "phi_serca_CICR_GluSynapse": 1.1929258197945853e-05,
    "k_leak_CICR_GluSynapse": 0.00012377018583732257,
    "V_RyR_CICR_GluSynapse": 0.0018239654504870481,
    "K_T_CICR_GluSynapse": 0.0015936426856495703,
    "n_T_CICR_GluSynapse": 9.889807475847107,
    "V_IP3R_CICR_GluSynapse": 0.00021326149476798612,
    "g_RyR_CICR_GluSynapse": 73.35820407442986,
    "g_IP3R_CICR_GluSynapse": 2093.211205675536,
    "tau_IP3R_CICR_GluSynapse": 20.371885106584717,
    "tau_effca_GB_GluSynapse": 50.00009455922666,
}

DHURUVA_PARAMS_V10 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 249.9229,
    "gamma_p_GB_GluSynapse": 191.8133,
    "a00": 2.8440, "a01": 1.2974, "a10": 4.9952, "a11": 3.5837,
    "a20": 1.0138, "a21": 1.0256, "a30": 1.0804, "a31": 3.3058,
    "delta_IP3_CICR_GluSynapse": 7.1016,
    "tau_IP3_CICR_GluSynapse": 553.6846,
    "tau_effca_GB_GluSynapse": 129.1823,
    # Li-Rinzel CICR params (d1, d2, d3, d5, a2, V_CICR, V_SERCA, K_SERCA,
    # V_leak, gamma_ER, tau_extrusion) are baked into the .mod file defaults
}

DHURUVA_PARAMS_V11 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 244.2704559170318,
    "gamma_p_GB_GluSynapse": 201.88922143839974,
    "a00": 1.3847731157411707, "a01": 2.649582088825211,
    "a10": 4.068764718250964, "a11": 4.2038365115376095,
    "a20": 3.9703813952674856, "a21": 1.4653156595015602,
    "a30": 3.006000244271462, "a31": 2.4354159698117264,
    "delta_IP3_CICR_GluSynapse": 2.9198673346235715,
    "tau_IP3_CICR_GluSynapse": 331.6473628827375,
    "tau_effca_GB_GluSynapse": 131.02615799535653,
    # Li-Rinzel CICR params (V_CICR, V_SERCA, tau_extrusion, d1, d2, d3, d5, a2,
    # K_SERCA, V_leak, gamma_ER) are baked into the .mod file defaults
}

DHURUVA_PARAMS_V12 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 82.19974298716522,
    "gamma_p_GB_GluSynapse": 241.7999404447529,
    "a00": 1.6583210931204957,
    "a01": 1.5136531069201258,
    "a10": 3.809033325093047,
    "a11": 4.5762615641623,
    "a20": 3.352173921151783,
    "a21": 3.3765942854214686,
    "a30": 1.0011444205701536,
    "a31": 1.2569437734232813,
    "delta_IP3_CICR_GluSynapse": 0.001386261047598346,
    "tau_IP3_CICR_GluSynapse": 304.9378710595687,
    "V_IP3R_CICR_GluSynapse": 0.027532028218439433,
    "V_RyR_CICR_GluSynapse": 0.002743740615742242,
    "V_SERCA_CICR_GluSynapse": 0.01329259228286096,
    "K_SERCA_CICR_GluSynapse": 0.006286133934321784,
    "V_leak_CICR_GluSynapse": 0.00014020222355700998,
    "tau_extrusion_CICR_GluSynapse": 34.209295623307625,
    "tau_effca_GB_GluSynapse": 170.14893667200263,
    "tau_ref_CICR_GluSynapse": 796.6915676506252,
    "K_h_ref_CICR_GluSynapse": 0.0007856151840756975,
}

DHURUVA_PARAMS_V9 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 67.1809,
    "gamma_p_GB_GluSynapse": 292.5239,
    "a00": 1.1844,
    "a01": 1.2404,
    "a10": 4.3000,
    "a11": 4.9448,
    "a20": 1.0464,
    "a21": 1.0548,
    "a30": 9.9931,
    "a31": 4.9954,
    "delta_IP3_CICR_GluSynapse": 1.6698,
    "tau_IP3_CICR_GluSynapse": 500.0067,
    "phi_serca_CICR_GluSynapse": 0.0003,
    "k_leak_CICR_GluSynapse": 0.00001,
    "T_threshold_CICR_GluSynapse": 0.3,
    "K_T_CICR_GluSynapse": 0.0024,
    "n_T_CICR_GluSynapse": 6.7190,
    "V_IP3R_CICR_GluSynapse": 0.0001,
    "g_RyR_CICR_GluSynapse": 24917.8768,
    "g_IP3R_CICR_GluSynapse": 2617.9421,
    "tau_RyR_CICR_GluSynapse": 36.0469,
    "tau_IP3R_CICR_GluSynapse": 150.6526,
    "tau_effca_GB_GluSynapse": 50.0160,
}

DHURUVA_PARAMS_V13 = {
    "enable_CICR_GluSynapse": 1,
    "gamma_d_GB_GluSynapse": 193.7068,
    "gamma_p_GB_GluSynapse": 198.2918,
    "a00": 1.4456,
    "a01": 4.9915,
    "a10": 1.0893,
    "a11": 3.4990,
    "a20": 4.9616,
    "a21": 1.7164,
    "a30": 3.7713,
    "a31": 4.4122,
    "delta_IP3_CICR_GluSynapse": 0.0005,
    "tau_IP3_CICR_GluSynapse": 200.8450,
    "V_IP3R_CICR_GluSynapse": 0.0626,
    "V_RyR_CICR_GluSynapse": 0.0044,
    "V_SERCA_CICR_GluSynapse": 0.0024,
    "K_SERCA_CICR_GluSynapse": 0.0057,
    "V_leak_CICR_GluSynapse": 0.0001,
    "tau_extrusion_CICR_GluSynapse": 99.7981,
    "tau_ref_CICR_GluSynapse": 322.4850,
    "K_h_ref_CICR_GluSynapse": 0.0004,
    "tau_effca_GB_GluSynapse": 278.318,
}

DHURUVA_PARAMS_V14 = {
    "gamma_d_GB_GluSynapse": 51.04911994819324,
    "gamma_p_GB_GluSynapse": 499.9990231629091,
    "a00": 1.0710194674877798,
    "a01": 3.127433724876557,
    "a10": 2.619606402608639,
    "a11": 3.534367333843786,
    "a20": 1.3090207037881711,
    "a21": 3.256342887636098,
    "a30": 1.8657221424019035,
    "a31": 4.526034613651234,
    "tau_effca_GB_GluSynapse": 52.01527226929264,
}

DHURUVA_PARAMS_V16 = {
    "gamma_d_GB_GluSynapse": 394.0494,
    "gamma_p_GB_GluSynapse": 974.1213,
    "a00": 7.1333,
    "a01": 1.4068,
    "a10": 5.4961,
    "a11": 1.8273,
    "a20": 2.1500,
    "a21": 1.1701,
    "a30": 1.6150,
    "a31": 1.7262,
    "tau_effca_GB_GluSynapse": 342.2716,
}

DHURUVA_PARAMS_V17 = {
    "gamma_d_GB_GluSynapse": 425.3030566335992,
    "gamma_p_GB_GluSynapse": 191.29497207693194,
    "a00": 1.0587667606965314,   # basal depression pre  (d_v0)
    "a01": 2.3146754819060655,   # basal depression post (d_v1)
    "a10": 3.5703223084695184,   # basal potentiation pre  (p_v0)
    "a11": 2.0375842834692235,   # basal potentiation post (p_v1)
    "a20": 1.9417469661205238,   # apical depression pre  (d_v0a)
    "a21": 2.318822399626272,    # apical depression post (d_v1a)
    "a30": 3.5257610110599797,   # apical potentiation pre  (p_v0a)
    "a31": 2.058889163377155,    # apical potentiation post (p_v1a)
    "tau_effca_GB_GluSynapse": 294.66864436137905,
}

DHURUVA_PARAMS_V18 = {
    "gamma_d_GB_GluSynapse": 421.3925785087944,
    "gamma_p_GB_GluSynapse": 191.8488624639991,
    "a00": 1.0900847272701293,   # basal depression pre  (d_v0)
    "a01": 2.1453182247355085,   # basal depression post (d_v1)
    "a10": 3.3453046820253,      # basal potentiation pre  (p_v0)
    "a11": 1.60351444535477,     # basal potentiation post (p_v1)
    "a20": 2.2868631133386157,   # apical depression pre  (d_v0a)
    "a21": 2.117828378107206,    # apical depression post (d_v1a)
    "a30": 3.5828013457800916,   # apical potentiation pre  (p_v0a)
    "a31": 1.5796432469516355,   # apical potentiation post (p_v1a)
    "tau_effca_GB_GluSynapse": 269.989959751138,
}

DHURUVA_PARAMS_V19 = {
    "gamma_d_GB_GluSynapse":   379.5256059628456,
    "gamma_p_GB_GluSynapse":   162.02006326916813,
    "a00": 1.0770144130928545,   # d_v0  — basal depression pre
    "a01": 2.0379105635626904,   # d_v1  — basal depression post
    "a10": 2.894174205875408,    # p_v0  — basal potentiation pre
    "a11": 1.3571983004135566,   # p_v1  — basal potentiation post
    "a20": 2.3097860482130477,   # d_v0a — apical depression pre
    "a21": 2.317590020162771,    # d_v1a — apical depression post
    "a30": 3.692021341986094,    # p_v0a — apical potentiation pre
    "a31": 1.927119096974426,    # p_v1a — apical potentiation post
    "tau_effca_GB_GluSynapse":  257.9025228124917,
}

def find_jobs_from_prefire_configs(simulations_path):
    """Discover (pair_dir, freq_dt) jobs by finding prefire_simulation_config.json files.

    Each config lives at:  <simulations_path>/<pair_dir>/<freq_dt>/prefire_simulation_config.json
    The pair_dir and freq_dt are derived directly from the file path — no file I/O needed.

    Returns:
        List of (pair_dir_str, freq_dt_str) tuples, sorted for reproducibility.
    """
    jobs = []
    simulations_path = Path(simulations_path)
    config_files = sorted(simulations_path.rglob("simulation_config.json"))
    logger.info(f"Found {len(config_files)} simulation_config.json files under {simulations_path}")
    for cfg_file in config_files:
        freq_dt_dir = cfg_file.parent          # e.g. .../10Hz_-10ms
        pair_dir    = freq_dt_dir.parent       # e.g. .../180164-197248
        # Only include if pair_dir is a direct child of simulations_path
        if pair_dir.parent != simulations_path:
            logger.warning(f"Unexpected nesting depth for {cfg_file}, skipping")
            continue
        jobs.append((str(pair_dir), freq_dt_dir.name))
    return jobs



def submit_slurm_jobs(jobs, output_dir, timeout_hours=12, mem_gb=8, cpus_per_task=10,
                       slurm_account="ctb-emuller", batch_size=10, parallel_sims=4,
                       fastforward=None, epg_variant="epg_dhuruva",
                       epg_full_pairs_file=DEFAULT_EPG_FULL_PAIRS_FILE, rho_ratio=None):
    """Submit jobs to SLURM as batch jobs.
    
    Args:
        jobs: List of (pair_dir, freq_dt_dir) tuples
        output_dir: Output directory for results
        timeout_hours: Job timeout in hours (default: 12 for overnight)
        mem_gb: Memory per job in GB
        cpus_per_task: CPUs per job
        slurm_account: SLURM account
        batch_size: Number of simulations per SLURM job
        fastforward: Fastforward begin point in ms (passed to pairrunner.py), or None
    """
    log_folder = Path(output_dir) / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)
    
    # Create batches
    batches = []
    for i in range(0, len(jobs), batch_size):
        batch = jobs[i:i + batch_size]
        batches.append(batch)
    
    logger.info(f"Submitting {len(batches)} SLURM jobs for {len(jobs)} simulations...")
    
    # Time format
    time_str = f"02:00:00"
    
    job_ids = []
    
    for batch_idx, batch in enumerate(batches):
        job_name = f"trace_{batch_idx:04d}"
        
        # Build the list of workdirs for this batch
        workdirs = [os.path.join(pair_dir, freq_dt) for pair_dir, freq_dt in batch]
        workdirs_str = " ".join(f'"{w}"' for w in workdirs)
        
        # Build parameter args string; optionally append --fastforward
        param_args = " ".join(f"--{name}={value}" for name, value in DHURUVA_PARAMS.items())
        if fastforward is not None:
            param_args += f" --fastforward={fastforward}"
        param_args += f" --epg-variant={epg_variant}"
        if rho_ratio is not None:
            param_args += f' --rho-ratio="{rho_ratio}"'
        
        # Create sbatch script
        sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={slurm_account}
#SBATCH --time={time_str}
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output={log_folder}/{job_name}_%j.out
#SBATCH --error={log_folder}/{job_name}_%j.err

source {PLASTYFIRE_DIR}/setupenv.sh
export PLASTYFIRE_EPG_FULL_PAIRS_FILE="{epg_full_pairs_file}"
{f'export {RHO_RATIO_ENV}="{rho_ratio}"' if rho_ratio is not None else f'unset {RHO_RATIO_ENV}'}

run_sim() {{
    local workdir=$1
    echo "Processing: $workdir"
    
    local pair_name=$(basename $(dirname "$workdir"))
    local freq_dt=$(basename "$workdir")
    local output_file="{output_dir}/$pair_name/$freq_dt/simulation_traces.pkl"
    
    # Skip if already exists
    if [ -f "$output_file" ]; then
        echo "Skipping $pair_name/$freq_dt - already exists"
        return 0
    fi
    
    mkdir -p "{output_dir}/$pair_name/$freq_dt"
    cd "$workdir"
    {sys.executable} {PLASTYFIRE_DIR}/plastyfire/pairrunner.py {param_args} --output-filename=simulation_traces_temp.pkl
    
    if [ -f "$workdir/simulation_traces_temp.pkl" ]; then
        mv "$workdir/simulation_traces_temp.pkl" "$output_file"
        echo "Completed: $pair_name/$freq_dt"
    else
        echo "Failed: $pair_name/$freq_dt"
    fi
}}

# Run each simulation in this batch
for workdir in {workdirs_str}; do
    run_sim "$workdir" &
    
    # Limit number of parallel jobs
    while [ $(jobs -r -p | wc -l) -ge {parallel_sims} ]; do
        wait -n
    done
done

wait
echo "Batch {batch_idx} completed"
"""
        
        # Submit via stdin
        result = subprocess.run(
            ["sbatch"],
            input=sbatch_script,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            logger.info(f"Submitted batch {batch_idx}: job {job_id} ({len(batch)} sims)")
        else:
            logger.error(f"Failed to submit batch {batch_idx}: {result.stderr}")
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Submitted {len(job_ids)} SLURM jobs")
    logger.info(f"Logs: {log_folder}/")
    logger.info(f"Results: {output_dir}/")
    logger.info(f"\nMonitor with: squeue -u $USER")
    
    return job_ids


def _cpu_worker(args):
    """Worker function for CPU multiprocessing mode (must be at module level for pickling).

    Calls `plastyfire.simulator.runconnectedpair_induction` directly (no subprocess),
    which internally uses `_runconnectedpair_process`.

    NOTE: `fit_params` and `fastforward` must be passed explicitly in the tuple — do NOT
    read module-level globals here, because multiprocessing spawns a fresh interpreter
    where reassignments made in main() are NOT visible.
    """
    pair_dir, freq_dt, output_dir, fit_params, fastforward, epg_variant, epg_full_pairs_file, rho_ratio = args
    workdir = os.path.join(pair_dir, freq_dt)
    pair_name = os.path.basename(pair_dir)
    output_file = os.path.join(output_dir, pair_name, freq_dt, "simulation_traces.pkl")

    if os.path.exists(output_file):
        return {"status": "skipped", "pair": pair_name}

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    try:
        import plastyfire.simulator as simulator
        os.environ["PLASTYFIRE_EPG_FULL_PAIRS_FILE"] = epg_full_pairs_file
        if rho_ratio is not None:
            os.environ[RHO_RATIO_ENV] = rho_ratio
        else:
            os.environ.pop(RHO_RATIO_ENV, None)
        results = simulator.runconnectedpair_induction(
            workdir, fit_params=fit_params, fastforward=fastforward, epg_variant=epg_variant)
        with open(output_file, "wb") as f:
            pickle.dump(results, f)
        return {"status": "success", "pair": pair_name, "freq_dt": freq_dt}
    except Exception as e:
        return {"status": "error", "pair": pair_name, "freq_dt": freq_dt, "error": str(e)}


def run_cpu_mode(jobs, output_dir, fit_params, workers=30, fastforward=None,
                 epg_variant="epg_dhuruva", epg_full_pairs_file=DEFAULT_EPG_FULL_PAIRS_FILE,
                 rho_ratio=None):
    """Run simulations locally using multiprocessing.

    `fit_params` and `fastforward` are passed explicitly so that worker
    subprocesses receive the correct values (module-level globals are not
    inherited by workers).
    """
    from multiprocessing import Pool
    import time

    logger.info(f"Running {len(jobs)} simulations with {workers} workers...")

    # Pack fit_params and fastforward into each job tuple so the worker receives them directly
    jobs_with_output = [
        (pair_dir, freq_dt, output_dir, fit_params, fastforward, epg_variant, epg_full_pairs_file, rho_ratio)
        for pair_dir, freq_dt in jobs
    ]

    start = time.time()
    with Pool(workers) as pool:
        results = pool.map(_cpu_worker, jobs_with_output)
    elapsed = time.time() - start

    success = sum(1 for r in results if r["status"] == "success")
    skipped = sum(1 for r in results if r["status"] == "skipped")
    failed  = sum(1 for r in results if r["status"] in ["failed", "error"])

    logger.info(f"\n{'='*60}")
    logger.info(f"Completed in {elapsed/60:.1f} minutes")
    logger.info(f"Success: {success}, Skipped: {skipped}, Failed: {failed}")



def main():
    parser = argparse.ArgumentParser(
        description="Submit L5TTPC trace simulations with Chindemi params",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all to SLURM (overnight run)
  python submit_l5ttpc_traces.py
  
  # Test with 2 pairs
  python submit_l5ttpc_traces.py --max-pairs 2
  
  # Run locally with multiprocessing
  python submit_l5ttpc_traces.py --execution-mode cpu --workers 30
        """
    )
    
    parser.add_argument("--params", choices=["v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15", "v16", "v17", "v18", "v19", "chindemi"], default="v1",
                        help="Parameter set to use: v1=DHURUVA_PARAMS, ..., v15=DHURUVA_PARAMS_V15, chindemi=CHINDEMI_PARAMS (default: v1)")
    parser.add_argument("--execution-mode", choices=["slurm", "cpu"], default="slurm",
                        help="Execution mode (default: slurm)")
    parser.add_argument("--max-pairs", type=int, default=None,
                        help="Max pairs to process (for testing)")
    parser.add_argument("--freq", type=str, default=None,
                        help="Filter by frequency, e.g., '10Hz'")
    parser.add_argument("--delays", type=str, nargs="+", default=None,
                        help="Filter by delays, e.g., '0' '-30' '10'")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Custom output directory")
    parser.add_argument("--simulations-dir", type=str, default=None,
                        help="Custom simulations directory (default depends on --params)")
    
    # CPU mode options
    parser.add_argument("--workers", type=int, default=30,
                        help="Workers for CPU mode (default: 30)")
    
    # SLURM options
    parser.add_argument("--timeout-hours", type=int, default=2,
                        help="SLURM timeout in hours (default: 2)")
    parser.add_argument("--mem-gb", type=int, default=8,
                        help="Memory per job in GB (default: 8)")
    parser.add_argument("--cpus-per-task", type=int, default=10,
                        help="CPUs per SLURM job (default: 10)")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Simulations per SLURM job (default: 10)")
    parser.add_argument("--parallel-sims", type=int, default=4,
                        help="Parallel simulations to pack inside one SLURM job (default: 4)")
    parser.add_argument("--fastforward", type=float, default=None,
                        help="Fastforward begin point in ms: run to this point first, snap synapses, then record")
    parser.add_argument("--epg-variant", "--variant",
                        choices=["epg", "epg_dhuruva", "epg_dhuruva_full", "epg_dhuruva_custom_ratio"],
                        default="epg_dhuruva",
                        help="Parameter generator module to use (default: epg_dhuruva)")
    parser.add_argument("--rho-ratio", "--rho_ratio", type=str, default=None,
                        help="Depressed:potentiated percentage split for epg_dhuruva_custom_ratio, e.g. 45:55")
    parser.add_argument("--epg-full-pairs-file", type=str, default=DEFAULT_EPG_FULL_PAIRS_FILE,
                        help="Pairs file used by epg_dhuruva_full to compute the pooled global median")
    parser.add_argument("--slurm-account", default="ctb-emuller",
                        help="SLURM account (default: ctb-emuller)")
    
    args = parser.parse_args()
    if args.epg_variant == "epg_dhuruva_custom_ratio":
        parse_rho_ratio(args.rho_ratio)
    elif args.rho_ratio is not None:
        parser.error("--rho-ratio can only be used with --epg-variant=epg_dhuruva_custom_ratio")
    
    # Select parameter set
    global DHURUVA_PARAMS
    default_simulations_dir = SIMULATIONS_DIR
    if args.params == "v19":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V19
        default_subdir = "DHURUVA_PARAMS_V19"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v18":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V18
        default_subdir = "DHURUVA_PARAMS_V18"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v17":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V17
        default_subdir = "DHURUVA_PARAMS_V17"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v16":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V16
        default_subdir = "DHURUVA_PARAMS_V16"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v15":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V15
        default_subdir = "DHURUVA_PARAMS_V15"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v14":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V14
        default_subdir = "DHURUVA_PARAMS_V14"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v13":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V13
        default_subdir = "DHURUVA_PARAMS_V13"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v12":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V12
        default_subdir = "DHURUVA_PARAMS_V12"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v11":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V11
        default_subdir = "DHURUVA_PARAMS_V11"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v10":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V10
        default_subdir = "DHURUVA_PARAMS_V10"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v9":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V9
        default_subdir = "DHURUVA_PARAMS_V9"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v8":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V8
        default_subdir = "DHURUVA_PARAMS_V8"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v7":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V7
        default_subdir = "DHURUVA_PARAMS_V7"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v6":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V6
        default_subdir = "DHURUVA_PARAMS_V6"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    elif args.params == "v5":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V5
        default_subdir = "DHURUVA_PARAMS_V5"
    elif args.params == "v4":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V4
        default_subdir = "DHURUVA_PARAMS_V4"
    elif args.params == "v3":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V3
        default_subdir = "DHURUVA_PARAMS_V3"
    elif args.params == "v2":
        DHURUVA_PARAMS = DHURUVA_PARAMS_V2
        default_subdir = "DHURUVA_PARAMS_V2"
    elif args.params == "chindemi":
        DHURUVA_PARAMS = CHINDEMI_PARAMS
        default_subdir = "CHINDEMI_PARAMS"
        default_simulations_dir = SIMULATIONS_DIR_STDP
    else:
        default_subdir = "DHURUVA_PARAMS"

    # Output directory
    output_dir = args.output_dir or os.path.join(OUTPUT_BASE_DIR, default_subdir)
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Find simulations by scanning for prefire_simulation_config.json files
    simulations_path = Path(args.simulations_dir if args.simulations_dir is not None else default_simulations_dir)
    jobs = find_jobs_from_prefire_configs(simulations_path)

    if args.freq:
        jobs = [(pd, fd) for pd, fd in jobs if args.freq in fd]

    if args.delays:
        target_delays = [f"{d}ms" for d in args.delays]
        jobs = [(pd, fd) for pd, fd in jobs
                if any(fd.endswith(f"_{delay}") for delay in target_delays)]

    if args.max_pairs:
        seen_pairs = []
        filtered = []
        for pd, fd in jobs:
            if pd not in seen_pairs:
                seen_pairs.append(pd)
            if len(seen_pairs) <= args.max_pairs:
                filtered.append((pd, fd))
        jobs = filtered

    logger.info(f"Total simulations: {len(jobs)}")

    metadata = {
        "params": DHURUVA_PARAMS,
        "jobs": jobs,
        "output_dir": output_dir,
        "fastforward": args.fastforward,
        "epg_variant": args.epg_variant,
        "rho_ratio": args.rho_ratio,
        "epg_full_pairs_file": args.epg_full_pairs_file,
    }
    with open(os.path.join(output_dir, "job_metadata.pkl"), "wb") as f:
        pickle.dump(metadata, f)

    if args.execution_mode == "slurm":
        submit_slurm_jobs(
            jobs, output_dir,
            timeout_hours=args.timeout_hours,
            mem_gb=args.mem_gb,
            cpus_per_task=args.cpus_per_task,
            batch_size=args.batch_size,
            slurm_account=args.slurm_account,
            parallel_sims=args.parallel_sims,
            fastforward=args.fastforward,
            epg_variant=args.epg_variant,
            rho_ratio=args.rho_ratio,
            epg_full_pairs_file=args.epg_full_pairs_file,
        )
    else:
        run_cpu_mode(jobs, output_dir, fit_params=DHURUVA_PARAMS, workers=args.workers,
                     fastforward=args.fastforward, epg_variant=args.epg_variant,
                     epg_full_pairs_file=args.epg_full_pairs_file, rho_ratio=args.rho_ratio)


if __name__ == "__main__":
    main()
