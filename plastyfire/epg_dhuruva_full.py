"""
Extra Parameter Generator variant using a global conductance median.

This variant matches `epg_dhuruva` except for rho assignment:
it samples `gmax0_AMPA` for all synapses across a configured pair set
(default: data/pairs_n100.txt), computes one global median over that pooled
distribution, and then assigns rho0 for each synapse by thresholding against
that single median instead of the per-pair median.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from neurom import NeuriteType

from plastyfire.epg_dhuruva import (
    MAX_SEED,
    BRANCH_TYPE_OFFSET,
    _get_covariance_matrix,
    _get_distributions,
    _normtodist,
    _get_ltpltd_params,
)


DEFAULT_PAIRS_FILE = (
    "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/data/pairs_n100.txt"
)


class ParamsGenerator(object):
    """Generate synapse parameters using a single pooled gmax median."""

    def __init__(self, circuit, node_pop, edge_pop, extra_recipe_path, k_u=0.2, k_gsyn=2):
        self.circuit = circuit
        self.node_pop = node_pop
        self.edge_pop = edge_pop
        self.extra_recipe = pd.read_csv(extra_recipe_path, index_col=[0, 1], skipinitialspace=True)
        self.extra_recipe.columns = self.extra_recipe.columns.str.strip()
        self.extra_recipe.index = self.extra_recipe.index.set_levels(
            [level.str.strip() for level in self.extra_recipe.index.levels]
        )
        for col in self.extra_recipe.columns:
            if self.extra_recipe[col].dtype == object:
                self.extra_recipe[col] = self.extra_recipe[col].str.strip()
        self.k_u = k_u
        self.k_gsyn = k_gsyn
        self.namelst = ["u", "d", "f", "nrrp", "gsyn", "spinevol"]
        self.paramlst = ["Use0_TM", "Dep_TM", "Fac_TM", "Nrrp_TM", "gmax0_AMPA", "volume_CR"]
        self.pairs_file = os.environ.get("PLASTYFIRE_EPG_FULL_PAIRS_FILE", DEFAULT_PAIRS_FILE)
        self._raw_cache = {}
        self._global_gsyn_median = None
        self._pairs = self._load_pairs()

    def _load_pairs(self):
        pairs = []
        with open(self.pairs_file) as f:
            for line in f:
                line = line.strip()
                if line:
                    pre_gid, post_gid = map(int, line.split())
                    pairs.append((pre_gid, post_gid))
        return pairs

    def _sample_raw_params_for_pair(self, pre_gid, post_gid):
        key = (pre_gid, post_gid)
        if key in self._raw_cache:
            return self._raw_cache[key]

        pre_mtype = self.circuit.nodes[self.node_pop].get(pre_gid, "mtype")
        post_mtype = self.circuit.nodes[self.node_pop].get(post_gid, "mtype")
        pathway_recipe = self.extra_recipe.loc[pre_mtype, post_mtype]
        distlst = [
            _get_distributions(
                pathway_recipe[f"{name}Dist"],
                pathway_recipe[name],
                pathway_recipe[f"{name}SD"],
            )
            for name in self.namelst
        ]

        syns = self.circuit.edges[self.edge_pop].pair_edges(pre_gid, post_gid, "afferent_section_type")
        raw_params = {}
        for syn_id, branch_type in syns.items():
            cov = _get_covariance_matrix(pathway_recipe)
            np.random.seed(np.mod(syn_id, MAX_SEED))
            sample_normal = stats.multivariate_normal.rvs(cov=cov)
            sample_params = map(_normtodist, distlst, sample_normal)
            params = dict(zip(self.paramlst, sample_params))
            raw_params[syn_id] = {"params": params, "branch_type": branch_type, "pathway_recipe": pathway_recipe}

        self._raw_cache[key] = raw_params
        return raw_params

    def _get_global_median(self):
        if self._global_gsyn_median is not None:
            return self._global_gsyn_median

        all_gsyn = []
        for pre_gid, post_gid in self._pairs:
            raw_params = self._sample_raw_params_for_pair(pre_gid, post_gid)
            if raw_params:
                all_gsyn.extend(raw["params"]["gmax0_AMPA"] for raw in raw_params.values())

        if not all_gsyn:
            raise ValueError(f"No synapses found while computing global median from {self.pairs_file}")

        self._global_gsyn_median = float(np.median(np.asarray(all_gsyn, dtype=np.float64)))
        return self._global_gsyn_median

    def generate_params(self, pre_gid, post_gid):
        """Generate parameters for one pair using the pooled global gmax median."""
        raw_params = self._sample_raw_params_for_pair(pre_gid, post_gid)
        gsyn_median = self._get_global_median()

        syn_params = {}
        for syn_id, data in raw_params.items():
            params = data["params"]
            branch_type = data["branch_type"]
            pathway_recipe = data["pathway_recipe"]
            rho0 = 1 if params["gmax0_AMPA"] >= gsyn_median else 0
            params.update(_get_ltpltd_params(params["Use0_TM"], params["gmax0_AMPA"], self.k_u, self.k_gsyn, rho0))
            params["gmax_NMDA"] = params["gmax0_AMPA"] * pathway_recipe["gsynSRSF"]
            if branch_type + BRANCH_TYPE_OFFSET == NeuriteType.basal_dendrite:
                params["loc"] = "basal"
            elif branch_type + BRANCH_TYPE_OFFSET == NeuriteType.apical_dendrite:
                params["loc"] = "apical"
            else:
                raise ValueError("Unknown neurite type")
            syn_params[syn_id] = params
        return syn_params
