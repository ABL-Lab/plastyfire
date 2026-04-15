"""
Extra Parameter Generator variant using a configurable depressed/potentiated ratio.

This variant matches `epg_dhuruva` except for rho assignment:
it sorts synapses within each connection by sampled `gmax0_AMPA` and assigns the
top requested fraction to potentiated (`rho0_GB = 1`) and the remainder to
depressed (`rho0_GB = 0`).
"""

import os

import numpy as np

from plastyfire.epg_dhuruva import (
    ParamsGenerator as BaseParamsGenerator,
    _get_ltpltd_params,
    _get_synapse_location,
)


RHO_RATIO_ENV = "PLASTYFIRE_EPG_RHO_RATIO"


def parse_rho_ratio(ratio_spec):
    """Parse a `depressed:potentiated` percentage string like `45:55`."""
    if ratio_spec is None:
        raise ValueError(
            f"Missing {RHO_RATIO_ENV}. Set it to a depressed:potentiated ratio such as 45:55."
        )

    parts = [part.strip() for part in str(ratio_spec).split(":")]
    if len(parts) != 2 or not all(parts):
        raise ValueError(
            f"Invalid rho ratio '{ratio_spec}'. Expected a depressed:potentiated ratio such as 45:55."
        )

    try:
        depressed_ratio, potentiated_ratio = (float(part) for part in parts)
    except ValueError as exc:
        raise ValueError(
            f"Invalid rho ratio '{ratio_spec}'. Both values must be numeric percentages."
        ) from exc

    if depressed_ratio < 0 or potentiated_ratio < 0:
        raise ValueError(
            f"Invalid rho ratio '{ratio_spec}'. Percentages must be non-negative."
        )
    if not np.isclose(depressed_ratio + potentiated_ratio, 100.0):
        raise ValueError(
            f"Invalid rho ratio '{ratio_spec}'. Depressed and potentiated percentages must sum to 100."
        )
    return depressed_ratio, potentiated_ratio


def _rounded_count(total, percentage):
    """Return the closest integer count for a requested percentage."""
    return int(np.floor((total * percentage / 100.0) + 0.5))


class ParamsGenerator(BaseParamsGenerator):
    """Generate synapse parameters using a configurable gmax-ranked rho split."""

    def __init__(self, circuit, node_pop, edge_pop, extra_recipe_path, k_u=0.2, k_gsyn=2):
        super().__init__(circuit, node_pop, edge_pop, extra_recipe_path, k_u=k_u, k_gsyn=k_gsyn)
        self.depressed_ratio, self.potentiated_ratio = parse_rho_ratio(os.environ.get(RHO_RATIO_ENV))

    def generate_params(self, pre_gid, post_gid):
        """Generate parameters for one pair using a requested depressed/potentiated ratio."""
        raw_params, pathway_recipe = self._sample_raw_params(pre_gid, post_gid)
        ranked_syn_ids = sorted(
            raw_params,
            key=lambda syn_id: (raw_params[syn_id]["params"]["gmax0_AMPA"], syn_id),
        )
        potentiated_count = min(
            len(ranked_syn_ids),
            max(0, _rounded_count(len(ranked_syn_ids), self.potentiated_ratio)),
        )
        potentiated_syn_ids = set(ranked_syn_ids[-potentiated_count:]) if potentiated_count else set()

        syn_params = {}
        for syn_id, data in raw_params.items():
            params = dict(data["params"])
            branch_type = data["branch_type"]
            rho0 = 1 if syn_id in potentiated_syn_ids else 0
            params.update(_get_ltpltd_params(params["Use0_TM"], params["gmax0_AMPA"], self.k_u, self.k_gsyn, rho0))
            params["gmax_NMDA"] = params["gmax0_AMPA"] * pathway_recipe["gsynSRSF"]
            params["loc"] = _get_synapse_location(branch_type)
            syn_params[syn_id] = params
        return syn_params
