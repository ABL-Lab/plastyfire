"""
Config file class
author: András Ecker; last update 05.2024
"""

import os
import yaml


class BaseConfig(object):
    """Class to store common config parameters about the simulations"""
    def __init__(self, config_path):
        """YAML config file based constructor"""
        self._config_path = config_path
        with open(config_path, "r") as f:
            self._config = yaml.load(f, Loader=yaml.SafeLoader)

    @property
    def config(self):
        return self._config

    @property
    def circuit_config(self):
        return self.config["circuit"]["config"]

    @property
    def pre_simwriter_circuit_config(self):
        """Circuit config used for cell/pair selection (find_pairs()); falls back to
        `circuit_config` if `pre_simwriter_config` isn't set. Lets pair selection run
        against the original (unmodified) edges, independent of `circuit_config`
        (which may point at plasticity-modified edges written by a separate step)."""
        return self.config["circuit"].get("pre_simwriter_config", self.circuit_config)

    @property
    def node_set(self):
        return self.config["circuit"]["node_set"]

    @property
    def target(self):
        return self.config["circuit"]["target"]

    @property
    def node_pop(self):
        return self.config["circuit"]["node_pop"]

    @property
    def edge_pop(self):
        return self.config["circuit"]["edge_pop"]

    @property
    def sims_dir(self):
        return self.config["sims_dir"]

    @property
    def env(self):
        return self.config["simulator"]["env"]

    @property
    def run(self):
        return self.config["simulator"]["run"]


class Config(BaseConfig):
    """Class to store config parameters about finding C_pre and C_post for all synapses"""
    @property
    def sim_config(self):
        return os.path.join(self.sims_dir, "simulation_config.json")

    @property
    def use_extra_recipe(self):
        return self.config["use_extra_recipe"]

    @property
    def out_dir(self):
        return self.config["out_dir"]

    @property
    def fit_params(self):
        return self.config["fit_params"]


class OptConfig(BaseConfig):
    """Class to store config parameters about the optimization of model parameters"""
    @property
    def label(self):
        return self.config["label"]

    @property
    def seed(self):
        return self.config["seed"]

    @property
    def npairs(self):
        return self.config["npairs"]

    @property
    def out_dir(self):
        return os.path.join(self.sims_dir, "fitting", "n%i" % self.npairs,
                            "seed%i" % self.seed, self.label, "simulations")

    @property
    def pre_mtype(self):
        return self.config["pregid_conf"]["mtype"]

    @property
    def post_mtype(self):
        return self.config["postgid_conf"]["mtype"]

    @property
    def max_dist(self):
        return self.config["geom_cons"]["max_dist"]

    @property
    def freq(self):
        return self.config["stimulus"]["freq"]

    @property
    def amp(self):
        return None if self.config["stimulus"]["amp"] == "find" else self.config["stimulus"]["amp"]

    @property
    def amp_min(self):
        return self.config["stimulus"]["amp_min"]

    @property
    def amp_max(self):
        return self.config["stimulus"]["amp_max"]

    @property
    def amp_lev(self):
        return self.config["stimulus"]["amp_lev"]

    @property
    def nspikes(self):
        return self.config["stimulus"]["nspikes"]

    @property
    def dt(self):
        return self.config["stimulus"]["dt"]

    @property
    def width(self):
        return self.config["stimulus"]["width"]

    @property
    def dt(self):
        return self.config["stimulus"]["dt"]

    @property
    def width(self):
        return self.config["stimulus"]["width"]

    @property
    def T(self):
        return self.config["stimulus"]["T"]

    @property
    def offset(self):
        return self.config["stimulus"]["offset"]

    @property
    def nreps(self):
        return self.config["stimulus"]["nreps"]

    @property
    def C01_duration(self):
        return self.config["stimulus"]["C01_duration"]

    @property
    def C02_duration(self):
        return self.config["stimulus"]["C02_duration"]

    @property
    def C01_T(self):
        return self.config["stimulus"]["C01_T"] if "C01_T" in self.config["stimulus"] else self.T

    @property
    def C02_T(self):
        return self.config["stimulus"]["C02_T"] if "C02_T" in self.config["stimulus"] else self.T

    @property
    def fastforward(self):
        return self.config["simulator"]["fastforward"] if "fastforward" in self.config["simulator"] else None

    @property
    def index_label(self):
        """Name of the index csv / yaml copy written next to the sims (default `label`). Set it to write
        extra protocols into an existing `label` folder without overwriting that folder's index"""
        return self.config.get("index_label", self.label)

    @property
    def prefire_only(self):
        """If True only the pairing (prefire) sim files are written: no simulation_config.json,
        prespikes.h5 or simulation.batch (C01/C02 test pulses of the full protocol)"""
        return self.config.get("prefire_only", False)

    @property
    def pairs_from(self):
        """Optional path to an existing `index_<label>.csv`: reuse its (pregid, postgid) pairs
        instead of running `find_pairs()` (keeps validation protocols on the same connections)"""
        return self.config.get("pairs_from")

    @property
    def pipette(self):
        """Optional extracellular stimulation mode (see `plastyfire.pipette`): instead of one connected
        presynaptic cell, a group of EXC cells recruited near a basal site and calibrated to a compound
        EPSP. Keys: path_dist [min, max] um, epsp_target mV, epsp_range [min, max] mV, max_pre"""
        return self.config.get("pipette")

    @property
    def calibration_circuit_config(self):
        """Circuit config used for the per-protocol current amplitude searches. Defaults to
        `circuit_config`, i.e. the emodels the protocol will actually be simulated with"""
        return self.config["circuit"].get("calibration_config", self.circuit_config)

    @property
    def protocols(self):
        """List of protocol dicts, each with keys: id, n_pre, pre_freq, n_post, freq, dt, dt_ref,
        nreps, T, width, amp (None -> search).

        If the yaml has no `protocols:` list, the classic freq x dt grid of `stimulus:` is expanded
        (n_pre = n_post = nspikes, ids "<f>Hz_<dt>ms"), which reproduces the original simwriter output.
        Otherwise every `protocols:` entry needs `id` and `dt` and inherits the rest from `stimulus:`.
        `dt` > 0 means pre before post. `dt_ref`: "ap" (default) measures dt to the first postsynaptic
        AP, "stim" to the onset of the first current pulse (Ebner 2019 convention); "ap_last" / "stim_last"
        to the last AP / pulse of the burst (post-pre bursts with dt to the closest AP, Nevian 2006).
        `amp_min`, `amp_max`, `amp_lev` override the amplitude search range of `stimulus:` per protocol.
        `burst_window_search: true` counts only APs inside the burst window during the amplitude search and
        allows later spikes (see `simwriter.burst_window_threshold_finder`).
        If n_pre == n_post and pre_freq == freq, pre spike i is paired with post AP i (Markram 1997);
        otherwise the n_pre pre spikes form their own train starting dt before the reference.

        Depolarising steps (Sjostrom 2007, Markram & Tsodyks 1996): `post_type: step` replaces the post
        pulse train by ONE `step_duration` ms current step per repetition (n_post: 1 = step, 0 = none);
        `pre_type: step` replaces the pre spike train by the spikes the presynaptic cell itself fires
        during a step of the same kind (n_pre: 1 = step, 0 = none). The step amplitude is `amp` (nA) if set,
        else the lowest amplitude on the [amp_min, amp_max] / amp_lev grid that fires >= `step_nspikes` APs
        in the step, accepted if <= `step_nspikes_max` (searched per cell, see
        `simwriter.step_rate_threshold_finder`). dt (dt_ref "stim") = pre step onset before post step onset.
        Defaults (post_type "pulse", pre_type "train") leave every other protocol unchanged.
        """
        stim = self.config["stimulus"]
        defaults = {"n_pre": stim["nspikes"], "n_post": stim["nspikes"], "dt_ref": "ap",
                    "nreps": stim["nreps"], "T": stim["T"], "width": stim["width"],
                    "amp": None if stim["amp"] == "find" else stim["amp"]}
        if "protocols" not in self.config:
            protos = [dict(defaults, id="%iHz_%ims" % (int(freq), int(dt)), freq=float(freq), dt=float(dt), grid=True)
                      for freq in stim["freq"] for dt in stim["dt"]]
        else:
            protos = []
            for entry in self.config["protocols"]:
                p = dict(defaults, freq=float(stim["freq"][0]), grid=False)
                p.update(entry)
                p["freq"], p["dt"] = float(p["freq"]), float(p["dt"])
                protos.append(p)
        for p in protos:
            p.setdefault("pre_freq", p["freq"])
            p["pre_freq"] = float(p["pre_freq"])
            assert p["dt_ref"] in ("ap", "stim", "ap_last", "stim_last"), \
                "dt_ref must be 'ap', 'stim', 'ap_last' or 'stim_last' (%s)" % p["id"]
            assert p.get("post_type", "pulse") in ("pulse", "step"), "post_type: 'pulse' or 'step' (%s)" % p["id"]
            assert p.get("pre_type", "train") in ("train", "step"), "pre_type: 'train' or 'step' (%s)" % p["id"]
            if "step" in (p.get("post_type"), p.get("pre_type")):
                assert p.get("step_duration", 0) > 0, "step protocols need step_duration (ms) (%s)" % p["id"]
                assert (p.get("post_type") != "step" or p["n_post"] in (0, 1)) and \
                    (p.get("pre_type") != "step" or p["n_pre"] in (0, 1)), \
                    "step: n_pre / n_post are 0 or 1 (one step per repetition) (%s)" % p["id"]
                assert "step_nspikes" in p, "step protocols need step_nspikes (%s)" % p["id"]
                assert not p.get("step_below_nA") or (p.get("post_type") == "step" and p.get("pre_type", "train") == "train"
                                                      and p.get("amp") is None), \
                    "step_below_nA: post step with a searched amp and a pre spike train (%s)" % p["id"]
        ids = [p["id"] for p in protos]
        assert len(ids) == len(set(ids)), "duplicate protocol ids: %s" % ids
        return protos
