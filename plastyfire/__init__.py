"""plastyfire"""

from plastyfire.version import __version__
from plastyfire.config import Config, OptConfig

def __getattr__(name):
    if name == "ParamsGenerator":
        from plastyfire.epg_dhuruva import ParamsGenerator
        return ParamsGenerator
    raise AttributeError(f"module 'plastyfire' has no attribute {name!r}")

from plastyfire.simwriter import SimWriter, OptSimWriter
from plastyfire.ephysutils import Experiment
from plastyfire.simulator import spike_threshold_finder, c_pre_finder, c_post_finder, runconnectedpair

