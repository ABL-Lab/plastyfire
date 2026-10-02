"""Run a gpu_v1x kernel with its fitter fit_v6 replaced by the module named in env FITTER (default fit_v7).

Usage (in a Slurm job, from plastyfire/): FITTER=fit_v7 python -u glusynapse_v2/rho_redesign/fit_launch.py KERNEL.py <kernel args>
The kernel's `import fit_v6` then returns the FITTER module, and all its patches (F.pack, unpack5, gpu_v6_rho, BatchV2,
MV, _CFG, _M) land on that module. No kernel file is edited (SPEC_FITPROC.md, part 2). run_stage4_fit.sh calls this when
env FITTER is set.
"""
import importlib, os, runpy, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]
if len(sys.argv) < 2 or not sys.argv[1].endswith(".py"):
    raise SystemExit("usage: fit_launch.py KERNEL.py <kernel args>")
kernel = sys.argv.pop(1)
fitter = os.environ.get("FITTER", "fit_v7")
sys.modules["fit_v6"] = importlib.import_module(fitter)
print(f"fit_launch: fit_v6 -> {fitter} ({sys.modules['fit_v6'].__file__}); kernel {kernel}", flush=True)
sys.argv[0] = kernel
runpy.run_path(kernel, run_name="__main__")
