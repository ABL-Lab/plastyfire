# Shared environment for the v3 GPU objective, sourced by every run_*_v3.sh from the plastyfire root.
# Own venv glusynapse_v2/.venv_v3 (the running fits use plastyfire/.venv, untouched): the same pinned numpy 1.26.4 / scipy /
# pandas / jax 0.7.1 / nvidia-*-cu12 versions as plastyfire/.venv, plus numba 0.65.1 + numba-cuda 0.30.2 (CC wheelhouse pair; numba-cuda 0.28.2 fails to import with numba 0.67: NPDatetime),
# which replaces the built-in numba.cuda (it could not target sm_90 with the CUDA 12.9 pip NVVM) and finds NVVM via cuda-pathfinder.
source glusynapse_v2/.venv_v3/bin/activate
unset CUDA_HOME
# numba 0.65.1 here vs 0.67 in plastyfire/.venv: own cache dir, so batch_v2 njit(cache=True) files of the running fits are not rewritten
export NUMBA_CACHE_DIR=/lustre09/project/6070394/dhuruva/plastyfire/glusynapse_v2/.numba_cache_v3
python -c "import numba, numba.cuda as c; print('numba', numba.__version__, 'numba.cuda', c.__file__)"
