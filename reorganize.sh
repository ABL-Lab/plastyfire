#!/usr/bin/env bash
# =============================================================================
# reorganize.sh — Plastyfire codebase cleanup script
# Generated: 2026-02-27
#
# Run from the plastyfire project root:
#   cd /home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire
#   bash reorganize.sh
#
# The script is idempotent — safe to rerun.
# All deletions are echoed before execution and gated behind DRY_RUN.
# Set DRY_RUN=0 to actually execute.
# =============================================================================

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
DRY_RUN="${DRY_RUN:-1}"   # default: dry run (0 = execute for real)

info()  { echo "[INFO]  $*"; }
warn()  { echo "[WARN]  $*"; }
doit()  {
    if [ "$DRY_RUN" = "1" ]; then
        echo "[DRY]   $*"
    else
        echo "[RUN]   $*"
        eval "$*"
    fi
}

echo "============================================================"
echo "  Plastyfire Reorganization"
if [ "$DRY_RUN" = "1" ]; then
    echo "  MODE: DRY RUN (set DRY_RUN=0 to apply changes)"
else
    echo "  MODE: LIVE — changes will be applied!"
fi
echo "============================================================"
echo ""

# ------------------------------------------------------------
# 1. Create new top-level directories
# ------------------------------------------------------------
info "=== 1. Creating new directories ==="
for d in tests tests/CICR_fitting scripts figures models data data/checkpoints; do
    if [ ! -d "$ROOT/$d" ]; then
        doit "mkdir -p '$ROOT/$d'"
    else
        info "  Already exists: $d/"
    fi
done
echo ""

# ------------------------------------------------------------
# 2. Move test_*.py files from root → tests/
# ------------------------------------------------------------
info "=== 2. Moving test files to tests/ ==="
TEST_FILES=(
    test_cpre.py
    test_cpre_cpost_analytical.py
    test_fitness_schedule.py
    test_jax.py
    test_pickle_reader.py
    test_v6_params.py
)
for f in "${TEST_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/tests/$f'"
    else
        warn "  Not found (already moved?): $f"
    fi
done

# Shell test scripts from CICR_fitting/
CICR_TEST_SHELLS=(
    CICR_fitting/test_10Hz.sh
    CICR_fitting/test_trace_cpu.sh
    CICR_fitting/test_trace_slurm.sh
)
for f in "${CICR_TEST_SHELLS[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/tests/CICR_fitting/$(basename "$f")'"
    else
        warn "  Not found (already moved?): $f"
    fi
done
echo ""

# ------------------------------------------------------------
# 3. Move debug/utility scripts from root → scripts/
# ------------------------------------------------------------
info "=== 3. Moving debug/utility scripts to scripts/ ==="
DEBUG_SCRIPTS=(
    debug_mismatch.py
    debug_prefire_process.py
    compare_fresh_vs_existing.py
    compare_sim_pickles.py
    check_hash.py
    inspect_pickle.py
    inspect_rho.py
    recover_cached_results.py
    find_ephys.py
)
for f in "${DEBUG_SCRIPTS[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/scripts/$f'"
    else
        warn "  Not found (already moved?): $f"
    fi
done
echo ""

# ------------------------------------------------------------
# 4. Move generated PNG outputs → figures/
# ------------------------------------------------------------
info "=== 4. Moving generated PNGs to figures/ ==="
PNG_FILES=(
    stdp_curve_10Hz.png
    stdp_curve_10Hz_v2.png
    stdp_curve_10Hz_v3.png
    stdp_curve_10Hz_v4.png
    stdp_curve_10Hz_v5.png
    stdp_curve_10Hz_v7.png
    output.png
    test_plot_debug.png
    effcai_comparison.png
    epsp_ratio_comparison.png
    plot_example.png
)
for f in "${PNG_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/figures/$f'"
    else
        warn "  Not found (already moved?): $f"
    fi
done
echo ""

# ------------------------------------------------------------
# 5. Move trained model/scaler artifacts → models/
# ------------------------------------------------------------
info "=== 5. Moving model artifacts to models/ ==="
MODEL_FILES=(
    fitness_model.pth
    scaler_x.joblib
    scaler_y.joblib
    target_cols.joblib
)
for f in "${MODEL_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/models/$f'"
    else
        warn "  Not found (already moved?): $f"
    fi
done
echo ""

# ------------------------------------------------------------
# 6. Move large data files out of the Python package dir → data/
# ------------------------------------------------------------
info "=== 6. Moving large data files out of plastyfire/ package ==="
LARGE_PKG_FILES=(
    plastyfire/simulation_872a4f3ab18a.pkl
    plastyfire/simulation_epsp_df.csv
    plastyfire/epsp_results.pkl
    plastyfire/topology_analysis.log
)
for f in "${LARGE_PKG_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/data/$(basename "$f")'"
    else
        warn "  Not found (already moved?): $f"
    fi
done

# Move checkpoint PKLs out of the package to data/checkpoints/
info "   Moving checkpoint PKLs..."
CHECKPOINT_FILES=(
    plastyfire/checkpoint.pkl
    plastyfire/checkpoint.pkl.tmp
    plastyfire/checkpoint3.pkl
    plastyfire/checkpoint_flippy.pkl
    plastyfire/checkpoint_full.pkl
    plastyfire/checkpoint_mrk97_07_mrk97_08.pkl
    plastyfire/checkpoint_mrk97_08.pkl
    plastyfire/checkpoint_weighted_ini.pkl
)
for f in "${CHECKPOINT_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "mv '$ROOT/$f' '$ROOT/data/checkpoints/$(basename "$f")'"
    else
        warn "  Not found (already moved?): $f"
    fi
done
echo ""

# ------------------------------------------------------------
# 7. Delete backup / dead files
# ------------------------------------------------------------
info "=== 7. Deleting backup and dead files ==="
DELETE_FILES=(
    plastyfire/evaluator.py.backup
    plastyfire/evaluator_old.py
    plastyfire/modelfitter.py.backup
    plastyfire/run_slurm_plasty_test.sh.backup
    # checkpoint.pkl.tmp already moved above
)
for f in "${DELETE_FILES[@]}"; do
    if [ -f "$ROOT/$f" ]; then
        doit "rm '$ROOT/$f'"
    else
        warn "  Not found (already deleted?): $f"
    fi
done

# Delete empty topology log duplicates
info "   Deleting empty topology_analysis.log files..."
for f in CICR_fitting/topology_analysis.log new_fitting/topology_analysis.log; do
    if [ -f "$ROOT/$f" ] && [ ! -s "$ROOT/$f" ]; then
        doit "rm '$ROOT/$f'"
    elif [ -f "$ROOT/$f" ]; then
        warn "  $f is NOT empty, skipping: $(du -h "$ROOT/$f" | cut -f1)"
    else
        warn "  Not found: $f"
    fi
done

echo ""

# ------------------------------------------------------------
# 8. Remove the root-level generated CSV (epsp_ratio_comparison.csv)
# ------------------------------------------------------------
info "=== 8. Moving generated CSV outputs to data/ ==="
if [ -f "$ROOT/epsp_ratio_comparison.csv" ]; then
    doit "mv '$ROOT/epsp_ratio_comparison.csv' '$ROOT/data/epsp_ratio_comparison.csv'"
fi
echo ""

# ------------------------------------------------------------
# 9. Remove duplicate files from CICR_fitting/ and new_fitting/
#    that are also in jax_fitting/ (the canonical latest pipeline)
# ------------------------------------------------------------
info "=== 9. Deduplicating older pipelines vs jax_fitting/ ==="

# Files that exist in both CICR_fitting/ and jax_fitting/ — delete from CICR_fitting/
# (these share identical or near-identical content; jax_fitting/ wins)
CICR_DUPS_IN_JAX=(
    cicr_common.py
    cicr_primed_junctional.py
    diagnose_model.py
    gb_only.py
    plot_cicr_diagnostic.py
    validate_jax_vs_neuron.py
)
info "   Removing duplicates from CICR_fitting/ (jax_fitting/ is canonical)..."
for f in "${CICR_DUPS_IN_JAX[@]}"; do
    src="$ROOT/CICR_fitting/$f"
    canon="$ROOT/jax_fitting/$f"
    if [ -f "$src" ] && [ -f "$canon" ]; then
        doit "rm '$src'"
    elif [ -f "$src" ] && [ ! -f "$canon" ]; then
        warn "  $f not in jax_fitting — keeping CICR_fitting/$f"
    else
        warn "  Not found: CICR_fitting/$f"
    fi
done

# Files that exist in both new_fitting/ and jax_fitting/ — delete from new_fitting/
NEW_DUPS_IN_JAX=(
    cicr_common.py
    cicr_primed_junctional.py
    diagnose_model.py
    gb_only.py
    plot_cicr_diagnostic.py
)
info "   Removing duplicates from new_fitting/ (jax_fitting/ is canonical)..."
for f in "${NEW_DUPS_IN_JAX[@]}"; do
    src="$ROOT/new_fitting/$f"
    canon="$ROOT/jax_fitting/$f"
    if [ -f "$src" ] && [ -f "$canon" ]; then
        doit "rm '$src'"
    elif [ -f "$src" ] && [ ! -f "$canon" ]; then
        warn "  $f not in jax_fitting — keeping new_fitting/$f"
    else
        warn "  Not found: new_fitting/$f"
    fi
done

# Root-level fit_params_CICR.py — older copy; new_fitting has a larger version but
# jax_fitting is canonical, so archive the root copy
if [ -f "$ROOT/fit_params_CICR.py" ]; then
    doit "mv '$ROOT/fit_params_CICR.py' '$ROOT/scripts/fit_params_CICR_old.py'"
    info "   (archived root fit_params_CICR.py → scripts/fit_params_CICR_old.py)"
fi

# Root-level README_CICR_phenom_fitting.md — older copy; archive it
if [ -f "$ROOT/README_CICR_phenom_fitting.md" ]; then
    doit "mv '$ROOT/README_CICR_phenom_fitting.md' '$ROOT/scripts/README_CICR_phenom_fitting_old.md'"
    info "   (archived root README_CICR_phenom_fitting.md → scripts/)"
fi
echo ""

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
echo "============================================================"
if [ "$DRY_RUN" = "1" ]; then
    echo "  DRY RUN COMPLETE — no files were actually changed."
    echo "  Review the actions above, then run:"
    echo "    DRY_RUN=0 bash reorganize.sh"
else
    echo "  REORGANIZATION COMPLETE!"
fi
echo "============================================================"
