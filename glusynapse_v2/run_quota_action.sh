#!/bin/bash
#SBATCH --account=def-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0:15:00
#SBATCH --output=/scratch/dhuruva/quota_scan/actions/%x_%j.out
# MEASURED 22132872 (plastyfitting tar, 368 GB): 37:19, 1.26 GB MaxRSS -> large tars 1600M, 1:00.
# Inode clean-up actions proposed in QUOTA_AUDIT.md. Run ONLY after user approval.
#   sbatch -J q_<ITEM> --export=ALL,ITEM=<item>,MODE=<mode> [--time=..] run_quota_action.sh
# MODE=dryrun (default): writes the item list + inode count to /scratch, changes nothing.
# MODE=tar:     list -> $A/<item>.tar, then <item>.count (tar entries) and <item>.expect (inodes).
# MODE=delete:  for tar items only if <item>.count == <item>.expect; A-items delete directly.
# MODE=copy / gc: see mmcr_venv / gc items.
# Sizing: no measured tar/delete job yet. Basis: quota scan 22127461 (find over 334K
# inodes: 2:16, MaxRSS 818 MB incl. large awk tables); tar streams (<100 MB RSS est).
# Time est. = 2 x bytes / 300 MB/s (write + verify read). Run ITEM=cooker first as the
# pilot, record seff here, then rescale the others.
# MEASURED round 1 (2026-10-01; all 1G/0:15 except copy 1G/0:30):
#   gc 22127977 0:49 455MB; empty_out 22127978 0:13 87MB; cooker tar 22127979 4:54 1014MB
#   (AT LIMIT; 10,619 inodes, ~15 GB, count==expect); cooker del 22127980 0:13 6MB;
#   sv del 22127981 0:20 95MB; tr del 22127982 0:24 96MB;
#   mmcr copy 22127983 14:20 1015MB (AT LIMIT; 56K files, 7 GB).
# => tar/copy/verify/swap: --mem=1536M; deletes: --mem=512M. Tar time = cooker 4:54 x
#    max(inode ratio, byte ratio) x 1.5, rounded up to 15 min (cooker ~51 MB/s incl. verify).
set -euo pipefail
D=/project/rrg-emuller/dhuruva
R=plastyfire/refitting_results/fitting/n120/seed20262009
S=$R/Sabrina_L5TTPC_L5TTPC_STDP/simulations
E=$R/Ebner2019_L23PC_L5TTPC/simulations
A=$D/archive_tars
O=/scratch/dhuruva/quota_scan/actions
: "${ITEM:?set ITEM}"; MODE=${MODE:-dryrun}
mkdir -p "$O"; L=$O/$ITEM.list
cd "$D"   # all list paths are relative to $D, so tars restore with: tar -xf X.tar -C $D

# raw prefire outputs of one hash in <pair>/<protocol>/ dirs (exact names: never matches
# -tr/-vseg/-ljp25 variants of a shorter hash)
hashlist() { local h=$1; shift
  find "$@" -mindepth 3 -maxdepth 3 \( -name "bluecellulab_results_$h" -o -name "pool_${h}_*" \
       -o -name "simulation_edges_$h.pkl" -o -name "simulation_$h.pkl" \); }

make_list() { case $ITEM in
  empty_out)  find $S $E -mindepth 3 -maxdepth 3 -type d -name out -empty -mmin +120 ;;   # -mmin: skip dirs a running simwriter just made
  sv)         hashlist delta-sv-prefire $S ;;
  tr)         hashlist delta-prefire-tr $S ;;
  cooker)     hashlist full_delta-cooker $S; hashlist delta-cooker $S
              find $S -mindepth 3 -maxdepth 3 -type f \( -name out.h5 -o -name soma.h5 -o -name rho.h5 \
                -o -name mcomplex.dat -o -name populations_offset.dat -o -name 'cx_S1nonbarrel_*.dat' \
                -o -name 'prefire_simulation_config_delta-cooker.json*' -o -name 'pydamus_*.log' \) -mmin +120
              hashlist full_delta-cooker-full $E; echo plastyfire/logs/nd_delta-cooker ;;
  optimizer)  hashlist optimizer $S ;;
  antic)      hashlist antic-delta-prefire $S ;;
  svtr)       hashlist delta-sv-prefire-tr $S ;;
  plastyfitting) find plastyfitting -mindepth 1 -maxdepth 1 ! -name .venv ;;
  am_extracted)  find plastyfire/analytical_method -mindepth 1 -maxdepth 1 -name 'extracted*' ;;
  v2x_rejected)  find plastyfire/glusynapse_v2/extracted -mindepth 1 -maxdepth 1 \( -name ebner_delta-sv \
                   -o -name markram_delta-sv -o -name markram_delta-cooker -o -name markram_delta-prefire-tr \
                   -o -name ebner_preview -o -name '*_delta-ljp25-prefire-vca' \) ;;
  delta_prefire) hashlist delta-prefire $S $E ;;
  gc|mmcr_venv)  : ;;
  *) echo "unknown ITEM $ITEM"; exit 2 ;;
esac; }

inodes() { tr '\n' '\0' < "$1" | du --inodes -cs --files0-from=- | tail -1 | cut -f1; }   # one du, one total

case $MODE in
dryrun)
  case $ITEM in
    gc) for r in MMCR_CL plastyfitting; do echo "$r"; git -C $r count-objects -v; done ;;
    mmcr_venv) echo "MMCR_CL/.venv inodes: $(du --inodes -s MMCR_CL/.venv | cut -f1)" ;;
    *) make_list > "$L"; echo "$ITEM: $(wc -l < "$L") top-level items, $(inodes "$L") inodes -> $L"
       head -5 "$L" ;;
  esac ;;
tar)
  if [ "$ITEM" = v2x_rejected ]; then   # ljp25 extracted npz: only after the raw offload succeeded
    OL=${OFFLOAD_LOG:-/scratch/dhuruva/plastyfire_offload/offload_ljp25_22128394.out}
    grep -q "left on project: 0;" "$OL" || { echo "offload not confirmed in $OL, abort"; tail -3 "$OL"; exit 3; }
    grep "left on project" "$OL"
  fi
  mkdir -p "$A"; make_list > "$L"; inodes "$L" > "$O/$ITEM.expect"
  if [ "$ITEM" = plastyfitting ]; then ~/.local/bin/uv pip freeze --python plastyfitting/.venv/bin/python > "$O/plastyfitting_requirements.txt"; fi
  tar -cf "$A/$ITEM.tar" -T "$L"
  tar -tf "$A/$ITEM.tar" | wc -l > "$O/$ITEM.count"
  echo "$ITEM: tar entries $(cat "$O/$ITEM.count"), inodes $(cat "$O/$ITEM.expect")"; ls -l "$A/$ITEM.tar" ;;
delete)
  case $ITEM in
    empty_out) make_list > "$L"; echo "$ITEM: $(wc -l < "$L") dirs"; xargs -a "$L" -d '\n' -r rmdir ;;   # rmdir refuses non-empty dirs
    sv|tr)     make_list > "$L"; echo "$ITEM: $(wc -l < "$L") items, $(inodes "$L") inodes"; xargs -a "$L" -d '\n' -r rm -rf ;;
    gc|mmcr_venv) echo "use MODE=$([ $ITEM = gc ] && echo gc || echo copy)"; exit 2 ;;
    *) [ "$(cat "$O/$ITEM.count")" = "$(cat "$O/$ITEM.expect")" ] || { echo "count != expect, abort"; exit 3; }
       xargs -a "$L" -d '\n' -r rm -rf
       if [ "$ITEM" = plastyfitting ]; then rm -rf plastyfitting/.venv; fi ;;
  esac; echo "deleted $ITEM" ;;
gc) git -C MMCR_CL gc; if [ -d plastyfitting/.git ]; then git -C plastyfitting gc; fi; echo gc done ;;
copy)   # mmcr_venv step 1: copy venv to /home (inodes not on /project); step 2 is manual (see audit)
  [ -e "$HOME/venvs/mmcr_cl" ] && { echo "$HOME/venvs/mmcr_cl exists, abort"; exit 3; }
  mkdir -p "$HOME/venvs"; cp -a MMCR_CL/.venv "$HOME/venvs/mmcr_cl"
  echo "src $(du --inodes -s MMCR_CL/.venv | cut -f1) dst $(du --inodes -s $HOME/venvs/mmcr_cl | cut -f1)" ;;
verify)  # mmcr_venv step 2: copy complete and usable?
  V=$HOME/venvs/mmcr_cl; [ -L MMCR_CL/.venv ] && { echo "already a symlink"; exit 3; }
  src=$(du --inodes -s MMCR_CL/.venv | cut -f1); dst=$(du --inodes -s "$V" | cut -f1)
  echo "src $src dst $dst"; [ "$src" = "$dst" ] || { echo "count mismatch, abort"; exit 3; }
  "$V/bin/python" -c 'import sys; print(sys.prefix)'
  echo "pyvenv.cfg:"; cat "$V/pyvenv.cfg"; grep -m1 "VIRTUAL_ENV=" "$V/bin/activate" || true
  "$V/bin/python" -c 'import numpy; print("numpy", numpy.__version__)' || echo "numpy not importable (non-fatal)" ;;
swap)    # mmcr_venv step 3: real dir -> symlink to the /home copy; rollback if python fails
  V=$HOME/venvs/mmcr_cl; cd "$D/MMCR_CL"
  [ -d .venv ] && [ ! -L .venv ] && [ ! -e .venv_old ] || { echo "unexpected .venv state, abort"; exit 3; }
  mv .venv .venv_old && ln -s "$V" .venv
  if .venv/bin/python -c 'import sys; print(sys.prefix); import numpy'; then
    rm -rf .venv_old; ls -l .venv; echo "swap done"
  else
    rm .venv && mv .venv_old .venv; echo "python via symlink failed, rolled back"; exit 4
  fi ;;
*) echo "unknown MODE $MODE"; exit 2 ;;
esac
echo "elapsed ${SECONDS}s"
