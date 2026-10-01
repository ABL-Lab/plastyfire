#!/bin/bash
# Re-run the whole plasticity pipeline on a NEW L5 TTPC emodel with one command (login-node safe: only writes job
# scripts and calls sbatch; no python). Steps and measured sizes: glusynapse_v2/EMODEL_PIPELINE.md.
#
#   EMODEL_NAME=<name> EMODEL_DIR=<dir with cADpyr_*.hoc> [CFG_L5L5=.. CFG_L23L5=.. CFG_L23L23=..] \
#   [PATHWAYS=l5l5,l23l5,l23l23] [PILOT=1] [SKIPFIT=1] [DRYRUN=1] [DECISIONS_LOGGED=1] \
#       bash glusynapse_v2/run_new_emodel_chain.sh
#
#   EMODEL_NAME  alnum tag of the new emodel. Names everything: param hash delta-<name>-prefire-vseg (L5->L5),
#                ...-vseg-rs (L2/3 pathways), extracted dirs glusynapse_v2/extracted/<dataset>_delta-<name>-prefire-vca
#                (L5) and <dataset>_delta-<name>-prefire-vseg-rs (L2/3), circuit configs data/dhuruva_<name>[_l23l5|_l23l23]_circuit_config.json.
#                Existing names (delta, ljp25, sv, ...) are refused; og-delta outputs are never touched.
#   EMODEL_DIR   biophysical_neuron_models_dir of the new emodel (cADpyr_L5TPC.hoc, _L2TPC, _L3TPC). Its mechanisms
#                (new .mod files with NEW SUFFIX) must be compiled into DEES_cell_packages/x86_64 (MECHANISMS_PATH is
#                hard-coded in plastyfire/simwriter.py and simulator_edges.py).
#   CFG_*        circuit config per pathway. Default: generated from the og-delta configs by replacing only
#                biophysical_neuron_models_dir (same edges files: dhuruva_modified_edges*.h5 are emodel-independent).
#   PATHWAYS     subset of l5l5,l23l5,l23l23 (default all three; the joint fit needs all three).
#   PILOT=1      2 pairs x 2 protocols per pathway, small jobs, sets SKIPFIT=1 and L5_SETS=mk,sj01 (nulls off).
#   SKIPBASIS=1  skip the EPSP-basis arrays (compare steps then skip the basis csv).
#   SKIPFIT=1    stop after extraction + basis (+ compare when the name starts with ogtest).
#   DRYRUN=1     write the job scripts to /scratch/dhuruva/<name>/jobs and print the plan; submit nothing.
#   L5_MK_CAL_CFG circuit config the Markram (Sabrina base) stimuli are calibrated on. Default = CFG_L5L5 (new emodel). The og-delta Markram
#                workdirs were calibrated on data/dhuruva_modified_ion_channels_circuit_config.json (older emodel); set it to that file
#                only for a bit-identical reproduction test of og-delta.
#   L5_SETS      L5->L5 datasets (default mk,sj01,sj03,sj03r50,sj07; sj03 does NOT include the r50 bursts (add sj03r50); sj04 optional, not in the best fit).
#   DECISIONS_LOGGED=1  required for a full (non-pilot) run: CLAUDE.md rule 7, log runs over 30 CPU-h in DECISIONS.md first.
#   Fit settings (defaults = best joint fit C1Ajn_s5, 3 pathways): FIT_SET FIT_FREE FIT_DROPT FIT_SEEDFITS FIT_FIX FIT_MAXITER.
#
# MEASURED pilot ogtest (2 pairs x 2 protocols, 22135633-45): basis_l5 17:41 10.0 GB; pf_l5_A 3:46 29.3 GB; basis_l23l5 3:51 5.5 GB; pf_l23l5 2:30 16.0 GB;
# sims, caches, extractions, compares < 1 min and < 2 GB (extractions: pilot requests kept, MaxRSS not reported below 2 GB).
# Raw outputs (workdirs, prefire traces, cache, basis, logs, job scripts) go to /scratch/dhuruva/<name>/; only the npz go to
# glusynapse_v2/extracted/ (project inode quota). Every step is --skip-existing, so a rerun with RESUME=1 continues.
set -euo pipefail
ROOT=/project/rrg-emuller/dhuruva/plastyfire
V2=$ROOT/glusynapse_v2
: "${EMODEL_NAME:?set EMODEL_NAME}"
: "${EMODEL_DIR:?set EMODEL_DIR}"
NAME=$EMODEL_NAME
PATHWAYS=${PATHWAYS:-l5l5,l23l5,l23l23}
PILOT=${PILOT:-0}; SKIPFIT=${SKIPFIT:-0}; DRYRUN=${DRYRUN:-0}; RESUME=${RESUME:-0}
L5_SETS=${L5_SETS:-mk,sj01,sj03,sj03r50,sj07}
NULLS=${NULLS:-1}
W=${WORKROOT:-/scratch/dhuruva}/$NAME
if [ "$PILOT" = 1 ]; then SKIPFIT=1; L5_SETS=${PILOT_L5_SETS:-mk,sj01}; NULLS=0; fi
case ",$L5_SETS," in *,mk,*) ;; *) L5_SETS=mk,$L5_SETS;; esac   # mk (base Sabrina folder) holds simulation_config.json for basis/cache

# ---------------------------------------------------------------- checks
[[ $NAME =~ ^[A-Za-z0-9]+$ ]] || { echo "EMODEL_NAME must be alphanumeric"; exit 1; }
case $NAME in delta|ljp25|sv|cooker|optimizer|antic|og) echo "EMODEL_NAME $NAME is reserved (existing outputs)"; exit 1;; esac
[ -d "$EMODEL_DIR" ] || { echo "no EMODEL_DIR $EMODEL_DIR"; exit 1; }
case $(basename "$EMODEL_DIR") in SSCx-AAD-delta-emodels) [[ $NAME == ogtest* ]] || { echo "og-delta dir only allowed for a plumbing test named ogtest*"; exit 1; };; esac
for h in cADpyr_L5TPC.hoc cADpyr_L2TPC.hoc cADpyr_L3TPC.hoc; do [ -f "$EMODEL_DIR/$h" ] || echo "WARNING: $EMODEL_DIR/$h missing"; done
[ -f /project/rrg-emuller/dhuruva/DEES_cell_packages/x86_64/libnrnmech.so ] || echo "WARNING: DEES_cell_packages/x86_64/libnrnmech.so missing (mechanisms)"
if [ "$RESUME" != 1 ] && ls -d $V2/extracted/*_delta-${NAME}-prefire* >/dev/null 2>&1; then
    echo "extracted dirs for $NAME already exist; set RESUME=1 to continue (skip-existing) or choose a new name"; exit 1; fi
has() { case ",$PATHWAYS," in *,$1,*) return 0;; *) return 1;; esac; }
hasset() { case ",$L5_SETS," in *,$1,*) return 0;; *) return 1;; esac; }

# ---------------------------------------------------------------- circuit configs (the only place the emodel is referenced)
mkcfg() {   # template out
    if [ -f "$2" ]; then grep -q "\"$EMODEL_DIR\"" "$2" || echo "WARNING: $2 does not point at $EMODEL_DIR"; return; fi
    sed -E "s#(\"biophysical_neuron_models_dir\": )\"[^\"]+\"#\1\"$EMODEL_DIR\"#" "$1" > "$2"
    grep -q "\"$EMODEL_DIR\"" "$2" || { echo "config generation failed for $2"; exit 1; }
}
D=$ROOT/data
CFG_L5L5=${CFG_L5L5:-$D/dhuruva_${NAME}_circuit_config.json}
CFG_L23L5=${CFG_L23L5:-$D/dhuruva_${NAME}_l23l5_circuit_config.json}
CFG_L23L23=${CFG_L23L23:-$D/dhuruva_${NAME}_l23l23_circuit_config.json}
mkdir -p $W/{jobs,logs,pairs,configs,cache,refitting_results/out}
mkcfg $D/dhuruva_delta_circuit_config.json       $CFG_L5L5
has l23l5  && mkcfg $D/dhuruva_delta_l23l5_circuit_config.json $CFG_L23L5
has l23l23 && mkcfg $D/dhuruva_delta_l23l23_circuit_config.json $CFG_L23L23
EDGES_L5L5=$D/dhuruva_modified_edges.h5; EDGES_L23L5=$D/dhuruva_modified_edges_l23l5.h5; EDGES_L23L23=$D/dhuruva_modified_edges_l23l23.h5
for e in $EDGES_L5L5 $EDGES_L23L5 $EDGES_L23L23; do [ -f $e ] || { echo "missing edges $e"; exit 1; }; done

# ---------------------------------------------------------------- naming
HP=delta-${NAME}-prefire                      # hash stem; og-delta is delta-prefire
H5=${HP}-vseg; H23=${HP}-vseg-rs; HN=${HP}-vseg-rs-nulls
X=$V2/extracted
SEEDDIR=$W/refitting_results/fitting/n120/seed20262009
SIMS_L5=$SEEDDIR/Sabrina_L5TTPC_L5TTPC_STDP/simulations
SIMS_L23L5=$SEEDDIR/Ebner2019_L23PC_L5TTPC/simulations
SIMS_L23L23=$SEEDDIR/Zilberter2009_L23PC_L23PC/simulations
OLDIDX=$ROOT/refitting_results/fitting/n120/seed20262009
declare -A YAML=([mk]=Sabrina_L5TTPC_L5TTPC_STDP [sj01]=Ebner2019_L5TTPC_L5TTPC [sj03]=Sjostrom2003_L5TTPC_L5TTPC [sj03r50]=Sjostrom2003b_L5TTPC_L5TTPC [sj07]=Sjostrom2007_L5TTPC_L5TTPC [sj04]=Sjostrom2004_L5TTPC_L5TTPC)
declare -A DIRN=([mk]=markram [sj01]=ebner [sj03]=sj03 [sj03r50]=sj03r50 [sj07]=sj07 [sj04]=sj04)
declare -A IDS_FULL=(
 [mk]="10Hz_-50ms,10Hz_-30ms,10Hz_-10ms,10Hz_5ms,10Hz_10ms,10Hz_30ms,10Hz_50ms"
 [sj01]="sjostrom_0.1hz_dt+10ms,sjostrom_0.1hz_dt-10ms,sjostrom_10hz_dt+10ms,sjostrom_10hz_dt-10ms,sjostrom_20hz_dt+10ms,sjostrom_20hz_dt-10ms,sjostrom_40hz_dt+10ms,sjostrom_40hz_dt-10ms,sjostrom_50hz_dt+10ms,sjostrom_50hz_dt-10ms"
 [sj03]="sjostrom_0.1hz_dt-25ms,sjostrom_0.1hz_dt-50ms,sjostrom_0.1hz_dt-100ms,sjostrom_0.1hz_dt-120ms,sjostrom_0.1hz_dt-200ms,sjostrom_0.1hz_dt+25ms,sjostrom_0.1hz_dt+50ms,sjostrom_pre_only_0.1hz,sjostrom_post_only_0.1hz"
 [sj03r50]="sjostrom_burst5x20hz_r50_dt-120ms,sjostrom_burst5x20hz_r50_dt-200ms"
 [sj07]="sjostrom07_step200ms_pair,sjostrom07_step200ms_pre_only,sjostrom07_step200ms_post_only,sjostrom07_step200ms_1.2nA_pair"
 [sj04]="sjostrom04_dltd_step250ms")
declare -A IDS_PILOT=(
 [mk]="10Hz_-10ms,10Hz_10ms" [sj01]="sjostrom_10hz_dt+10ms,sjostrom_10hz_dt-10ms"
 [sj03]="sjostrom_pre_only_0.1hz,sjostrom_post_only_0.1hz" [sj03r50]="sjostrom_burst5x20hz_r50_dt-120ms,sjostrom_burst5x20hz_r50_dt-200ms"
 [sj07]="sjostrom07_step200ms_pair,sjostrom07_step200ms_post_only" [sj04]="sjostrom04_dltd_step250ms")
ids() { if [ "$PILOT" = 1 ]; then echo "${IDS_PILOT[$1]}"; else echo "${IDS_FULL[$1]}"; fi; }
L23L5_IDS_PILOT="sjostrom_50hz_dt+10ms,letzkus_1ap_dt+10ms"
L23L23_IDS_PILOT="zilberter_1ap_dt+10ms,zilberter_1ap_dt-10ms"
NULL_IDS="letzkus_nopost,letzkus_3ap_200hz_dt-500ms"

# ---------------------------------------------------------------- pair lists (same pairs as og-delta: comparability)
NP=0; [ "$PILOT" = 1 ] && NP=2
if [ "$PILOT" = 1 ]; then head -$((NP+1)) $ROOT/configs/Sjostrom2007_subset24_pairs.csv > $W/pairs/l5.csv; CSV_L5=$W/pairs/l5.csv
else CSV_L5=$ROOT/configs/Sjostrom2007_subset24_pairs.csv; fi
mkpairs() {   # old index csv -> $W/pairs/<tag>.csv (pair,pregid,postgid; first NP unique pairs when pilot)
    awk -F, -v lim=$NP 'NR==1{print "pair,pregid,postgid"} NR>1 && !s[$1]++ && (lim==0 || ++c<=lim){print $1","$2","$3}' "$2" > $W/pairs/$1.csv
}
has l23l5  && mkpairs l23l5  $OLDIDX/index_Ebner2019_L23PC_L5TTPC.csv
has l23l23 && mkpairs l23l23 $OLDIDX/index_Zilberter2009_L23PC_L23PC.csv
plist() { tail -n +2 "$1" | cut -d, -f1 | paste -sd,; }
npairs() { tail -n +2 "$1" | wc -l; }

# ---------------------------------------------------------------- derived yamls: sims_dir -> scratch, circuit.config -> new emodel,
# optional protocol subset (comma ids; the Sabrina grid yaml has no protocols list: pilot narrows its dt list instead)
mkyaml() {   # src_yaml_name cfg ids out_tag
    local src=$ROOT/configs/$1.yaml out=$W/configs/$4.yaml
    sed -e "s#^sims_dir:.*#sims_dir: \"$W/refitting_results\"#" -e "s#^    config:.*#    config: \"$2\"#" "$src" \
    | { if [ "$PILOT" = 1 ] && [ "$1" = Sabrina_L5TTPC_L5TTPC_STDP ]; then sed -e 's#^    dt: \[.*#    dt: [-10.0, 10.0]  \# pilot#'; else cat; fi; } \
    | awk -v ids="$3" 'BEGIN{n=split(ids,a,","); for(i=1;i<=n;i++) keep[a[i]]=1; all=(ids=="")}
        /^protocols:/ {inp=1; print; next}
        inp && /^[^ #\t]/ {inp=0}
        inp { if ($0 ~ /^ *#/) next
              if ($0 ~ /^ *- /) { cur=0; if (all) cur=1; else if (match($0,/id: *[^,} ]+/)) { id=substr($0,RSTART,RLENGTH); sub(/id: */,"",id); if (id in keep) cur=1 } }
              if (cur) print; next }
        {print}' > $out
    grep -q "$W/refitting_results" $out || { echo "yaml $out: sims_dir not replaced"; exit 1; }
    echo $out
}

# ---------------------------------------------------------------- job submission
JOBS=$W/jobs; LOGS=$W/logs
IDN=0; EST=0; LAST=""; DRY=1
PRE="cd $ROOT; source $ROOT/.venv/bin/activate; export PYTHONPATH=$ROOT"
join() { local o=""; for x in "$@"; do [ -n "$x" ] && o="$o:$x"; done; echo "${o#:}"; }
sub() {   # sub STEP CPUS MEM TIME DEP EST_CPUH [ARRAY] [ACCOUNT] [extra sbatch args...]; script body on stdin; id in $LAST
    local step=$1 cpus=$2 mem=$3 tl=$4 dep=$5 est=$6 arr=${7:-} acct=${8:-rrg-emuller}
    local extra=("${@:9}") f=$JOBS/$1.sh
    { printf '#!/bin/bash\nset -euo pipefail\n'; cat; } > $f
    EST=$(awk -v a="$EST" -v b="$est" 'BEGIN{print a+b}')
    if [ "$DRY" = 1 ]; then IDN=$((IDN+1)); LAST=D$IDN
    else
        local a=(--parsable -J "${NAME}_$step" -A "$acct" -c "$cpus" --mem="$mem" -t "$tl" -D "$ROOT")
        if [ -n "$arr" ]; then a+=(--array="$arr" -o "$LOGS/%x_%A_%a.out"); else a+=(-o "$LOGS/%x_%j.out"); fi
        [ -n "$dep" ] && a+=(-d "afterok:$dep")
        LAST=$(sbatch "${a[@]}" "${extra[@]}" "$f")
        printf '%-26s %-10s %3s cpu %5s %s est %6.1f CPU-h %s\n' "$step" "$LAST" "$cpus" "$mem" "$tl" "$est" "${dep:+afterok:$dep}" >> $W/jobids.txt
    fi
}

simwriter() {   # step dep cpus mem time est yaml_path pairs_csv
    sub $1 $3 $4 $5 "$2" $6 <<EOF
cd $ROOT/plastyfire; source $ROOT/.venv/bin/activate; export PYTHONPATH=$ROOT
python -u simwriter.py --config $7 --pairs-from $8 --workers $3
EOF
}

prefire_body() {   # sims cache cfg hash pairs protos workers extra
    echo "python -u run_de_fit2_pool.py --cooker --param-hash $4 --results-dir $1 --cache $2 --circuit-config $3 --trace-vars cai_CR,shaft_cai,ica_VDCC,v_seg --skip-existing --pairs \$PAIRS --workers $7 --protocols $6 ${8:-}"
}
prefire() {   # step dep cpus mem time est_total sims cache cfg hash allpairs protos workers [extra]; array of 3 (pairs i::3) when full
    local arr="" split="PAIRS=${11}"
    if [ "$PILOT" != 1 ]; then arr=0-2; split="PAIRS=\$(echo ${11} | tr , '\n' | awk -v t=\$SLURM_ARRAY_TASK_ID '(NR-1)%3==t' | paste -sd,)"; fi
    sub $1 $3 $4 $5 "$2" $6 "$arr" rrg-emuller <<EOF
$PRE
$split
echo "pairs \$PAIRS"
$(prefire_body $7 $8 $9 ${10} x ${12} ${13} "${14:-}")
EOF
}

basis() {   # step dep cpus mem time est_total sims out cfg csv
    if [ "${SKIPBASIS:-0}" = 1 ]; then LAST=""; return; fi
    local n; n=$(npairs ${10})
    sub $1 $3 $4 $5 "$2" $6 "0-$((n-1))" <<EOF
$PRE
SIMS=$7; OUT=$8
PAIR=\$(ls \$SIMS | grep -E '^[0-9]+-[0-9]+\$' | sort | sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p")
[ -z "\$PAIR" ] && { echo "no pair for task \$SLURM_ARRAY_TASK_ID"; exit 0; }
PRE_G=\${PAIR%-*}; POST_G=\${PAIR#*-}; CSV=\$OUT/basis_\${PRE_G}_\${POST_G}.csv
[ -f "\$CSV" ] && { echo "exists \$CSV"; exit 0; }
CFG=""
for d in \$(ls \$SIMS/\$PAIR | sort); do
    [ -f \$SIMS/\$PAIR/\$d/simulation_config.json ] && [ -f \$SIMS/\$PAIR/\$d/prespikes.h5 ] && { CFG=\$SIMS/\$PAIR/\$d/simulation_config.json; break; }
done
[ -z "\$CFG" ] && { echo "no simulation_config.json + prespikes.h5 for \$PAIR"; exit 1; }
mkdir -p \$OUT
python -u run_basis_pair_edges.py --pre-gid \$PRE_G --post-gid \$POST_G --sim-config \$CFG --output-csv \$CSV \
    --num-trials 5 --workers 12 --circuit-config $9
EOF
}

cache() {   # step dep cpus mem time est sims cfg edges out
    sub $1 $3 $4 $5 "$2" $6 <<EOF
$PRE
python -u precompute_cpre_cpost.py --params defit2 --results-dir $W/refitting_results --sims-dir $7 \
    --circuit-config $8 --edges-h5 $9 --output ${10} --workers $3
EOF
}

extract_line() {   # sims cache hash outdir pairs protos workers [window flag; "" = no --window]
    echo "python -u $V2/extract_v2.py --param-hash $3 --sims $1 --cache $2 --protocols $6 $(w=${8:---window}; [ "$w" = NOWIN ] || echo "$w") --workers $7 --skip-existing --out $4 --pairs $5"
}

# ---------------------------------------------------------------- the chain
build() {
    : > $W/jobids.txt
    local L5_PF="" L5_BASIS="" L5_CACHE="" FIT_DEPS="" CMP_DEPS=""
    local P_L5 P_L23L5 P_L23L23
    # ===== L5 -> L5 (24 subset24 pairs) =====
    if has l5l5; then
        P_L5=$(plist $CSV_L5)
        local prev="" y s sim_last=""
        declare -A SIMJ
        # simwriters run in sequence (shared single_cells pkls). sizes: sj04 22127910 / sj07 22086838 measured; the others scaled (pilot measures)
        for s in mk sj01 sj03 sj03r50 sj07 sj04; do
            hasset $s || continue
            local calcfg=$CFG_L5L5; if [ $s = mk ]; then calcfg=${L5_MK_CAL_CFG:-$CFG_L5L5}; fi
            y=$(mkyaml ${YAML[$s]} $calcfg "$([ $s = mk ] && echo "" || ids $s)" l5_$s)
            local cpus=8 mem=2G tl=00:30:00 est=0.8
            [ $s = mk ] && { tl=00:15:00; est=0.3; }
            [ $s = sj07 ] && { tl=00:15:00; est=0.7; }
            [ $s = sj04 ] && { tl=00:15:00; est=0.3; }
            [ "$PILOT" = 1 ] && { cpus=2; mem=2G; tl=00:15:00; est=0.1; }
            simwriter sim_l5_$s "$prev" $cpus $mem $tl $est $y $CSV_L5
            prev=$LAST; SIMJ[$s]=$LAST
        done
        # c_pre/c_post cache (24 pairs; scaled from 22085683: 120 pairs 16 workers 5:11 19.8 GiB)
        if [ "$PILOT" = 1 ]; then cache cache_l5 "${SIMJ[mk]}" 2 4G 00:15:00 0.05 $SIMS_L5 $CFG_L5L5 $EDGES_L5L5 $W/cache/sabrina_n120_${NAME}.pkl
        else cache cache_l5 "${SIMJ[mk]}" 8 13G 00:15:00 0.3 $SIMS_L5 $CFG_L5L5 $EDGES_L5L5 $W/cache/sabrina_n120_${NAME}.pkl; fi
        L5_CACHE=$LAST
        # EPSP basis, one array task per pair (12 CPU = 12 trial procs; L23->L5 22086449-568: <=11.6 GiB, <=14:49, mean ~6 min)
        if [ "$PILOT" = 1 ]; then basis basis_l5 "${SIMJ[mk]}" 12 15G 00:30:00 2.4 $SIMS_L5 $W/basis_l5l5 $CFG_L5L5 $CSV_L5
        else basis basis_l5 "${SIMJ[mk]}" 12 15G 00:30:00 29 $SIMS_L5 $W/basis_l5l5 $CFG_L5L5 $CSV_L5; fi
        L5_BASIS=$LAST
        # prefire groups, sized as the measured og-delta equivalents
        local pfdeps=""
        local A_IDS="" B_IDS=""
        hasset mk   && A_IDS=$(ids mk)
        hasset sj01 && A_IDS=${A_IDS:+$A_IDS,}$(ids sj01)
        hasset sj03 && B_IDS=$(ids sj03)
        hasset sj03r50 && B_IDS=${B_IDS:+$B_IDS,}$(ids sj03r50)
        local lastsimA=${SIMJ[sj01]:-${SIMJ[mk]}} lastsimB=${SIMJ[sj03r50]:-${SIMJ[sj03]:-${SIMJ[sj01]:-${SIMJ[mk]}}}}
        # A: Markram + Sjostrom 2001 (22093947: 3 tasks x 16 workers, 103-110 GiB, 20:53-27:45 -> 140G, 0:45)
        if [ "$PILOT" = 1 ]; then prefire pf_l5_A "$(join $L5_CACHE $lastsimA)" 8 41G 00:15:00 0.7 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$A_IDS" 8
        else prefire pf_l5_A "$(join $L5_CACHE $lastsimA)" 16 140G 00:45:00 19 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$A_IDS" 16; fi
        pfdeps=$LAST
        # B: Sjostrom 2003 (9) + r50 (2) (22093946: 32 workers, 131-164 GiB, 19:57-27:50 -> 205G, 0:45)
        if [ -n "$B_IDS" ]; then
            local nb; nb=$(echo $B_IDS | tr , '\n' | wc -l)
            if [ "$PILOT" = 1 ]; then prefire pf_l5_B "$(join $L5_CACHE $lastsimB)" $((2*nb)) $((2*nb*5+1))G 00:30:00 $(awk -v n=$nb 'BEGIN{print 2*n*0.2}') $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$B_IDS" $((2*nb))
            else prefire pf_l5_B "$(join $L5_CACHE $lastsimB)" 32 205G 00:45:00 37.6 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$B_IDS" 32; fi
            pfdeps=$(join $pfdeps $LAST)
        fi
        # C: Sjostrom 2007 (22089771: 16 workers, 78.99 GB at the limit, 21:44 -> 98G, 0:45)
        if hasset sj07; then
            if [ "$PILOT" = 1 ]; then prefire pf_l5_sj07 "$(join $L5_CACHE ${SIMJ[sj07]})" 4 21G 00:30:00 0.6 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$(ids sj07)" 4
            else prefire pf_l5_sj07 "$(join $L5_CACHE ${SIMJ[sj07]})" 16 98G 00:45:00 5.8 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$(ids sj07)" 16; fi
            pfdeps=$(join $pfdeps $LAST)
        fi
        # D: Sjostrom 2004 dLTD (22128324: 12 workers, 56.2 GB, 11:41 -> 70G, 0:30; pilot 22127911 2 workers 9.99 GB -> 13G, 0:15)
        if hasset sj04; then
            if [ "$PILOT" = 1 ]; then prefire pf_l5_sj04 "$(join $L5_CACHE ${SIMJ[sj04]})" 2 13G 00:15:00 0.2 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$(ids sj04)" 2
            else prefire pf_l5_sj04 "$(join $L5_CACHE ${SIMJ[sj04]})" 12 70G 00:30:00 2.3 $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $CFG_L5L5 $H5 "$P_L5" "$(ids sj04)" 12; fi
            pfdeps=$(join $pfdeps $LAST)
        fi
        # extraction, all L5 datasets in one job (22126528: 4 workers, 121 files 1:51, 78.9 GB; markram 115 GB peak 22074074 -> 145G)
        local EX="export PLASTYFIRE_EDGES=$EDGES_L5L5"
        for s in mk sj01 sj03 sj03r50 sj07 sj04; do
            hasset $s || continue
            EX="$EX
$(extract_line $SIMS_L5 $W/cache/sabrina_n120_${NAME}.pkl $H5 $X/${DIRN[$s]}_${HP}-vca "$P_L5" "$(ids $s)" $([ "$PILOT" = 1 ] && echo 2 || echo 4) $([ $s = mk ] && echo NOWIN || echo --window))"
        done
        if [ "$PILOT" = 1 ]; then sub ext_l5 2 60G 00:15:00 "$(join $pfdeps $L5_CACHE)" 0.4 <<EOF
$PRE
$EX
EOF
        else sub ext_l5 4 145G 00:30:00 "$(join $pfdeps $L5_CACHE)" 1.5 <<EOF
$PRE
$EX
EOF
        fi
        local ext_l5=$LAST
        FIT_DEPS=$(join $FIT_DEPS $ext_l5 $L5_BASIS)
        # compare with the existing og-delta extraction (only meaningful when EMODEL is og-delta, i.e. a test name ogtest*)
        if [[ $NAME == ogtest* ]]; then
            local CP="" s
            for s in mk sj01 sj03 sj03r50 sj07 sj04; do hasset $s && CP="$CP --pair $X/${DIRN[$s]}_${HP}-vca=$X/${DIRN[$s]}_delta-prefire-vca"; done
            sub cmp_l5 1 2G 00:15:00 "$(join $ext_l5 $L5_BASIS)" 0.05 <<EOF
$PRE
python -u $V2/compare_emodel_npz.py $CP $([ "${SKIPBASIS:-0}" = 1 ] || echo "--basis $W/basis_l5l5=$ROOT/basis_results_edges_sabrina_n120_delta")
EOF
        fi
    fi
    # ===== L2/3 -> L5 (120 pairs; + Letzkus null controls) =====
    if has l23l5; then
        P_L23L5=$(plist $W/pairs/l23l5.csv)
        local idsm="" idsp=""
        [ "$PILOT" = 1 ] && idsm=$L23L5_IDS_PILOT
        local ym; ym=$(mkyaml Ebner2019_L23PC_L5TTPC $CFG_L5L5 "$idsm" l23l5)
        # simwriter with fresh calibration: 22043887 17:56, 36.4 GB on 32 workers (incl. pair search); 22128319 (L2/3->L2/3) 18:58, 19.8 GB on 16
        if [ "$PILOT" = 1 ]; then simwriter sim_l23l5 "" 2 2G 00:15:00 0.1 $ym $W/pairs/l23l5.csv
        else simwriter sim_l23l5 "" 16 25G 01:00:00 8 $ym $W/pairs/l23l5.csv; fi
        local sim_l23l5=$LAST sim_null=$LAST
        if [ "$NULLS" = 1 ]; then
            local yn; yn=$(mkyaml Ebner2019_L23PC_L5TTPC_nulls $CFG_L5L5 "$NULL_IDS" l23l5_nulls)
            simwriter sim_l23l5_nulls "$sim_l23l5" 1 1G 00:15:00 0.01 $yn $W/pairs/l23l5.csv   # 22123227 (cached stimuli): 20 s, 315 MB
            sim_null=$LAST
        fi
        if [ "$PILOT" = 1 ]; then cache cache_l23l5 "$sim_l23l5" 2 4G 00:15:00 0.05 $SIMS_L23L5 $CFG_L23L5 $EDGES_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl
        else cache cache_l23l5 "$sim_l23l5" 16 25G 00:15:00 1.4 $SIMS_L23L5 $CFG_L23L5 $EDGES_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl; fi   # 22085683: 5:11, 19.8 GiB
        local cache_l23l5=$LAST
        if [ "$PILOT" = 1 ]; then basis basis_l23l5 "$sim_l23l5" 12 15G 00:30:00 2.4 $SIMS_L23L5 $W/basis_l23l5 $CFG_L23L5 $W/pairs/l23l5.csv
        else basis basis_l23l5 "$sim_l23l5" 12 15G 00:30:00 170 $SIMS_L23L5 $W/basis_l23l5 $CFG_L23L5 $W/pairs/l23l5.csv; fi
        local basis_l23l5=$LAST
        local pfids; pfids=$( [ "$PILOT" = 1 ] && echo $L23L5_IDS_PILOT || echo "sjostrom_50hz_dt+10ms,letzkus_1ap_dt+10ms,letzkus_3ap_200hz_dt+10ms,letzkus_3ap_200hz_dt-10ms")
        # prefire main (22089155: 32 workers, 189.7 GiB, 37:18 -> 255G, 1:00); single job (no array)
        local pfm pfn=""
        if [ "$PILOT" = 1 ]; then
            sub pf_l23l5 4 21G 00:15:00 "$(join $cache_l23l5 $sim_l23l5)" 0.3 <<EOF
$PRE
PAIRS=$P_L23L5
$(prefire_body $SIMS_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl $CFG_L23L5 $H23 x $pfids 4)
EOF
        else
            sub pf_l23l5 32 255G 01:00:00 "$(join $cache_l23l5 $sim_l23l5)" 20 <<EOF
$PRE
PAIRS=$P_L23L5
$(prefire_body $SIMS_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl $CFG_L23L5 $H23 x $pfids 32)
EOF
        fi
        pfm=$LAST
        local exdeps=$pfm
        if [ "$NULLS" = 1 ]; then   # 22123818: 16 workers, 125 GB at the limit, 39:18 -> 150G, 1:00
            sub pf_l23l5_nulls 16 150G 01:00:00 "$(join $cache_l23l5 $sim_null)" 10.5 <<EOF
$PRE
PAIRS=$P_L23L5
$(prefire_body $SIMS_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl $CFG_L23L5 $HN x $NULL_IDS 16)
EOF
            pfn=$LAST; exdeps=$(join $pfm $pfn)
        fi
        local OUT23=$X/ebner_l23l5_${HP}-vseg-rs
        local EX23="export PLASTYFIRE_EDGES=$EDGES_L23L5
$(extract_line $SIMS_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl $H23 $OUT23 $P_L23L5 $pfids $([ "$PILOT" = 1 ] && echo 2 || echo 4))"
        [ "$NULLS" = 1 ] && EX23="$EX23
$(extract_line $SIMS_L23L5 $W/cache/ebner_l23l5_${NAME}_rs.pkl $HN $OUT23 $P_L23L5 $NULL_IDS 4)"
        if [ "$PILOT" = 1 ]; then sub ext_l23l5 2 50G 00:15:00 "$(join $exdeps $cache_l23l5)" 0.2 <<EOF
$PRE
$EX23
EOF
        else sub ext_l23l5 4 99G 00:15:00 "$(join $exdeps $cache_l23l5)" 0.6 <<EOF
$PRE
$EX23
EOF
        fi
        local ext_l23l5=$LAST
        FIT_DEPS=$(join $FIT_DEPS $ext_l23l5 $basis_l23l5)
        if [[ $NAME == ogtest* ]]; then
            sub cmp_l23l5 1 2G 00:15:00 "$(join $ext_l23l5 $basis_l23l5)" 0.05 <<EOF
$PRE
python -u $V2/compare_emodel_npz.py --pair $OUT23=$X/ebner_l23l5_delta-prefire-vseg-rs $([ "${SKIPBASIS:-0}" = 1 ] || echo "--basis $W/basis_l23l5=$ROOT/basis_results_edges_ebner_l23l5_delta_rs")
EOF
        fi
    fi
    # ===== L2/3 -> L2/3 (120 pairs, 13 Zilberter protocols) =====
    if has l23l23; then
        P_L23L23=$(plist $W/pairs/l23l23.csv)
        local idz="" ; [ "$PILOT" = 1 ] && idz=$L23L23_IDS_PILOT
        local yz; yz=$(mkyaml Zilberter2009_L23PC_L23PC $CFG_L5L5 "$idz" l23l23)
        local zpf; zpf=$( [ "$PILOT" = 1 ] && echo $L23L23_IDS_PILOT || echo "zilberter_1ap_dt+10ms,zilberter_1ap_dt-10ms,zilberter_pre_only,zilberter_5ap_10hz_dt+10ms,zilberter_5ap_20hz_dt+10ms,zilberter_5ap_20hz_dt-10ms,zilberter_train10_50hz_dt+4ms_last,zilberter_train10_50hz_dt-4ms_last,zilberter_train10_50hz_dt-10ms_last,zilberter_train10_50hz_dt+5ms,zilberter_train10_50hz_post_only,zilberter_train4_50hz_dt+4ms_last,zilberter_train8_50hz_dt+4ms_last")
        # 22128319: 18:58, 19.8 GB on 16 -> 25G, 0:30
        if [ "$PILOT" = 1 ]; then simwriter sim_l23l23 "" 2 2G 00:15:00 0.1 $yz $W/pairs/l23l23.csv
        else simwriter sim_l23l23 "" 16 25G 00:30:00 5.1 $yz $W/pairs/l23l23.csv; fi
        local sim_z=$LAST
        if [ "$PILOT" = 1 ]; then cache cache_l23l23 "$sim_z" 2 4G 00:15:00 0.05 $SIMS_L23L23 $CFG_L23L23 $EDGES_L23L23 $W/cache/zilberter_l23l23_${NAME}_rs.pkl
        else cache cache_l23l23 "$sim_z" 16 13G 00:15:00 0.7 $SIMS_L23L23 $CFG_L23L23 $EDGES_L23L23 $W/cache/zilberter_l23l23_${NAME}_rs.pkl; fi   # 22128321: 2:30, 10.2 GB
        local cache_z=$LAST
        if [ "$PILOT" = 1 ]; then basis basis_l23l23 "$sim_z" 12 7G 00:15:00 0.5 $SIMS_L23L23 $W/basis_l23l23 $CFG_L23L23 $W/pairs/l23l23.csv
        else basis basis_l23l23 "$sim_z" 12 7G 00:15:00 24 $SIMS_L23L23 $W/basis_l23l23 $CFG_L23L23 $W/pairs/l23l23.csv; fi   # 22128322: max 1:39, 5.33 GB
        local basis_z=$LAST
        if [ "$PILOT" = 1 ]; then
            sub pf_l23l23 4 21G 00:15:00 "$(join $cache_z $sim_z)" 0.3 <<EOF
$PRE
PAIRS=$P_L23L23
$(prefire_body $SIMS_L23L23 $W/cache/zilberter_l23l23_${NAME}_rs.pkl $CFG_L23L23 $H23 x $zpf 4)
EOF
        else   # 22129908: 32 workers, 248 GB at the limit, 1:13:57 -> 310G, 1:45
            sub pf_l23l23 32 310G 01:45:00 "$(join $cache_z $sim_z)" 39.4 <<EOF
$PRE
PAIRS=$P_L23L23
$(prefire_body $SIMS_L23L23 $W/cache/zilberter_l23l23_${NAME}_rs.pkl $CFG_L23L23 $H23 x $zpf 32)
EOF
        fi
        local pf_z=$LAST
        local OUTZ=$X/zilberter_l23l23_${HP}-vseg-rs
        local EXZ="export PLASTYFIRE_EDGES=$EDGES_L23L23
$(extract_line $SIMS_L23L23 $W/cache/zilberter_l23l23_${NAME}_rs.pkl $H23 $OUTZ $P_L23L23 $zpf $([ "$PILOT" = 1 ] && echo 2 || echo 4))"
        if [ "$PILOT" = 1 ]; then sub ext_l23l23 2 50G 00:15:00 "$(join $pf_z $cache_z)" 0.2 <<EOF
$PRE
$EXZ
EOF
        else sub ext_l23l23 4 124G 00:45:00 "$(join $pf_z $cache_z)" 1.4 <<EOF
$PRE
$EXZ
EOF
        fi
        local ext_z=$LAST
        FIT_DEPS=$(join $FIT_DEPS $ext_z $basis_z)
        if [[ $NAME == ogtest* ]]; then
            sub cmp_l23l23 1 2G 00:15:00 "$(join $ext_z $basis_z)" 0.05 <<EOF
$PRE
python -u $V2/compare_emodel_npz.py --pair $OUTZ=$X/zilberter_l23l23_delta-prefire-vseg-rs $([ "${SKIPBASIS:-0}" = 1 ] || echo "--basis $W/basis_l23l23=$ROOT/basis_results_edges_zilberter_l23l23_delta_rs")
EOF
        fi
    fi
    # ===== joint fit (best settings: C1Ajn_s5 recipe, 3 pathways; seeded from the og-delta best fit + an unseeded basin check) =====
    if [ "$SKIPFIT" != 1 ]; then
        has l5l5 && has l23l5 && has l23l23 || { echo "joint fit needs PATHWAYS=l5l5,l23l5,l23l23; skipping fit"; return 0; }
        local L5DIRS="" s G=paired_l5,sjostrom07
        for s in mk sj01 sj03 sj03r50 sj07 sj04; do hasset $s && L5DIRS="$L5DIRS,$X/${DIRN[$s]}_${HP}-vca"; done
        L5DIRS=${L5DIRS#,}; hasset sj04 && G=$G,sjostrom04
        local EXTRA="l23:paired_l23l5:$OUT23:$W/basis_l23l5:$ROOT/ebner/pair_geometry_L23PC_L5TTPC.csv;l23l23:paired_l23l23:$OUTZ:$W/basis_l23l23"
        local seedfits=${FIT_SEEDFITS:-$V2/rho_redesign/results/v4_C1Ajn_s5.json}
        local tagn seed seedenv fset='{"vamp_mode": 1}' ffree=theta_V,rho_gamma fdrop='letzkus_3ap_200hz_dt-10ms@distal|control,letzkus_3ap_200hz_dt-10ms@distal|nmdar_block' fixline=""
        fset=${FIT_SET:-$fset}; ffree=${FIT_FREE:-$ffree}; fdrop=${FIT_DROPT:-$fdrop}
        [ -n "${FIT_FIX:-}" ] && fixline="export FIX='$FIT_FIX'"
        for tagn in s5:5:1 s6:6:0; do
            IFS=: read -r tg seed seeded <<< "$tagn"
            seedenv=""; [ "$seeded" = 1 ] && seedenv="export SEEDFITS='$seedfits'; export SEEDSET='{}'"
            # 3 pathways: 22134473 MAXITER=2 pilot 6:01, 63.28 GB on 3g.40gb; full fit 25 min x 2.5 records + 50% -> 79G, 1:45
            sub fit_$tg 1 79G 01:45:00 "$FIT_DEPS" 1 "" def-emuller --gpus-per-node=nvidia_h100_80gb_hbm3_3g.40gb:1 <<EOF
cd $ROOT
export TAG=${NAME}_j$tg SEED=$seed FITGAMMA=1 FITPY=fit_v4n.py MAXITER=${FIT_MAXITER:-150}
export SET='$fset' FREE='$ffree'
export DROPT='$fdrop'
$fixline
export L5DIRS='$L5DIRS' L5GROUPS='$G' EXTRA='$EXTRA' ANALYTICAL_BASIS_DIR=$W/basis_l5l5
$seedenv
bash $V2/rho_redesign/run_fit_v4.sh
EOF
        done
    fi
}

# phase 1: plan (writes the job scripts, submits nothing) -> estimate; phase 2: submit
DRY=1; build
echo "plan: name $NAME, pathways $PATHWAYS, pilot $PILOT, L5 sets $L5_SETS, estimated ${EST} CPU-h (jobs in $JOBS)"
if [ "$DRYRUN" = 1 ]; then ls $JOBS; exit 0; fi
if [ "$PILOT" != 1 ] && awk -v e="$EST" 'BEGIN{exit !(e>30)}' && [ "${DECISIONS_LOGGED:-0}" != 1 ]; then
    echo "estimated $EST CPU-h > 30: log the run in glusynapse_v2/DECISIONS.md first (CLAUDE.md rule 7), then set DECISIONS_LOGGED=1"; exit 1; fi
DRY=0; EST=0; build
echo "submitted (ids in $W/jobids.txt):"; cat $W/jobids.txt
echo "check: jobstat; npz in $X/*_${HP}-*; compare logs $LOGS/${NAME}_cmp_*.out"
