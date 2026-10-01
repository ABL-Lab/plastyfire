#!/bin/bash
#SBATCH --job-name=quota_scan
#SBATCH --account=rrg-emuller
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0:15:00
#SBATCH --output=/scratch/dhuruva/quota_scan/%x_%j.out
# Inode audit scan of /project/rrg-emuller/dhuruva (QUOTA_AUDIT.md).
# PILOT: no comparable finished job. Sizing basis: ~335K inodes; one find pass
# plus awk/sort on a ~40 MB listing. Est. well under 1 GB RSS and <30 min.
# MEASURED (pilot 22127461, 2026-10-01, 2G/1:00): elapsed 2:16, MaxRSS 818 MB,
# CPU eff 29% (I/O bound on find), 333,957 inodes. Next runs: 1G, 0:15 (set above).
# Read-only: lists and counts, never modifies /project. All output goes to /scratch.
# The rejected ljp25 hash (being offloaded by job 22127439) is excluded from the
# aggregates and only counted on its own line.
set -euo pipefail
D=/project/rrg-emuller/dhuruva
O=/scratch/dhuruva/quota_scan
SKIP=delta-ljp25-prefire-vseg
mkdir -p "$O"
cd "$D"

# 1) one filesystem pass: type, size, relative path
find . -xdev -mindepth 1 -printf '%y\t%s\t%P\n' > "$O/listing.tsv"
echo "total inodes: $(wc -l < "$O/listing.tsv")"
echo "ljp25 inodes (excluded below): $(grep -c "$SKIP" "$O/listing.tsv" || true)"
grep -v "$SKIP" "$O/listing.tsv" > "$O/listing_noljp.tsv"
L="$O/listing_noljp.tsv"

# 2) subtree inode counts per directory prefix, depths 1-5 (>=200 inodes)
awk -F'\t' '{n=split($3,a,"/"); p=""; for(i=1;i<=n-1 && i<=5;i++){p=(i==1)?a[1]:p"/"a[i]; c[i"\t"p]++}}
     END{for(k in c) if(c[k]>=200) print c[k]"\t"k}' "$L" | sort -t$'\t' -k1,1nr > "$O/dirs_depth1-5.tsv"

# 3) per dataset x prefire hash: bluecellulab_results_<hash>/ contents + pool_<hash>_*.log
#    + simulation_edges_<hash>.pkl / simulation_<hash>.pkl, with bytes
awk -F'\t' '$3 ~ /^plastyfire\/refitting_results\// {
     n=split($3,a,"/"); ds=a[6]; h=""; f=a[n]
     for(i=1;i<=n;i++) if(a[i] ~ /^bluecellulab_results_/){h=substr(a[i],22); break}
     if(h=="" && f ~ /^pool_/){ m=split(substr(f,6),t,"_"); s=(t[1]=="full")?"full_"t[2]:t[1]; h="pool:"s }   # hashes have no "_" except the full_ prefix
     if(h=="" && f ~ /^simulation_(edges_)?.*\.pkl$/){ s=f; sub(/^simulation_(edges_)?/,"",s); sub(/\.pkl$/,"",s); h="pkl:"s }
     if(h=="") h="(shared)"
     c[ds"\t"h]++; b[ds"\t"h]+=$2 }
     END{for(k in c) printf "%d\t%.1f\t%s\n", c[k], b[k]/1e9, k}' "$L" \
  | sort -t$'\t' -k3,3 -k1,1nr > "$O/hash.tsv"   # cols: inodes GB dataset hash

# 4) per top-level project x extension (files only)
awk -F'\t' '$1=="f"{n=split($3,a,"/"); f=a[n]; e=(f ~ /\./)?f:"<noext>"; sub(/^.*\./,".",e); c[a[1]"\t"e]++}
     END{for(k in c) if(c[k]>=100) print c[k]"\t"k}' "$L" | sort -t$'\t' -k1,1nr > "$O/ext.tsv"

# 5) caches / build dirs / tiny-file classes per top-level project
awk -F'\t' '{p=$3; split(p,a,"/"); k=""
     if(p ~ /(^|\/)\.git\//) k="git"
     else if(p ~ /(^|\/)\.?venv[^\/]*\//) k=(p ~ /__pycache__/)?"venv_pycache":"venv"
     else if(p ~ /__pycache__/) k="pycache"
     else if(p ~ /\.ipynb_checkpoints/) k="ipynb_checkpoints"
     else if(p ~ /(^|\/)x86_64(\/|$)/) k="nrn_x86_64"
     else if(p ~ /node_modules/) k="node_modules"
     else if(p ~ /(^|\/)(\.cache|pip-cache|\.conda|conda-pkgs)(\/|$)/) k="cache"
     else if(p ~ /\.SUCCESS$/) k="SUCCESS"
     else if(p ~ /\.(batch|sh)$/) k="sh_batch"
     else if(p ~ /\.log$/ || p ~ /(^|\/)logs?\//) k="logs"
     if(k!="") c[a[1]"\t"k]++}
     END{for(k in c) print c[k]"\t"k}' "$L" | sort -t$'\t' -k1,1nr > "$O/caches.tsv"

# 6) empty dirs named out/ under refitting_results; neurodamus delta-cooker artefacts
awk -F'\t' '{n=split($3,a,"/"); par=a[1]; for(i=2;i<n;i++) par=par"/"a[i]; if(n>1) haskid[par]=1
     if($1=="d" && a[n]=="out" && $3 ~ /^plastyfire\/refitting_results\//) outd[$3]=1}
     END{e=0; for(d in outd) if(!(d in haskid)) e++; print e"\tempty out/ dirs in refitting_results"}' "$L" > "$O/misc.tsv"
awk -F'\t' '$3 ~ /^plastyfire\/refitting_results\// && $3 !~ /bluecellulab_results_/ {n=split($3,a,"/"); f=a[n]
     if(f ~ /^(out\.h5|soma\.h5|rho\.h5|mcomplex\.dat|populations_offset\.dat|cx_S1nonbarrel.*\.dat|prefire_simulation_config_delta-cooker\.json(\.SUCCESS)?|pydamus_.*\.log)$/){c++; b+=$2}}
     END{printf "%d\t%.1f GB\tneurodamus delta-cooker outputs (excl. pool/pkl)\n", c, b/1e9}' "$L" >> "$O/misc.tsv"

rm -f "$O/listing_noljp.tsv"   # keep listing.tsv for follow-up queries
ls -l "$O"
echo "elapsed ${SECONDS}s"
