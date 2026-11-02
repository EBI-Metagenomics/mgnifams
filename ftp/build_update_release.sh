#!/usr/bin/env bash
#SBATCH --job-name=mgnifams_ftp_stage
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
# Stage the files of an update release for the FTP (README.md, "Publishing an update release", step 2), outside the
# public tree. Hard links to the previous release's HMMs and seeds are made later, inside the FTP tree.
#
#   sbatch build_update_release.sh <update outdir> <previous release dir on the FTP> <stage dir>
#
# Writes <stage dir>/{full_msa.tar,families.tsv.gz,pipeline_info/,md5sums.txt,sizes.txt}. Read-only on its inputs.
set -euo pipefail
OUT=$1 PREV=$2 ST=$3
FAM=$OUT/update_families
mkdir -p "$ST/prev" "$ST/all/full_msa" "$ST/pipeline_info"

total=$(zcat "$FAM/families.tsv.gz" | tail -n +2 | wc -l)
echo "families: $total"

# 1. The previous alignments of the families this update did not refresh
zcat "$FAM/families.tsv.gz" | awk -F'\t' 'NR>1 && $2=="not_updated"{print "full_msa/"$1".sto.gz"}' > "$ST/not_updated.txt"
if [ -s "$ST/not_updated.txt" ]; then
    tar xf "$PREV"/full_msa.tar* -C "$ST/prev" -T "$ST/not_updated.txt"
fi
kept=$(find "$ST/prev" -name '*.sto.gz' | wc -l)
echo "kept from the previous release: $kept of $(wc -l < "$ST/not_updated.txt")"
[ "$kept" -eq "$(wc -l < "$ST/not_updated.txt")" ]

# 2. full_msa.tar: one member per family, sorted, without owner names (symlinks followed with -h)
find "$FAM/full_msa" "$ST/prev/full_msa" -name '*.sto.gz' -exec ln -s -t "$ST/all/full_msa" {} + 2>/dev/null || true
linked=$(find "$ST/all/full_msa" -name '*.sto.gz' | wc -l)
echo "alignments linked: $linked"
[ "$linked" -eq "$total" ]
tar cf "$ST/full_msa.tar" -h --sort=name --owner=0 --group=0 --numeric-owner -C "$ST/all" full_msa
members=$(tar tf "$ST/full_msa.tar" | grep -c '\.sto\.gz$')
echo "full_msa.tar members: $members"
[ "$members" -eq "$total" ]

# 3. Metadata: families.tsv.gz, software versions, the last run's parameters without local paths
cp "$FAM/families.tsv.gz" "$ST/"
cp "$OUT"/pipeline_info/*versions*.yml "$ST/pipeline_info/"
params=$(ls "$OUT"/pipeline_info/params_*.json | sort | tail -1)
python3 -I - "$params" "$ST/pipeline_info/$(basename "$params")" <<'PY'
import json, sys
params = json.load(open(sys.argv[1]))
local = {k for k, v in params.items() if isinstance(v, str) and "/" in v and not v.startswith(("http://", "https://"))}
print("dropped local paths:", ", ".join(sorted(local)))
json.dump({k: v for k, v in params.items() if k not in local}, open(sys.argv[2], "w"), indent=4)
PY

# 4. Checksums: the new files, plus the unchanged HMMs and seeds from the previous release
cd "$ST"
md5sum families.tsv.gz full_msa.tar > md5sums.txt
grep -E ' (mgnifams_hmm\.lib\.gz|seed_msa\.tar\.gz)$' "$PREV/md5sums.txt" >> md5sums.txt
stat -c '%s %n' families.tsv.gz full_msa.tar > sizes.txt
cat md5sums.txt sizes.txt
rm -r "$ST/all"
echo "done"
