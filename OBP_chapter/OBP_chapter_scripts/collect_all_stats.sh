#!/bin/bash
# collect_all_stats.sh
#
# Walks every stage of the snpAD pipeline, parses the bcftools stats PSC
# lines it already wrote at each stage, and produces ONE tidy long-format
# table (Sample x Stage), ready to load into Excel or R.
#
# Run this once, after the pipeline has finished (i.e. after 6_concat).
# No bcftools/stats regeneration happens here - it only reads files that
# already exist.
#
# Usage:
#   bash collect_all_stats.sh
#
# Output (all kept - nothing temporary):
#   ${BASE}/OBP_perscaffold_stats.tsv     - one row per Sample x Scaffold x Stage
#   ${BASE}/OBP_corrected_vcf_het.tsv     - one row per Sample, corrected VCF-het
#   ${BASE}/OBP_all_stages_summary.tsv    - final tidy table (the two above combined + collapsed)

set -euo pipefail
shopt -s nullglob

BASE="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output"
MODERN_PREFIX="O-"
HIST_PREFIX="OBP"
PERSCAFFOLD="${BASE}/OBP_perscaffold_stats.tsv"
CORRECTED="${BASE}/OBP_corrected_vcf_het.tsv"
OUT="${BASE}/OBP_all_stages_summary.tsv"

echo -e "Stage\tScaffold\tFile\tSample\tCallable_sites\tRefHom\tNonRefHom\tHets\tMissing\tTransitions\tTransversions\tTsTv\tIndels\tAverageDepth\tSingletons\tHapRef\tHapAlt" > "$PERSCAFFOLD"

# ---------------------------------------------------------------------------
# Parse every .stats file matching a glob into the shared per-scaffold table.
# Scaffold is pulled out of the path itself (works whether it's a directory
# name, as in stages 1-2, or embedded in the filename, as in stages 3-6),
# so this one function covers every stage.
# ---------------------------------------------------------------------------
collect_stage () {
    local stage_label="$1"
    local glob="$2"

    for f in $glob; do
        [[ -s "$f" ]] || continue
        scaffold=$(grep -oE 'scaffold_[0-9]+' <<< "$f" | head -n1)
        [[ -n "$scaffold" ]] || scaffold="genome_wide"
        file=$(basename "$f")

        awk -F'\t' -v stage="$stage_label" -v scaffold="$scaffold" -v file="$file" '
        /^PSC/ {
            total = $4 + $5 + $6
            tstv = ($8 > 0) ? $7/$8 : "NA"
            print stage"\t"scaffold"\t"file"\t"$3"\t"total"\t"$4"\t"$5"\t"$6"\t"$14"\t"$7"\t"$8"\t"tstv"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13
        }' "$f" >> "$PERSCAFFOLD"
    done
}

echo "Collecting per-scaffold stats..."
collect_stage "1_Raw_gVCF"               "${BASE}/1_VCFs/*/*.stats"
echo "Finished 1_Raw_gVCF"
collect_stage "2_Filtered_gVCF"          "${BASE}/2_FiltGVCFs/*/*.stats"
echo "Finished 2_Filtered_gVCF"
collect_stage "3_Merged"                 "${BASE}/3_MergeGVCFs_perChrm/*.filt.g.vcf.stats"
echo "Finished 3_Merged"
collect_stage "4_GroupMiss80"            "${BASE}/4_GroupMissFiltered/*.groupmiss80.vcf.stats"
echo "Finished 4_GroupMiss80"
collect_stage "5_Bisnp_with_singletons"  "${BASE}/5_bisnps/miss80/*.groupmiss80.bisnps.vcf.stats"
echo "Finished 5_Bisnp_with_singletons"
collect_stage "6_Bisnp_no_singletons"    "${BASE}/5_bisnps/miss80/*.groupmiss80.bisnps.MAC2.vcf.stats"
echo "Finished 6_Bisnp_no_singletons"
collect_stage "7_Concat_genome_wide"     "${BASE}/6_concat/*.concat.groupmiss80.bisnps.vcf.stats"
echo "Finished 7_Concat_genome_wide"

n_rows=$(($(wc -l < "$PERSCAFFOLD") - 1))
echo "Parsed ${n_rows} per-scaffold-per-sample rows -> ${PERSCAFFOLD}"

# Files get appended in whatever order the shell glob returns them, which is
# alphabetical (scaffold_100 sorts before scaffold_2). Re-sort numerically by
# scaffold (sort -V handles the embedded number correctly), keeping the header
# on top. Writes back to the same filename - no file left behind.
echo "Sorting per-scaffold table numerically by scaffold..."
{
    head -n1 "$PERSCAFFOLD"
    tail -n +2 "$PERSCAFFOLD" | sort -t$'\t' -k1,1 -k2,2V -k4,4
} > "${PERSCAFFOLD}.sorted"
mv "${PERSCAFFOLD}.sorted" "$PERSCAFFOLD"

# ---------------------------------------------------------------------------
# Collapse to one row per Sample x Stage (sum across scaffolds; the
# already-genome-wide stage just collapses to itself).
# ---------------------------------------------------------------------------
echo "Summarising per sample per stage..."

awk -F'\t' -v modern="$MODERN_PREFIX" -v hist="$HIST_PREFIX" '
BEGIN { OFS = "\t" }
NR == 1 { next }
{
    key = $1 SUBSEP $4   # Stage, Sample
    n_scaffolds[key]++
    callable[key]      += $5
    refhom[key]        += $6
    nonrefhom[key]     += $7
    hets[key]          += $8
    missing[key]       += $9
    transitions[key]   += $10
    transversions[key] += $11
    indels[key]        += $13
    depth_wsum[key]    += $14 * $5
    singletons[key]    += $15
}
END {
    print "Sample","Group","Stage","nScaffolds","Total_callable_sites",
          "Total_RefHom","Total_NonRefHom","Total_Hets","Total_missing",
          "Total_transitions","Total_transversions","TsTv","Total_indels",
          "Mean_depth","Total_singletons","Variant_sites","Het_per_callable"

    for (k in callable) {
        split(k, a, SUBSEP)
        stage = a[1]; sample = a[2]

        group = (sample ~ "^"modern) ? "Modern" : (sample ~ "^"hist) ? "Historical" : "Unknown"

        tstv = (transversions[k] > 0) ? transitions[k]/transversions[k] : "NA"
        het_per_callable = (callable[k] > 0) ? hets[k]/callable[k] : "NA"
        mean_depth = (callable[k] > 0) ? depth_wsum[k]/callable[k] : "NA"

        print sample, group, stage, n_scaffolds[k], callable[k], refhom[k],
              nonrefhom[k], hets[k], missing[k], transitions[k],
              transversions[k], tstv, indels[k], mean_depth, singletons[k],
              nonrefhom[k]+hets[k], het_per_callable
    }
}
' "$PERSCAFFOLD" | { 
    read -r header
    printf "%s\n" "$header"
    # Custom group order (Modern before Historical) instead of alphabetical -
    # prepend Stage/rank/Sample as sort keys, sort, then drop them again.
    awk -F'\t' 'BEGIN{OFS="\t"; rank["Modern"]=0; rank["Historical"]=1} { r=($2 in rank)?rank[$2]:2; print $3, r, $1, $0 }' \
        | sort -t$'\t' -k1,1 -k2,2n -k3,3 \
        | cut -f4-
} > "$OUT"
