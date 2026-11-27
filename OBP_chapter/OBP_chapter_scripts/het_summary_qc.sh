
# To get summary het values across all individuals, run the following code after all jobs are done:

SUMMARY="heterozygosity_summary.tsv"
:> $SUMMARY

echo -e "sample\tn_boot\tmean_H\tsd_H\tmean_sites\tsd_sites" > "${SUMMARY}"

shopt -s nullglob

for file in *ml; do
  sample=$(basename "${file}" _filter_wim.ml)
  
  grep -v "SHAPE" "${file}" | awk -v sample="${sample}" '
    {
      tot = $1 + $2 + $3
      if (tot > 0) {
        h = $2 / tot
        n++

        sumH   += h
        sumsqH += h*h

        sumT   += tot
        sumsqT += tot*tot
      }
    }
    END {
      if (n == 0) {
        printf("%s\t0\tNA\tNA\tNA\tNA\n", sample)
        exit
      }

      meanH = sumH / n
      if (n > 1) {
        varH = (sumsqH - n*meanH*meanH) / (n-1)
        if (varH < 0) varH = 0
        sdH = sqrt(varH)
      } else {
        sdH = 0
      }

      meanT = sumT / n
      if (n > 1) {
        varT = (sumsqT - n*meanT*meanT) / (n-1)
        if (varT < 0) varT = 0
        sdT = sqrt(varT)
      } else {
        sdT = 0
      }

      printf("%s\t%d\t%.10f\t%.10f\t%.0f\t%.0f\n",
             sample, n, meanH, sdH, meanT, sdT)
    }' >> "${SUMMARY}"
done

echo "Wrote: ${SUMMARY}"


# SCRIPT - PART 2 CHECKS ====================================================================================
#!/bin/bash
set -euo pipefail

# Calculating the Overlap of combined with modern and historical lists 
combined_total=$(wc -l < entire_genome_combined_nosex.list)
modern_total=$(wc -l < entire_genome_modern_nosex.list)
historical_total=$(wc -l < entire_genome_historical_nosex.list)

combined_in_modern=$(comm -12 entire_genome_combined_nosex.list entire_genome_modern_nosex.list | wc -l)
combined_in_historical=$(comm -12 entire_genome_combined_nosex.list entire_genome_historical_nosex.list | wc -l)
modern_historical_intersect=$(comm -12 entire_genome_modern_nosex.list entire_genome_historical_nosex.list | wc -l)

combined_in_both=$(comm -12 <(comm -12 entire_genome_combined_nosex.list entire_genome_modern_nosex.list) \
                           <(comm -12 entire_genome_combined_nosex.list entire_genome_historical_nosex.list) | wc -l)

combined_modern_only=$(comm -23 entire_genome_combined_nosex.list entire_genome_historical_nosex.list | wc -l)
combined_historical_only=$(comm -23 entire_genome_combined_nosex.list entire_genome_modern_nosex.list | wc -l)

echo "combined_total	$combined_total"
echo "modern_total	$modern_total"
echo "historical_total	$historical_total"
echo "combined_in_modern	$combined_in_modern"
echo "combined_in_historical	$combined_in_historical"
echo "combined_in_both	$combined_in_both"
echo "combined_modern_only	$combined_modern_only"
echo "combined_historical_only	$combined_historical_only"
echo "modern_historical_intersect	$modern_historical_intersect"

# Calculating the per-indivudal representation in the combined list
sed 's#.*/##; s/.rmdup.bam$//' ${LIST_COMBINED} > beagle_sample_names.txt # Sample names

(
echo -e "sample\tindex\tinformative_sites\tmissing_sites\ttotal_sites\tinformative_fraction"
paste beagle_sample_names.txt <(
    zcat entire_genome_combined_nosex.beagle.gz | awk '
    BEGIN { OFS="\t" }
    NR==1 { next }
    {
        nind = (NF-3)/3
        for (i=1; i<=nind; i++) {
            c1 = 3 + (i-1)*3 + 1
            c2 = c1 + 1
            c3 = c1 + 2
            if ($c1 == 0.333333 && $c2 == 0.333333 && $c3 == 0.333333) {
                missing[i]++
            } else {
                informative[i]++
            }
            total[i]++
        }
    }
    END {
        for (i=1; i<=nind; i++) {
            frac = (total[i] > 0 ? informative[i]/total[i] : 0)
            print i, informative[i]+0, missing[i]+0, total[i]+0, frac
        }
    }'
)
) > beagle_per_individual_counts.tsv


# Calculate usable sites per indivudal ("How many sites from the combined list are actually usable per individual?”)
echo -e "sample\ttotal_sites" > per_individual_callable_sites.tsv

for f in *FINAL*.saf.idx; do
    name=$(basename "$f" .saf.idx)

    total=$(realSFS print "$f" 2>/dev/null | awk 'NF>=2 && $1 ~ /scaffold/ {n++} END{print n}')

    echo -e "${name}\t${total}" >> per_individual_callable_sites.tsv
done

exit 0
