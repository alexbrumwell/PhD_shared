bcftools query -f '[%CHROM\t%SAMPLE\t%GT\t%DP\n]' merged.vcf.gz | \
awk '
BEGIN { OFS="\t" }
{
    scaffold=$1
    sample=$2
    gt=$3
    dp=$4

    total[scaffold]++

    if (gt == "./." || gt == ".|." || gt == ".") {
        missing[scaffold]++
    } else {
        called[scaffold]++

        if (dp != "." && dp != "") {
            sum_dp[scaffold] += dp
            n_dp[scaffold]++
        }
    }
}
END {
    print "Scaffold", "Mean_DP", "N_genotypes_with_DP", "Called_genotypes", "Total_genotypes", "Missing_fraction"

    for (s in total) {
        mean = (n_dp[s] > 0 ? sum_dp[s] / n_dp[s] : "NA")
        miss = missing[s] / total[s]

        print s, mean, n_dp[s], called[s], total[s], miss
    }
}' > scaffold_meanDP_missingness.tsv

# Sort scaffolds
sort -V scaffold_meanDP_missingness.tsv > tmp.tsv && mv tmp.tsv scaffold_meanDP_missingness.tsv
# Inspect suspiciouis scaffolds (>15x and >0.8 missingness)
awk 'NR==1 || $2 > 10 && $6 > 0.8' scaffold_meanDP_missingness.tsv
