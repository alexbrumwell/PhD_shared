# Calculate depth percentiles per scaffold
echo -e "Sample\tN_sites\tMean_DP\tDP_1pct\tDP_99pct" > all_depth_percentiles.tsv

for VCF in *.vcf.gz; do

    SAMPLE=$(basename "$VCF" .vcf.gz)

    bcftools query -f '%CHROM\t%POS\t[%DP]\n' "$VCF" | \
    awk -v sample="$SAMPLE" '
    BEGIN { OFS = "\t" }
    $3 != "." && $3 != "" {
        dp = $3
        depth[dp]++
        n++
        sum += dp
    }
    END {
        p1_target = n * 0.01
        p99_target = n * 0.99

        cumulative = 0
        p1 = "NA"
        p99 = "NA"

        for (d = 0; d <= 1000000; d++) {
            if (d in depth) {
                cumulative += depth[d]

                if (p1 == "NA" && cumulative >= p1_target) p1 = d

                if (p99 == "NA" && cumulative >= p99_target) {
                    p99 = d
                    break
                }
            }
        }

        mean = sum / n

        print sample, n, mean, p1, p99
    }' >> all_depth_percentiles.tsv

done

# Depth percentiles from BAM files 
#!/bin/bash
#SBATCH --job-name=bam_depth_pct
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=bam_depth_pct_%j.log
#SBATCH --error=bam_depth_pct_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk
#SBATCH --qos=normal

module load samtools

BAM_DIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP"

cd $BAM_DIR

echo -e "Sample\tN_covered_sites\tMean_depth_covered_sites\tDP_1pct\tDP_99pct" > "$OUT"

for BAM in "$BAM_DIR"/*.bam; do

    [ -e "$BAM" ] || continue

    SAMPLE=$(basename "$BAM" .bam)

    echo "Processing $SAMPLE"

    samtools depth -@ "$SLURM_CPUS_PER_TASK" "$BAM" | \
    awk -v sample="$SAMPLE" '
    BEGIN { OFS = "\t" }
    {
        dp = $3
        depth[dp]++
        n++
        sum += dp
    }
    END {
        if (n == 0) {
            print sample, 0, "NA", "NA", "NA"
            exit
        }

        p1_target = n * 0.01
        p99_target = n * 0.99

        cumulative = 0
        p1 = "NA"
        p99 = "NA"

        for (d = 1; d <= 1000000; d++) {
            if (d in depth) {
                cumulative += depth[d]

                if (p1 == "NA" && cumulative >= p1_target) p1 = d

                if (p99 == "NA" && cumulative >= p99_target) {
                    p99 = d
                    break
                }
            }
        }

        mean = sum / n

        print sample, n, mean, p1, p99
    }' >> "$OUT"

done
