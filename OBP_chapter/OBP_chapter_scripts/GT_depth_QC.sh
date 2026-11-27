#!/bin/bash
#SBATCH --job-name=GT_depth_QC
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=jobs/GT_depth_QC_%j.out
#SBATCH --error=jobs/GT_depth_QC_%j.err

set -euo pipefail

module load bcftools/1.22

cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/2_FiltGVCFs/scaffold_30

OUT="genotype_depth_counts.tsv"

echo -e "VCF\tSample\tDepth\tGenotype_class\tCount" > "$OUT"

for VCF in *gz ; do
    [[ -z "$VCF" ]] && continue

    echo "Processing $VCF"

    bcftools query -f '[%SAMPLE\t%GT\t%DP\n]' "$VCF" | \
    awk -v vcf="$(basename "$VCF")" '
    BEGIN { OFS="\t" }
    {
        sample=$1
        gt=$2
        dp=$3

        if (gt=="." || gt=="./." || gt==".|." || dp=="." || dp=="") next

        if (gt=="0/0" || gt=="0|0") class="HomRef"
        else if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") class="Het"
        else if (gt=="1/1" || gt=="1|1") class="HomAlt"
        else next

        count[vcf, sample, dp, class]++
    }
    END {
        for (key in count) {
            split(key, a, SUBSEP)
            print a[1], a[2], a[3], a[4], count[key]
        }
    }' >> "$OUT"

done 

exit 0
# Then after in R:
library(tidyverse)

df <- read_tsv("genotype_depth_counts.tsv",
               col_types = cols(
                 VCF = col_character(),
                 Sample = col_character(),
                 Depth = col_integer(),
                 Genotype_class = col_character(),
                 Count = col_double()
               ))

df <- df %>%
  filter(!is.na(Depth)) %>%
  mutate(
    Genotype_class = factor(Genotype_class, levels = c("HomRef", "Het", "HomAlt"))
  )
#Genotype counts
ggplot(df, aes(x = Depth, y = Count, colour = Genotype_class)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Sample, scales = "free_y") +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Depth",
    y = "Number of genotypes",
    colour = "Genotype class",
    title = "Genotype counts by depth"
  ) 
# Proportions
df_prop <- df %>%
  group_by(Sample, Depth) %>%
  mutate(
    Total_at_depth = sum(Count),
    Proportion = Count / Total_at_depth
  ) %>%
  ungroup()

ggplot(df_prop, aes(x = Depth, y = Proportion, colour = Genotype_class)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Sample) +
  theme_bw() +
  labs(
    x = "Depth",
    y = "Proportion of genotypes at depth",
    colour = "Genotype class",
    title = "Genotype class proportions by depth"
  )

# Zoom in
ggplot(df_prop %>% filter(Depth <= 50),
       aes(x = Depth, y = Proportion, colour = Genotype_class)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Sample) +
  theme_bw() +
  labs(
    x = "Depth",
    y = "Proportion of genotypes at depth",
    colour = "Genotype class",
    title = "Genotype proportions by depth, depth ≤ 50"
  )
# Only hets
ggplot(df %>% filter(Genotype_class == "Het", Depth <= 50),
       aes(x = Depth, y = Count)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Sample, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Depth",
    y = "Number of heterozygous genotypes",
    title = "Heterozygous genotype discovery by depth"
  )
