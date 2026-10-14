# Using ANGSD to generate GL for scaffold 1

## 1. OPTIONAL: Create sites lists

I am actually skipping this step because I want to make the variant caller comparisons more comparable by avoiding a sites list.

```bash
#!/bin/bash
#SBATCH --job-name=ANGSD_site_lists
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/data/ANGSD_GL"
OUTDIR="${DIR}/site_lists"
mkdir -p ${OUTDIR}

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected.list" 
METADATA="/projects/erode/people/qhk515/01_OBP/data/OBP_metadata.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"

cd ${OUTDIR}

# Make BAM list of selected samples
awk 'NR==FNR {samples[$1]; next} $1 in samples {print $2}' "$LIST_SELECTED_SAMPLES" "$METADATA" > selected_bams.list
N_SELECTED=$(wc -l < selected_bams.list)

# List 1: High quality sites Making group lists (requieres at least data from 2 individuals to be included)
LIST1="highqual"

echo "$(timestamp) [*] Starting ANGSD for $LIST1 "
angsd -b selected_bams.list -r scaffold_1 -ref ${REF} -anc ${REF} -out scaf1_${LIST1} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 30 -doCounts 1 -minInd 2 \
    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

# List 2: NO TRANSITIONS - High quality sites Making group lists (requieres at least data from 2 individuals to be included)
LIST2="highqual_notrans"

echo "$(timestamp) [*] Starting ANGSD for $LIST2 "
angsd -b selected_bams.list -r scaffold_1 -ref ${REF} -anc ${REF} -out scaf1_${LIST2} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 30 -doCounts 1 -minInd 2 \
    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

# List 3: NO MISSINGESS - High quality sites Making group lists (requieres at least data from 2 individuals to be included)
LIST3="highqual_notrans_nomissingess"

echo "$(timestamp) [*] Starting ANGSD for ${LIST3} "
angsd -b selected_bams.list -r scaffold_1 -ref ${REF} -anc ${REF} -out scaf1_${LIST3} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
    -C 50 -baq 0 -minMapQ 30 -minQ 30 -doCounts 1 -minInd ${N_SELECTED} \
    -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

echo "$(timestamp) [OK] Site lists created"

# Extract sites lists
for LIST in "$LIST1" "$LIST2" "$LIST3"; do
    zcat scaf1_${LIST}.mafs.gz | awk 'NR > 1 {print $1, $2}' > scaf1_${LIST}.list
    angsd sites index scaf1_${LIST}.list
done

# Convert site lists to BED files so that I can filter VCFs later
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' scaf1_highqual.list > scaf1_highqual.bed
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' scaf1_highqual_notrans.list > scaf1_highqual_notrans.bed
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' scaf1_highqual_notrans_nomissingess.list > scaf1_highqual_notrans_nomissingess.bed

echo "$(timestamp) [OK] Site lists filtered to just positions created"


```

## 2. Run ANGSD per sample

```bash
# SCRIPT PART 2 ====================================================================================
#!/bin/bash
#SBATCH --job-name=ANGSD_GL_scaf1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%5
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/data/ANGSD_GL"
OUTDIR="${DIR}/scaffold_1"
METADATA="/projects/erode/people/qhk515/01_OBP/data/OBP_metadata.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
mkdir -p $DIR $OUTDIR

cd $OUTDIR

TORUN=$1

# Get sample
read -r chrom sample < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TORUN")

BAM=$(awk -v sample="$sample" '$1 == sample {print $2}' "$METADATA")
NAME=$(basename "${BAM}" .rmdup.bam)
# Check BAM exists
if [[ ! -f "$BAM" ]]; then
    echo "[${timestamp}] ERROR: BAM file not found: $BAM"
    exit 1
fi

# Run ANGSD: Requierements (at least >4x coverage and <50x coverage per site to be included in the analysis)
echo "$(timestamp) [*] Starting ANGSD for ${NAME}"
echo "$(timestamp) [*] Region: ${chrom}"
echo "$(timestamp) [*] BAM: ${BAM}"

angsd -i ${BAM} -ref ${REF} -r scaffold_1 -anc ${REF} \
    -out ${NAME}_scaffold1_tv_minDP4_maxDP50 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
    -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
    -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1 \
    -doGlf 2 -doMajorMinor 1 -doMaf 2 
		#-dopost 1 -dobcf 1 --ignore-RG 0 -dogeno 1


echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

# Now convert calculate SFS and convert to wimSFS format (winsfs can only fold afterwards but this shouldn't affect heterozygosity calculations)
winsfs "${NAME}_scaffold1_tv_minDP4_maxDP50.saf.idx" > "${NAME}_scaffold1_tv_minDP4_maxDP50.winsfs.ml"

echo "$(timestamp) [OK] winsfs completed for ${NAME}"
echo "$(timestamp) [OK] Output: ${NAME}_scaffold1_tv_minDP4_maxDP50.winsfs.ml"

```

After, to run only the 28 OBP samples I am interested in:

```bash
cd /projects/erode/people/qhk515/01_OBP/data/ANGSD_GL


for chunk in /projects/erode/people/qhk515/01_OBP/data/ToRun_chunk_00.txt
do
    n=$(wc -l < "$chunk")

    chunk_base=$(basename "$chunk")

    echo "$chunk_base"

    sbatch --array=1-28 ANGSD_GL_scaf1.sh "$chunk"

done
```

Afterwards:

```bash

# To get summary het values across all individuals, run the following code after all jobs are done:
SUMMARY="heterozygosity_summary.tsv"

echo -e "sample\tHomRef\tHet\tHomAlt\tTotal_sites\tHeterozygosity" > "${SUMMARY}"

shopt -s nullglob

for file in *_scaffold1_tv_minDP4_maxDP50.winsfs.ml; do
  sample=$(basename "${file}" _scaffold1_tv_minDP4_maxDP50.winsfs.ml)

  grep -v "SHAPE" "${file}" | awk -v sample="${sample}" '
    {
      hom_ref = $1
      het     = $2
      hom_alt = $3

      total = hom_ref + het + hom_alt

      if (total > 0) {
        H = het / total
        printf("%s\t%.6f\t%.6f\t%.6f\t%.0f\t%.10f\n",
               sample, hom_ref, het, hom_alt, total, H)
      } else {
        printf("%s\tNA\tNA\tNA\t0\tNA\n", sample)
      }
    }
  ' >> "${SUMMARY}"
done

echo "Wrote: ${SUMMARY}"
```

# 2.5: Check sites lists

```bash
# SCRIPT - PART 3 CHECKS ====================================================================================
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

for f in *_nosex_filtlist.saf.idx; do
    name=$(basename "$f" .saf.idx)

    total=$(realSFS print "$f" 2>/dev/null | awk 'NF>=2 && $1 ~ /scaffold/ {n++} END{print n}')

    echo -e "${name}\t${total}" >> per_individual_callable_sites.tsv
done

exit 0
```