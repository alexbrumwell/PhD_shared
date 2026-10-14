# Downsampling variant calling pipeline

This is the markdown for the downsampling pipeline.
I decided to downsample everything to 6x because that is the minimum depth threshold I am using for OBP samples in general, so it makes sense to use that sample threshold to downsample everything to.

## 1. OPTIONAL: Filtering BAM file

In order to avoid removing a higher fraction of reads than necessary, I tested filtering for mapping quality in the BAM files prior to downsampling (if not, I might calcaulte a higher mean coverage which would result in a bigger fraction of reads and then afterwards filtering for quality and getting lower average depth than the target downsampling coverage).

To test, I will look at:

- mean depth after downsampling
- callable sites
- missingness
- heterozygosity
- Ti/Tv
- 

```bash
#!/bin/bash
#SBATCH --job-name=Filter_BAM 
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

timestamp() { date '+%F %T'; }

# === Load modules ===
module load samtools/1.21

DIR_BAM="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP" # Direcotry where BAM files are located
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
GENOME_SIZE=$(awk '{sum += $2} END {print sum}' "${REF}.fai")

OUTDIR="${DIR_BAM}/filt_bam_q20"
mkdir -p "$OUTDIR"


SUMMARY="${OUTDIR}/summary_stats.tsv"
: > $SUMMARY

echo -e "Sample\tStage\tTotal_Reads\tMapped_Reads\tMapping_Percentage" > "$SUMMARY"

# Script to filter BAM files for quality
for BAM in ${DIR_BAM}/*rmdup.bam; do
    sample=$(basename $BAM .rmdup.bam)
    OUTBAM="${OUTDIR}/${sample}.q20.bam"

    echo "[$(timestamp)] Processing $sample"


    samtools view -@ "$SLURM_CPUS_PER_TASK" -b -q 10 $BAM -o $OUTBAM
    samtools index -@ "$SLURM_CPUS_PER_TASK" $OUTBAM
    samtools flagstat -@ "$SLURM_CPUS_PER_TASK" $OUTBAM > "${OUTBAM%.bam}.flagstat.txt"
		samtools stats -@ "$SLURM_CPUS_PER_TASK" $OUTBAM > "${OUTBAM%.bam}.stats.txt"
		
	echo "[$(timestamp)] Calculating summary stats for $sample"

	BEFORE_FLAGSTAT="${BAM}.flagstats.txt"
	AFTER_FLAGSTAT="${OUTBAM%.bam}.flagstat.txt"

    for STAGE in before after; do

    if [[ "$STAGE" == "before" ]]; then
        THISBAM="$BAM"
        FLAGSTAT="$BEFORE_FLAGSTAT"
    else
        THISBAM="$OUTBAM"
        FLAGSTAT="$AFTER_FLAGSTAT"
    fi

    total_reads=$(awk 'NR==1 {print $1}' "$FLAGSTAT")
    mapped_reads=$(awk '/ mapped \(/ {print $1; exit}' "$FLAGSTAT")

    if [[ "$total_reads" -gt 0 ]]; then
        mapping_percentage=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN {printf "%.2f", (m/t)*100}')
    else
        mapping_percentage="NA"
    fi

    mean_depth=$(samtools depth -a -@ "$SLURM_CPUS_PER_TASK" "$THISBAM" | \
        awk -v genome="$GENOME_SIZE" '{sum += $3} END {printf "%.4f", sum/genome}')

    echo -e "${sample}\t${STAGE}\t${total_reads}\t${mapped_reads}\t${mapping_percentage}\t${mean_depth}" >> "$SUMMARY"

	done

done 

echo "[$(timestamp)] Done"

```

> #### Results
>
> | Sample          | Cov_Before | Cov_After | Cov_Lost | %Cov_Lost |
> | --------------- | ---------: | --------: | -------: | --------: |
> | O-76339-001     | 7.5522     | 6.0814    | 1.4708   | 19.47     |
> | O-76349-001     | 7.9413     | 6.4072    | 1.5341   | 19.32     |
> | O-76350-001     | 7.5660     | 6.1380    | 1.4280   | 18.87     |
> | O-76521-001     | 7.5265     | 6.0704    | 1.4561   | 19.35     |
> | O-76523-001     | 7.5908     | 6.1259    | 1.4649   | 19.30     |
> | O-76526-001     | 7.5756     | 6.0911    | 1.4845   | 19.60     |
> | O-76538-001     | 7.5525     | 6.0895    | 1.4630   | 19.37     |
> | O-81171-001     | 7.5833     | 6.6602    | 0.9231   | 12.17     |
> | O-81175-001     | 7.5298     | 6.5339    | 0.9959   | 13.23     |
> | O-81179-001     | 7.7679     | 6.6824    | 1.0855   | 13.97     |
> | O-81182-001     | 7.5424     | 6.6794    | 0.8630   | 11.44     |
> | O-81186-001     | 7.3618     | 6.4631    | 0.8987   | 12.21     |
> | O-81188-001     | 7.5358     | 6.5019    | 1.0339   | 13.72     |
> | O-81192-001     | 7.5684     | 6.5986    | 0.9698   | 12.81     |
> | O-81193-001     | 7.5684     | 6.7100    | 0.8584   | 11.34     |
> | O-81196-001     | 7.5167     | 6.4971    | 1.0196   | 13.56     |
> | O-81198-001     | 6.4975     | 5.6852    | 0.8123   | 12.50     |
> | O-81201-001     | 7.5389     | 6.6368    | 0.9021   | 11.97     |
> | O-81203-001     | 7.5620     | 6.6266    | 0.9354   | 12.37     |
> | OBPI_501        | 24.2892    | 19.7443   | 4.5449   | 18.71     |
> | OBPI_502        | 9.9758     | 7.7021    | 2.2737   | 22.79     |
> | OBPI_503        | 7.2480     | 5.2779    | 1.9701   | 27.18     |
> | OBPI_504        | 17.6808    | 14.4276   | 3.2532   | 18.40     |
> | OBPI_505        | 2.4213     | 1.8718    | 0.5495   | 22.69     |
> | OBPI_506        | 4.8643     | 3.8663    | 0.9980   | 20.52     |
> | OBPI_507        | 5.7727     | 4.4795    | 1.2932   | 22.40     |
> | OBPI_508        | 3.9134     | 3.0113    | 0.9021   | 23.05     |
> | OBPII_ReHoI_501 | 5.0992     | 4.0026    | 1.0966   | 21.51     |
> | OBPII_ReHoI_502 | 6.7085     | 5.4863    | 1.2222   | 18.22     |
> | OBPII_ReHoI_503 | 3.9404     | 3.3008    | 0.6396   | 16.23     |
> | OBPII_ReHoI_504 | 7.6062     | 6.1317    | 1.4745   | 19.38     |
> | OBPII_ReHoI_505 | 2.6354     | 2.1131    | 0.5223   | 19.82     |
> | OBPII_ReHoI_506 | 4.0322     | 3.1628    | 0.8694   | 21.56     |
> | OBPII_ReHoI_507 | 9.0666     | 7.2934    | 1.7732   | 19.56     |
> | OBPII_ReHoI_508 | 10.2505    | 8.2636    | 1.9869   | 19.38     |

Because there seemed to be a big drop in mean coverage. I decided to not use these files as the mapping quality filtering might be too strict. Although I am thinking of using these files later for some tests. 

> Note: variant calling with these files

## 0. Create downsampled BAM files for target coverage
I downsample to 6x:


```bash
#!/bin/bash
#SBATCH --job-name=Downsample 
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

timestamp() { date '+%F %T'; }

# === Load modules ===
module load samtools/1.21

# === MODIFY ===
TARGET_COVERAGES=("6")  # target coverages to downsample to

DIR_BAM="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP" # Direcotry where BAM files are located
LIST_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_bam_samples_selected_6x.list" # List of sample names (without .rmdup.bam extension) to process

OUTDIR="${DIR_BAM}/Downsampled"  # Outdirectory where downsampled BAM files will be saved
OUTDIR_6X="${OUTDIR}/6x"
#OUTDIR_5X="${OUTDIR}/5x"
mkdir -p "$OUTDIR" "$OUTDIR_6X"


# ===== Output summary files =====
QC_OUT="${OUTDIR}/downsampling_qc.tsv" # Stores downsampling information and observed coverage after downsampling

LIST_DOWN6X="${OUTDIR}/Down_6x.list" # Lists BAMs to use for the 6x and 5x datasets
#LIST_DOWN5X="${OUTDIR}/Down_5x.list"

# Empty old output files before starting
: > "$QC_OUT"
: > "$LIST_DOWN6X"


# Add header to the QC table
echo -e "Sample\tInput_BAM\tOriginal_coverage\tTarget_coverage\tFraction\tOutput_BAM\tObserved_coverage" > "$QC_OUT"


# ===== Function to calculate mean genome-wide depth =====
calc_cov() {
    local bam="$1"

    # samtools depth -a reports depth at every position, including zero-depth sites.
    # awk sums depth across all positions and divides by the number of positions.
    samtools depth -@ "$SLURM_CPUS_PER_TASK" -a "$bam" | awk '{sum += $3; n++} END {if (n > 0) print sum/n; else print "NA"}'
}

# ===== Function to downsample one BAM to one target coverage =====
downsample_if_needed() {
    local bam="$1"        # Input BAM
    local sample="$2"     # Sample name
    local coverage="$3"   # Original observed coverage
    local target="$4"     # Target coverage, e.g. 7 or 5
    local outdir="$5"     # Output directory
    local list_out="$6"   # File to store output BAM path

    local output_bam="${outdir}/${sample}.down${target}x.bam"

    # Only downsample if original coverage is greater than the target coverage
    if awk -v c="$coverage" -v t="$target" 'BEGIN {exit !(c > t)}'; then

        # Calculate the fraction of reads to keep
        # Example: if coverage is 10x and target is 5x, fraction = 5/10 = 0.5
        frac=$(awk -v t="$target" -v c="$coverage" 'BEGIN {print t/c}')

        echo "[$(timestamp)] Downsampling ${sample}: ${coverage}x -> ${target}x; fraction=${frac}"

        # Only run downsampling if output BAM does not already exist
        if [[ ! -f "$output_bam" ]]; then

            # Downsample BAM using samtools view.
            # -b writes BAM output.
            # -s keeps a random fraction of reads.
            samtools view \
                -@ "$SLURM_CPUS_PER_TASK" \
                -b \
                -s "${frac}" \
                "$bam" \
                -o "$output_bam"

            # Index the downsampled BAM
            samtools index -@ "$SLURM_CPUS_PER_TASK" "$output_bam"
						samtools stats -@ "$SLURM_CPUS_PER_TASK" $output_bam > ${output_bam%.bam}.stats.txt
        else
            echo "[$(timestamp)] Output already exists: $output_bam"
        fi

        # Calculate observed coverage after downsampling
        observed=$(calc_cov "$output_bam")

        # Add downsampled BAM path to the relevant BAM list
        echo "$output_bam" >> "$list_out"

    else
        # If original coverage is already below the target, keep the original BAM
        echo "[$(timestamp)] ${sample} has coverage ${coverage}x, which is <= ${target}x. Keeping original BAM."

        output_bam="$bam"
        frac="NA"
        observed="$coverage"

        # Add original BAM path to the relevant BAM list
        echo "$output_bam" >> "$list_out"
    fi

    # Write one line to the QC summary table
    echo -e "${sample}\t${bam}\t${coverage}\t${target}\t${frac}\t${output_bam}\t${observed}" >> "$QC_OUT"
}


# ===== Main loop over samples =====
while read -r bam; do

    # Check that the BAM exists before processing
    if [[ ! -f "$bam" ]]; then
        echo "Error: BAM file not found: $bam" >&2
        continue
    fi

    # Construct the expected BAM path
    sample=$(basename $bam .rmdup.bam)


    # Calculate original coverage
    echo "[$(timestamp)] Calculating coverage for $sample..."
    coverage=$(calc_cov "$bam")

    # Downsample to 6x and 5x
    downsample_if_needed "$bam" "$sample" "$coverage" "6" "$OUTDIR_6X" "$LIST_DOWN6X"
    #downsample_if_needed "$bam" "$sample" "$coverage" "5" "$OUTDIR_5X" "$LIST_DOWN5X"

done < "$LIST_SAMPLES"

# ===== Final log messages =====
echo "[$(timestamp)] Downsampling complete."
echo "Original coverage table: $LIST_COVERAGE"
echo "Downsampling QC table: $QC_OUT"
echo "6x BAM list: $LIST_DOWN6X"
```

And this will create our two sets of downsampled BAM files, as well as QC text files to ensure that we have downsampled correctly.


## 1. ANGSD 

Due to the low coverage of some historical samples, we used ANGSD to generate genotype likelihood files.

First create a list of BAM files at or over the target downsampled coverage (i.e. samples that have been downsampled to 7x or at 7x).  I am using N=26 total samples. 5 samples hisotircal samples have been downsampled because they had a meandcoverage of >9x. The remaining 21 samples havec a mean coverage of about 7x (although 2 samples are at 6.5X).

#### 1.1 Standard ANGSD method - INDIVIDUAL

Then, I wrote what I thought would be the best approach to reduce group biases and coverage biases when running ANGSD.

It differs to Lukes script by:

Only uses scaffolds >100kb

individual max depth rather than a general cap at 50X to avoid low complexity sites (i.e. repetitive)

I based this script on Fontsere et al. 2025 (Whooping Crane paper).7.1.B Standard ANGSD method - INDIVIDUAL

```bash
#!/bin/bash
#SBATCH --job-name=Het_ANGSD_ind
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-26%10
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/Down_6x.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
#FILTERED_META="/projects/erode/people/qhk515/01_OBP/data/OBP_2x_meandepth.list"

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Down6x/het_ANGSD_ind"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST_COMBINED}")
NAME=$(basename "${BAM}" .bam)
#maxdepth=$(awk -v s="$NAME" '$1 == s {print $3}' "$FILTERED_META")

# Step 1: Calculate high quality sites per individual 
echo "$(timestamp) [*] Creating sites list for individual ${NAME}"

angsd -i ${BAM} -ref ${REF} -rf ${AUTOSOMAL_SCAFFOLDS} -anc ${REF} \
-out ${NAME}_min4x_down6x \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans 1 -trim 0 -C 50 \
-baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 15 -doCounts 1 -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 \
-nThreads ${SLURM_CPUS_PER_TASK}

zcat ${NAME}_min4x_down6x.mafs.gz| awk '{print $1,$2}' | sed 1d > ${NAME}_min4x_down6x.list && angsd sites index ${NAME}_min4x_down6x.list

# Step 2: Run ANGSD
echo "$(timestamp) [*] Creating GL for individual ${NAME}"

angsd  -i ${BAM} -rf ${AUTOSOMAL_SCAFFOLDS} -ref $REF -sites ${NAME}_min4x_down6x.list -anc $REF -out ${NAME}_min4x_down6x_FINAL -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1  -noTrans 1 -C 50 -baq 0 -minMapQ 30 -minQ 20  -setMinDepth 4 -setMaxDepth 15 -doCounts 1 -nThreads $SLURM_CPUS_PER_TASK -GL 2 -doSaf 1

# Now convert calculate SFS and convert to wimSFS format 
winsfs ${NAME}_min4x_down6x_FINAL.saf.idx > ${NAME}_min4x_down6x_FINAL.ml

echo "$(timestamp) [OK] winSFS completed for ${NAME}"


exit 0

# ------------------------
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
```


#### 7.3 Standard ANGSD method -  PER GROUP COMMON SITES 

I am using the same parameters

```bash
#!/bin/bash
#SBATCH --job-name=Het_ANGSD_sharedlist80_pergroup
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --output=jobs/%x_%j.out
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# Combines: (1) Luke's population-level site-discovery step, but now requiring
# data present in >=80% of individuals WITHIN EACH GROUP separately (Modern,
# Historical) and taking the intersection of the two passing site lists --
# mirroring the per-group missingness + bedtools intersect logic used in the
# snpAD/bcftools VCF pipeline -- rather than a single global -minInd cutoff
# applied across the combined cohort; and (2) script 2's per-individual SAF
# step with individual-tailored max-depth cutoffs, looped serially over all
# individuals against the shared site list.

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bedtools/2.30.0 

# === Inputs ====
LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/Down_6x.list"
AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"

MODERN_PREFIX="O-"
HIST_PREFIX="OBP"

MISS_FRAC=0.8   # require data in >=80% of individuals, per group

OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Down6x/ANGSD_sharedlist80_pergroup"
mkdir -p ${OUTDIR}

cd ${OUTDIR}

# =========================
# 0. Split combined BAM list into per-group BAM lists
# =========================
# Matches sample names against the same prefixes used in the snpAD pipeline
# (MODERN_PREFIX="O-", HIST_PREFIX="OBP"). Assumes BAM basenames (minus
# .rmdup.bam) start with one of these prefixes, as in script 2's NAME= line.

MODERN_LIST="${OUTDIR}/Modern_bam.list"
HIST_LIST="${OUTDIR}/Historical_bam.list"

# Skip BAM list splitting if already done
echo "$(timestamp) [*] Creating per-group BAM lists"
: > "${MODERN_LIST}"
:  > "${HIST_LIST}"

    while read -r BAM; do
        NAME=$(basename "${BAM}" .bam)
        if [[ "${NAME}" == ${HIST_PREFIX}* ]]; then
            echo "${BAM}" >> "${HIST_LIST}"
        elif [[ "${NAME}" == ${MODERN_PREFIX}* ]]; then
            echo "${BAM}" >> "${MODERN_LIST}"
        else
            echo "WARNING: ${NAME} did not match ${HIST_PREFIX}* or ${MODERN_PREFIX}* -- not assigned to a group" >&2
        fi
    done < "${LIST_COMBINED}"

    n_all=$(wc -l < "${LIST_COMBINED}")
    n_grouped=$(( $(wc -l < "${HIST_LIST}") + $(wc -l < "${MODERN_LIST}") ))
    if [[ "${n_grouped}" -ne "${n_all}" ]]; then
        echo "WARNING: ${n_all} total BAMs but only ${n_grouped} matched a group prefix -- check sample naming" >&2
    fi

[[ -s "${MODERN_LIST}" ]] || { echo "ERROR: ${MODERN_LIST} empty -- no Modern samples found" >&2; exit 1; }
[[ -s "${HIST_LIST}" ]]   || { echo "ERROR: ${HIST_LIST} empty -- no Historical samples found" >&2; exit 1; }

# =========================
# 1. Per-group site discovery, each requiring data in >=80% of THAT group's individuals
# =========================
declare -A GROUP_LIST=( [Modern]="${MODERN_LIST}" [Historical]="${HIST_LIST}" )
declare -A GROUP_SITELIST

for GRP in Modern Historical; do
    LIST="${GROUP_LIST[${GRP}]}"
    N_IND=$(wc -l < "${LIST}")
    MININD=$(awk -v n="${N_IND}" -v p="${MISS_FRAC}" 'BEGIN{v=n*p; iv=int(v); print (v>iv)?iv+1:iv}')
    echo "$(timestamp) [*] ${GRP}: ${N_IND} individuals -> requiring data in >=${MININD} individuals (${MISS_FRAC})"

    # Skip ANGSD if site list already exists
    if [[ -s "${GRP}_80pct.sites.bed" ]]; then
        echo "$(timestamp) [*] ${GRP}_80pct.sites.bed already exists, skipping ANGSD"
    else
        echo "$(timestamp) [*] Starting ANGSD site discovery for ${GRP}"
        angsd -b "${LIST}" -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -out "${GRP}_80pct" \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -rmTrans 1 \
            -C 50 -baq 0 -minMapQ 30 -minQ 20 -doCounts 1 -minInd ${MININD} -setMinDepthInd 4 -setMaxDepthInd 14 \
            -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -nThreads ${SLURM_CPUS_PER_TASK}

        echo "$(timestamp) [OK] ANGSD completed for ${GRP}, extracting site list"

        zcat "${GRP}_80pct.mafs.gz" | awk 'NR>1 {print $1"\t"$2-1"\t"$2}' > "${GRP}_80pct.sites.bed"
    fi

    n_sites=$(wc -l < "${GRP}_80pct.sites.bed")
    echo "$(timestamp) [OK] ${GRP}: ${n_sites} sites pass >=${MISS_FRAC} presence"

    GROUP_SITELIST[${GRP}]="${OUTDIR}/${GRP}_80pct.sites.bed"
done

# =========================
# 2. Intersect: sites passing in BOTH groups (this IS the final shared site list)
# =========================
echo "$(timestamp) [*] Intersecting Modern and Historical passing site lists"

# Create a genome file from your reference or scaffolds list (skip if exists)
if [[ ! -f genome.txt ]]; then
    if [[ ! -f ${REF}.fai ]]; then
        samtools faidx ${REF}
    fi
    cut -f1,2 ${REF}.fai > genome.txt
fi

bedtools intersect -sorted -g genome.txt \
    -a "${GROUP_SITELIST[Modern]}" \
    -b "${GROUP_SITELIST[Historical]}" \
    > entire_genome_combined_80pct.bed

awk '{print $1"\t"$3}' entire_genome_combined_80pct.bed > entire_genome_combined_80pct.list

n_shared=$(wc -l < entire_genome_combined_80pct.list)
echo "$(timestamp) [OK] ${n_shared} shared sites pass in both groups"

if [[ "${n_shared}" -eq 0 ]]; then
    echo "WARNING: no shared passing sites across groups" >&2
    exit 1
fi

angsd sites index entire_genome_combined_80pct.list

echo "$(timestamp) [OK] Shared per-group 80%-presence site list ready: ${OUTDIR}/entire_genome_combined_80pct.list"

# === Part 2: per-individual SAF restricted to the shared list, looped serially over ALL individuals (both groups) ===
SITES_80PCT="${OUTDIR}/entire_genome_combined_80pct.list"

for BAM in $(cat ${LIST_COMBINED}); do
    NAME=$(basename "${BAM}" .down6x.bam)

    # Skip individual if SAF/ML already exists
    if [[ -f "${NAME}_sharedlist80_min4x_down6x.ml" ]]; then
        echo "$(timestamp) [*] ${NAME}_sharedlist80_min4x.ml already exists, skipping"
        continue
    fi

    maxdepth=14

    echo "$(timestamp) [*] Starting ANGSD for individual ${NAME} (maxdepth=${maxdepth})"

    angsd -i ${BAM} -rf ${AUTOSOMAL_SCAFFOLDS} -ref ${REF} -anc ${REF} -sites ${SITES_80PCT} -out ${NAME}_sharedlist80_min4x_down6x \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
        -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth ${maxdepth} -doCounts 1 \
        -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1

    echo "$(timestamp) [OK] ANGSD completed for ${NAME}"

    winsfs ${NAME}_sharedlist80_min4x_down6x.saf.idx > ${NAME}_sharedlist80_min4x_down6x.ml

    echo "$(timestamp) [OK] winSFS completed for ${NAME}"
done


exit 0

# To get summary het values across all individuals, run the following code after the job is done:

SUMMARY="heterozygosity_summary.tsv"

echo -e "sample\tn_boot\tmean_H\tsd_H\tmean_sites\tsd_sites" > "${SUMMARY}"

shopt -s nullglob

for file in *.ml; do
  sample=$(basename "${file}" _sharedlist80_min4x_down6x.ml)

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


```

## 2. Variant calling w/ snpAD

### 2.1 Variant call using snpAD

> To submit the scripts in a chain:
>
```bash
# SCRIPT TO CHAIN SCRIPTS FOR PIPELINE 
DIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x"
chunk="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/ToRun_chunk_00.txt"

# Submit script 0
n=$(wc -l < "$chunk")
cd ${DIR}/1_VCFs
 job0=$(sbatch --parsable --time=16:00:00 --array=1-${n}%100 1_snpad.sh $chunk )
 echo "Submitted script1.sh as job ${job0}"

job0="45719123"

# Submit script 1
n=$(wc -l < "$chunk")
cd ${DIR}/2_FilterGVCFs
 job1=$(sbatch --parsable --dependency=afterok:${job0} --time=03:00:00 --array=1-${n}%100 2_Filter.sh $chunk )
 echo "Submitted script1.sh as job ${job1}"

 # Submit script 2 after script 1 finishes successfully
 cd ${DIR}/3_MergeGVCFs_perChrm
 job2=$(sbatch --parsable --dependency=afterok:${job1} --time=07:00:00 --array=1-357%100 3_Merge.sh)
 echo "Submitted script2.sh as job ${job2}"

 # Submit script 3 after script 2 finishes successfully
 cd ${DIR}/4_GroupMissFiltered
 job3=$(sbatch --parsable --dependency=afterok:${job2} --time=07:00:00 --array=1-357%100 4_Miss_filter.sh)
 echo "Submitted script3.sh as job ${job3}"

 # Submit script 4 after script 3 finishes successfully
 cd ${DIR}/5_bisnps
 job4=$(sbatch --parsable --time=05:00:00 --dependency=afterok:${job3} 5_bisnps_nosingletons.sh)
 echo "Submitted script4.sh as job ${job4}"

 # Submit script 5 after script 4 finishes successfully
 cd ${DIR}/6_concat
 job5=$(sbatch --parsable --dependency=afterok:${job4} --time=01:00:00 6_concat.sh)
 echo "Submitted script5.sh as job ${job5}"

 ```

We can re-use the script used for calling the normal BAM files

First make a text file of the scaffolds and samples of interest. 

```bash
#!/bin/bash
#SBATCH --job-name=0_listcreation
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=11
#SBATCH --time=00:00:05
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

LIST_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/Down_6x.list"
LIST_SCAF="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/

: > ToRun.txt

while IFS= read -r scaffold; do
  while IFS= read -r sample; do
    ID=$(basename $sample .down6x.bam)
    printf '%s\t%s\n' "$scaffold" "$ID" >> ToRun.txt
  done < $LIST_SAMPLES
done < $LIST_SCAF

# Filter for whichever samples and scaffolds to downsample
#grep -w "scaffold_1" ToRun.txt > tmp.txt && mv tmp.txt ToRun.txt

# Now split the text file into max. 9999 lines for SLURM submission
split -l 9999 -d --additional-suffix=.txt ToRun.txt ToRun_chunk_

```

Next run a script to create multiple SLURM scripts (1 per indivudal and scaffold)

```bash
#!/bin/bash
#SBATCH --job-name=1_snpad
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

timestamp() { date '+%F %T'; }

module purge
module load snpad/0.3.11 bcftools/1.22 libdeflate/1.21 htslib perl/5.38.0 gsl/2.5

TORUN=$1

DIR_BAM="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/Downsampled/6x" # Direcotry where BAM files are located

read -r chrom sample < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TORUN")

# Get BAM file
BAM_FILE="${DIR_BAM}/${sample}.down6x.bam" # Extracts the BAM file path for the given sample
if [[ ! -f "$BAM_FILE" ]]; then
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

# Set directories
DIR_OUTPUT1="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x" # Output directory for VCF files
mkdir -p "${DIR_OUTPUT1}" # Create output directory if it doesn't exist

DIR_OUTPUT2=${DIR_OUTPUT1}"/1_VCFs" # Output directory for VCF files
mkdir -p "${DIR_OUTPUT2}" # Create output directory if it doesn't exist

# Define reference
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"  # Reference genome

# Define output files and directory
VCF_FILE=${DIR_OUTPUT2}/${chrom}/${sample}_${chrom}.g.vcf
SNPAD_FILE=${DIR_OUTPUT2}/${chrom}/${sample}_${chrom}.snpAD
PRIORS_FILE=${DIR_OUTPUT2}/${chrom}/${sample}.priors.txt
ERRORS_FILE=${DIR_OUTPUT2}/${chrom}/${sample}.errors.txt
mkdir -p ${DIR_OUTPUT2}/${chrom}/


# Check if the script already ran succesfully for this sample and scaffold by checking if the stats file already exists
if [[ -s "${VCF_FILE}.stats" ]] ; then
    echo "$(timestamp) Raw gVCF has already been created before! No need to run script. Exiting."
    exit 0
else
    # Info for SLURM log file
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${TIMESTAMP}] Processing: ${chrom} for sample: ${sample}"

    # Run Bam2snpAD to preprocess the BAM file
    Bam2snpAD -r ${chrom} -f ${REF} -Q 30 ${BAM_FILE} > ${SNPAD_FILE}

    # Log file message
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    if [ $? -eq 0 ]; then
        echo "[${TIMESTAMP}] Preprocessing done. Next: Calculating priors and errors"
    else
        echo "ERROR: Preprocessing failed" 
        exit 1
    fi

    # Run snpAD to calculate priors and errors
    snpAD -c ${SLURM_CPUS_PER_TASK} -o ${PRIORS_FILE} -O ${ERRORS_FILE} ${SNPAD_FILE}

    # Log file message
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    if [ $? -eq 0 ]; then
        echo "[${TIMESTAMP}] 1st snpAD done. Next: Running snpADCall"
    else
        echo "ERROR: snpAD failed" 
        exit 1
    fi

    # Check if there are zeros in the prior file and run snpADCall
    # When one of the priors is zero snpADCall crashs, and in previous runs this was happening very often in the scaffolds, and I think in very low cov samples)  
    ZERO_PRIOR=".00000000075556111967"  # Set zero prior to 1/genome size (1323519665 for OBP)

    if awk -F, '{for(i=1;i<=NF;i++) if($i==0) exit 1}' "${PRIORS_FILE}";
    then
    # No zero found
        snpADCall -N ${sample} -e ${ERRORS_FILE} -p "$(cat ${PRIORS_FILE})" ${SNPAD_FILE} > ${VCF_FILE}
    else
    # Zero found (use --set_zero_priors with 1/genomesize as recomended)
    snpADCall -N ${sample} --set_zero_priors ${ZERO_PRIOR} -e ${ERRORS_FILE} -p "$(cat ${PRIORS_FILE})" ${SNPAD_FILE} > ${VCF_FILE}
    fi


    # Log file message
    if [ $? -eq 0 ]; then    
        TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
        echo "[${TIMESTAMP}] snpAD calling done. Next: indexing VCF files"
    else
        echo "ERROR: variant calling failed" 
        exit 1
    fi

    # Index the VCF file
    bgzip "${VCF_FILE}" && bcftools index "${VCF_FILE}.gz"
    bcftools stats -s - "${VCF_FILE}.gz" > "${VCF_FILE}.stats"


    # Log file message
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    if [ $? -eq 0 ]; then
        echo "[${TIMESTAMP}] Indexing complete. Script fully run."
    else
        echo "ERROR: indexing failed"
        exit 1
    fi


    # Info for SLURM log file
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${TIMESTAMP}] Finished ${chrom} for ${sample}" 
fi
```

Now you can submit multiple version of the previous script doing: 

```bash
cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/


for chunk in ToRun_chunk_*.txt
do
    n=$(wc -l < "$chunk")

    chunk_base=$(basename "$chunk")

    echo "$chunk_base"

    if [[ "$chunk_base" == "ToRun_chunk_00.txt" ]]; then
        sbatch --time=05:00:00 --array=1-${n}%100 1_snpad.sh "$chunk"
    else
        sbatch --time=03:00:00 --array=1-${n}%100 1_snpad.sh "$chunk"
    fi
done


3,19-26,46-52,70-75,102-103,130,177


```

Check summary stats after:

```bash
# Print out summary stats per scaffold and indivudal
cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/1_VCFs

shopt -s nullglob

OUT="counts_per_individual_scaf.tsv"

echo -e "Sample\tScaffold\tTotal_sites\tRefHom\tNonRefHom\tnHets" > "$OUT"

for num in scaffold_*; do
    scaffold="scaffold_${num}"

    for f in scaffold_${num}/*stats ; do
        file=$(basename "$f")

        sample=$(echo "$file" | sed "s/_scaffold_${num}.*stats//")
        
        refhomcount=$(grep "PSC" "$f" | tail -n 1 | awk -F'\t' '{print $4}')
        nonrefhomcount=$(grep "PSC" "$f" | tail -n 1 | awk -F'\t' '{print $5}')
        het=$(grep "PSC" "$f" | tail -n 1 | awk -F'\t' '{print $6}')

        total=$((refhomcount + nonrefhomcount + het))

        echo -e "${sample}\t${scaffold}\t${total}\t${refhomcount}\t${nonrefhomcount}\t${het}" >> "$OUT"

    done

    echo "Finished ${scaffold}"

done
 # Quickly check how many VCFs have less than 10 SNPs (i.e. are not useful)
 awk '$6 < 10 {print}' counts_per_individual_scaf.tsv | head

 awk -F'\t' 'NR>1 && $6 > 10 {print $2} ' counts_per_individual_scaf.tsv  | uniq

 # Sum of SNPs per scaffold across indivudals
 awk 'NR>1 {sum[$2] += $5 + $6} END {for (s in sum) print s, sum[s]}' counts_per_individual_scaf.tsv | sort -V
 

# Now sum the scaffolds per indivudal
awk -F'\t' '
BEGIN { OFS="\t" }
NR==1 {
    print "Sample", "Total_sites", "RefHom", "NonRefHom", "nHets", "Variant_sites"
    next
}
{
    total[$1]+=$3
    refhom[$1]+=$4
    nonrefhom[$1]+=$5
    hets[$1]+=$6
}
END {
    for (s in total) {
        print s, total[s], refhom[s], nonrefhom[s], hets[s], nonrefhom[s] + hets[s]
    }
}
' "$OUT" | sort > bcftools_summary_per_sample.tsv



# MORE DETAILED VERSION
# Print out summary stats per scaffold and indivudal from the bcftools stats file
cd .

OUT="bcftools_summary_per_scaffold.tsv"
:>$OUT  
shopt -s nullglob

echo -e "Scaffold\tFile\tSample\tCallable_sites\tRefHom\tNonRefHom\tHets\tMissing\tTransitions\tTransversions\tTsTv\tIndels\tAverageDepth\tSingletons\tHapRef\tHapAlt" > "$OUT"

shopt -s nullglob

for dir in $(printf '%s\n' scaffold_* | sort -V); do
    scaffold=$(basename "$dir")
    for f in "$dir"/*.stats; do echo "Processing: $f" >&2
        awk -F'\t' -v scaffold="$scaffold" -v file="$(basename "$f")" '
        /^PSC/ {
            total = $4 + $5 + $6
            if ($8 > 0)
                tstv = $7 / $8
            else
                tstv = "NA"
            print scaffold"\t"file"\t"$3"\t"total"\t"$4"\t"$5"\t"$6"\t"$14"\t"$7"\t"$8"\t"tstv"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13
        }
        ' "$f" >> "$OUT"
    done
done

# To summarise scaffolds per individual
IN="bcftools_summary_per_scaffold.tsv"
OUT="summary_stats_per_individual.tsv"

awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR == 1 {
    next
}

{
    sample = $3

    n_scaffolds[sample]++

    callable[sample] += $4
    refhom[sample]   += $5
    nonrefhom[sample] += $6
    hets[sample]     += $7
    missing[sample]  += $8
    transitions[sample] += $9
    transversions[sample] += $10
    indels[sample]   += $12
    depth_weighted_sum[sample] += $13 * $4
    singletons[sample] += $14
}

END {
    print "Sample",
          "nScaffolds",
          "Total_callable_sites",
          "Total_RefHom",
          "Total_NonRefHom",
          "Total_Hets",
          "Total_missing",
          "Total_transitions",
          "Total_transversions",
          "TsTv",
          "Total_indels",
          "Mean_depth_across_scaffolds",
          "Total_singletons",
          "Variant_sites",
          "Het_per_callable"

    for (s in callable) {

        if (transversions[s] > 0)
            tstv = transitions[s] / transversions[s]
        else
            tstv = "NA"

        if (callable[s] > 0)
            het_per_callable = hets[s] / callable[s]
        else
            het_per_callable = "NA"

        mean_depth = depth_weighted_sum[s] / callable[s]

        print s,
              n_scaffolds[s],
              callable[s],
              refhom[s],
              nonrefhom[s],
              hets[s],
              missing[s],
              transitions[s],
              transversions[s],
              tstv,
              indels[s],
              mean_depth,
              singletons[s],
              nonrefhom[s] + hets[s],
              het_per_callable
    }
}
' "$IN" | {
    read -r header
    printf "%s\n" "$header"
    sort -t$'\t' -k1,1
} > "$OUT"
```

### 2.2 Filtering gVCF

```bash
#!/bin/bash
#SBATCH --job-name=2_Filter
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21

timestamp() { date '+%F %T'; }

INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/1_VCFs"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/2_FilterGVCFs"
mkdir -p $OUTDIR

TORUN=$1

# Get sample and chrom
read -r chrom sample < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TORUN")

# Get maxDP
maxDP=14

# Ensure subdir exists
mkdir -p "${OUTDIR}/${chrom}"

INVCF="${INDIR}/${chrom}/${sample}_${chrom}.g.vcf.gz"
OUTVCF="${OUTDIR}/${chrom}/${sample}_${chrom}.filt.g.vcf.gz"

echo "$(timestamp) Filtering ${sample} ${chrom}"
echo "Input: $INVCF"
echo "Max depth: ${maxDP}"

# Filter
bcftools +setGT "${INVCF}" -Ou -- -t q -n '.' -i "FMT/DP<4 || FMT/GQ<30 || FMT/DP>${maxDP} " | bcftools view -Oz -o "$OUTVCF"
#bcftools filter -e "FMT/DP<4 || FMT/GQ<30 || FMT/DP>${maxDP} " "${INVCF}" -Oz -o "$TMPVCF"

bcftools index -f "$OUTVCF"
bcftools stats -s - "$OUTVCF" > "${OUTVCF%.gz}.stats"

echo "$(timestamp) Finished filtering ${sample} ${chrom}"

```

To submit the previous script:

```bash
cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/


for chunk in ToRun_chunk_*.txt
do
    n=$(wc -l < "$chunk")

    chunk_base=$(basename "$chunk")

    echo "$chunk_base"

    if [[ "$chunk_base" == "ToRun_chunk_00.txt" ]]; then
        sbatch --dependency=afterok:45174963 --time=04:00:00 --array=1-${n}%100 2_Filter.sh "$chunk"
    else
        sbatch --time=02:00:00 --array=1-${n}%100 2_Filter.sh "$chunk"
    fi
done
```

Then run the summary statistics script. 

### 2.3 Merge individuals per chromosome

Now we have to merge the downsampled inviduals and regular individuals.

To submit, simply write `sbatch --array=1-nlength <script>`

```bash
#!/bin/bash
#SBATCH --job-name=3_Merge
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21
module load libdeflate/1.21 htslib

timestamp() { date '+%F %T'; }

# === Directories ===
DIR_DOWN_VCF="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/2_FilterGVCFs"
DIR_NORMAL_VCF="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/2_FiltGVCFs"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/3_MergeGVCFs_perChrm"
mkdir -p ${OUTDIR}

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHROM_LIST )


echo "$(timestamp) Merging gVCFs for $chrom ..."

# === Create VCF list ===
VCF_LIST="${OUTDIR}/vcfs_merge.${chrom}.list"
: > "${VCF_LIST}"

missing=0

while IFS= read -r sample || [[ -n "${sample}" ]]; do
    [[ -z "${sample}" ]] && continue

    v_down="${DIR_DOWN_VCF}/${chrom}/${sample}_${chrom}.filt.g.vcf.gz"
    v_normal="${DIR_NORMAL_VCF}/${chrom}/${sample}_${chrom}.filt.gvcf.gz"

    if [[ -s "${v_down}" ]]; then
        v="${v_down}"
        echo "$(timestamp) Using downsampled VCF for ${sample}: ${v}"
    elif [[ -s "${v_normal}" ]]; then
        v="${v_normal}"
        echo "$(timestamp) Using normal VCF for ${sample}: ${v}"
    else
        echo "WARNING: Missing both downsampled and normal VCF for ${sample} on ${chrom}"
        echo "  Checked: ${v_down}"
        echo "  Checked: ${v_normal}"
        missing=$((missing + 1))
        continue
    fi

    if [[ ! -s "${v}.csi" && ! -s "${v}.tbi" ]]; then
        echo "$(timestamp) Indexing missing index for ${v}"
        bcftools index -f "${v}"
    fi

    echo "${v}" >> "${VCF_LIST}"

done < "${LIST_SELECTED_SAMPLES}"

n_vcfs=$(wc -l < "${VCF_LIST}")

if [[ "${n_vcfs}" -eq 0 ]]; then
    echo "ERROR: No valid VCFs found for ${chrom}"
    exit 0
fi

echo "$(timestamp) Found ${n_vcfs} VCFs for ${chrom}; missing ${missing}"

# === Merge ===
OUTVCF="${OUTDIR}/OBP_merged_down6x.26.${chrom}.filt.g.vcf.gz"

echo "$(timestamp) Merging starting"

bcftools merge \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -l "${VCF_LIST}" \
    -Ou \
| bcftools +fill-tags \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "${OUTVCF}" \
    -- -t AC,AN,MAF,F_MISSING


bcftools index -f "${OUTVCF}"
bcftools stats -s - "${OUTVCF}" > "${OUTVCF%.gz}.stats"

echo "$(timestamp) Done: ${OUTVCF}"
```

Then run the summary statistics script. 

### 2.4 Missingess filter

It is best to apply a per-group filters rather than a global 80% filter.

```bash
#!/bin/bash
#SBATCH --job-name=4_Missingness_filter
#SBATCH --cpus-per-task=4
#SBATCH --time=07:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1-357%100  # set to number of lines in CHROM_LIST
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21 bedtools/2.30.0

timestamp() { date '+%F %T'; }

# =========================
# Parameters
# =========================
MISS_THRESH=0.2   # F_MISSING > this fails, i.e. requires >=80% call rate per group

# =========================
# Directories / inputs
# =========================
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/3_MergeGVCFs_perChrm"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/4_GroupMissFiltered"
TMPBASE="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/4_GroupMissFiltered/tmp"
BED_DIR="${OUTDIR}/beds"

mkdir -p "${OUTDIR}"

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

MODERN_PREFIX="O-"
HIST_PREFIX="OBP"

MODERN_LIST="${OUTDIR}/Modern.txt"
HIST_LIST="${OUTDIR}/Historical.txt"

mkdir -p "${OUTDIR}" "${TMPBASE}" "${BED_DIR}"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_LIST}")
[[ -n "${chrom}" ]] || { echo "No chrom for array index ${SLURM_ARRAY_TASK_ID}"; exit 0; }


IN_VCF="${INDIR}/OBP_merged_down6x.26.${chrom}.filt.g.vcf.gz"
[[ -s "${IN_VCF}" ]] || { echo "ERROR: missing input ${IN_VCF}" >&2; exit 1; }

# =========================
# 0. RUN BEFORE SUBMITTING SCRIPT
# =========================
#        echo "$(timestamp) Generating group lists from ${IN_VCF} sample names"
#        bcftools query -l "${IN_VCF}" > "${OUTDIR}/all_samples.txt"
#        grep "^${HIST_PREFIX}"   "${OUTDIR}/all_samples.txt" > "${HIST_LIST}"
#        grep "^${MODERN_PREFIX}" "${OUTDIR}/all_samples.txt" > "${MODERN_LIST}"

#        n_all=$(wc -l < "${OUTDIR}/all_samples.txt")
#        n_grouped=$(( $(wc -l < "${HIST_LIST}") + $(wc -l < "${MODERN_LIST}") ))
#        if [[ "${n_grouped}" -ne "${n_all}" ]]; then
#            echo "WARNING: ${n_all} total samples but only ${n_grouped} matched a prefix -- check sample naming:" >&2
#            comm -23 <(sort "${OUTDIR}/all_samples.txt") <(cat "${HIST_LIST}" "${MODERN_LIST}" | sort)
#        fi

[[ -s "${MODERN_LIST}" ]] || { echo "ERROR: ${MODERN_LIST} empty/missing after generation" >&2; exit 1; }
[[ -s "${HIST_LIST}" ]]   || { echo "ERROR: ${HIST_LIST} empty/missing after generation" >&2; exit 1; }

TMPDIR="${TMPBASE}/${chrom}"
rm -rf "${TMPDIR}"
mkdir -p "${TMPDIR}"


OUT_VCF="${OUTDIR}/OBP_merged_down6x.26.${chrom}.groupmiss80.vcf.gz"
OUT_BED="${BED_DIR}/${chrom}.shared_groupmiss80.bed"

echo "$(timestamp) chrom=${chrom} MISS_THRESH=${MISS_THRESH}"
echo "IN_VCF : ${IN_VCF}"
echo "OUT_VCF: ${OUT_VCF}"


# =========================
# 1. Per-group F_MISSING -> passing sites written straight to a merged BED
# =========================
for GRP in Modern Historical; do
    if [[ "${GRP}" == "Modern" ]]; then LIST="${MODERN_LIST}"; else LIST="${HIST_LIST}"; fi

    bcftools view -S "${LIST}" --force-samples -Ou "${IN_VCF}" \
      | bcftools +fill-tags --threads "${SLURM_CPUS_PER_TASK}" -Ou -- -t F_MISSING \
      | bcftools view -i "F_MISSING<= ${MISS_THRESH}" \
      | bcftools query -f '%CHROM\t%POS0\t%POS\n' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - \
      > "${TMPDIR}/${GRP}.pass.bed"

    n_int=$(wc -l < "${TMPDIR}/${GRP}.pass.bed")
    echo "$(timestamp) ${GRP}: ${n_int} merged intervals pass <= ${MISS_THRESH} missingness"
done

# =========================
# 2. Intersect: regions passing in BOTH groups (this IS the final mask)
# =========================
bedtools intersect \
    -a "${TMPDIR}/Historical.pass.bed" \
    -b "${TMPDIR}/Modern.pass.bed" \
    > "${OUT_BED}"

n_shared=$(wc -l < "${OUT_BED}")
echo "$(timestamp) ${n_shared} shared intervals pass in both groups -> ${OUT_BED}"

if [[ "${n_shared}" -eq 0 ]]; then
    echo "WARNING: no shared passing regions for ${chrom}"
exit 0
fi

# =========================
# 3. Apply shared BED mask back to the full, all-sample VCF
#    (bcftools recognises the .bed extension as 0-based, half-open)
# =========================
bcftools view -T "${OUT_BED}" -Oz -o "${OUT_VCF}" "${IN_VCF}"
bcftools index -f "${OUT_VCF}"
bcftools stats -s - "${OUT_VCF}" > "${OUT_VCF%.gz}.stats"


echo "$(timestamp) Done: ${OUT_VCF}"
```

##### OLD - Global 80% filter

```bash
#!/bin/bash
#SBATCH --job-name=4_Miss_filter
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --array=1-357%100
#SBATCH --output=jobs/%x_%A_%a.log
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21
module load libdeflate/1.21 htslib

timestamp() { date '+%F %T'; }

# === Directories ===
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/3_MergeGVCFs_perChrm"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/4_stricter_tests"
mkdir -p "${OUTDIR}"

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

chrom=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CHROM_LIST}")

INVCF="${INDIR}/OBP_merged_down6x.26.${chrom}.filt.g.vcf.gz"

if [[ ! -s "${INVCF}" ]]; then
    echo "ERROR: Input VCF not found: ${INVCF}"
    exit 1
fi

echo "$(timestamp) Filtering ${chrom}"

# === 80% called filter ===
# For 26 samples, this means that there should be at least 20 samples 

OUT80="${OUTDIR}/OBP_merged_down6x.26.${chrom}.miss80.vcf.gz"

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'F_MISSING<=0.20' \
    -Oz \
    -o "${OUT80}" \
    "${INVCF}"

bcftools index -f "${OUT80}"
bcftools stats -s - "${OUT80}" > "${OUT80%.gz}.stats"


# === 100% called filter ===
# F_MISSING == 0 means no missing genotypes across all samples.

OUT100="${OUTDIR}/OBP_merged_down6x.26.${chrom}.miss100.vcf.gz"

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'F_MISSING==0' \
    -Oz \
    -o "${OUT100}" \
    "${INVCF}"

bcftools index -f "${OUT100}"
bcftools stats -s - "${OUT100}" > "${OUT100%.gz}.stats"


# === Quick site counts ===
echo -e "Chromosome\tInput_sites\tMissing90_sites\tMissing100_sites" > "${OUTDIR}/missingness_site_counts.tsv"

echo -e "${chrom}\t$(bcftools view -H "${INVCF}" | wc -l)\t$(bcftools view -H "${OUT80}" | wc -l)\t$(bcftools view -H "${OUT100}" | wc -l)" \
    >> "${OUTDIR}/missingness_site_counts.tsv"

echo "$(timestamp) Done filtering ${chrom}"
echo "$(timestamp) 80% VCF:  ${OUT80}"
echo "$(timestamp) 100% VCF: ${OUT100}"
```

### 2.5 Bisnps

```bash
#!/bin/bash
#SBATCH --job-name=5_bisnps_nosingletons
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --output=jobs/%x_%j.log
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.21 libdeflate/1.21 htslib

# === Directories ===
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/4_GroupMissFiltered"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/5_bisnps/miss80"
mkdir -p "${OUTDIR}"

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

cat ${CHROM_LIST} | while read chrom; do 
    [[ -n "${chrom}" ]] || continue
 
    # Check if the group-missingness-filtered VCF for this scaffold exists
    
INVCF="${INDIR}/OBP_merged_down6x.26.${chrom}.groupmiss80.vcf.gz"
    if [[ -s "${INVCF}" ]]; then

OUTVCF="${OUTDIR}/OBP_merged_down6x.26.${chrom}.groupmiss80.bisnps.vcf.gz"
OUTVCF_NOSINGLETONS="${OUTDIR}/OBP_merged_down6x.26.${chrom}.miss80.bisnps.MAC2.vcf.gz"


echo "$(timestamp) Extracting biallelic snps for $chrom ..."

# Extract biallelic SNPs and make sure to only keep sites where there is at least one individual with a SNP
        bcftools view --threads "$SLURM_CPUS_PER_TASK" -m2 -M2 -v snps "${INVCF}" -Ou \
          | bcftools +fill-tags --threads "${SLURM_CPUS_PER_TASK}" -Ou -- -t AC,AN,MAF,F_MISSING \
          | bcftools view --threads "${SLURM_CPUS_PER_TASK}" -e 'INFO/AC[0]=0' -Oz -o "${OUTVCF}"
        bcftools index -f "${OUTVCF}"
        bcftools stats -s - "${OUTVCF}" > "${OUTVCF%.gz}.stats"

echo "$(timestamp) Creating VCF with no singletons "

# MAC2 (biallelic; AC[0] is ALT allele count)
bcftools view --threads "${SLURM_CPUS_PER_TASK}" -i 'INFO/AC[0] > 1' -Oz -o "${OUTVCF_NOSINGLETONS}" "${OUTVCF}"
bcftools index -f "${OUTVCF_NOSINGLETONS}"
bcftools stats -s - "${OUTVCF_NOSINGLETONS}" > "${OUTVCF_NOSINGLETONS%.gz}.stats"

echo "$(timestamp) Finished extracting bisnps for ${chrom}"
    else
        echo "$(timestamp) ${INVCF} was SKIPPED because it was not created (did not pass the missingness filter)."
    fi
 
done
```

### 2.6 Concat 

```bash
#!/bin/bash 
#SBATCH --job-name=6_concat
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --output=jobs/%x_%j.log
#SBATCH --error=jobs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge
module load perl/5.38.0 gsl/2.5 bcftools/1.22
module load libdeflate/1.21 htslib
module load vcftools

# === Directories and files ===
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/5_bisnps/miss80"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x/6_concat"
mkdir -p $OUTDIR

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

OUTVCF="${OUTDIR}/OBP_merged_down6x.26.concat.groupmiss80.bisnps.vcf.gz"


echo "[$(timestamp)] Concatenating scaffold VCFs..."

###############################################################################
# Step 1: Build ordered list of per-scaffold VCFs
###############################################################################
VCF_LIST=()
MISSING=0

while IFS= read -r CHRM; do
    CHRM=${CHRM%$'\r'}

    # Skip blank lines and comments
    [[ -z "${CHRM//[[:space:]]/}" || "${CHRM}" =~ ^[[:space:]]*# ]] && continue

    INVCF="${INDIR}/OBP_merged_down6x.26.${CHRM}.groupmiss80.bisnps.vcf.gz"

    if [[ -f "${INVCF}" ]]; then
        VCF_LIST+=("${INVCF}")
    else
        echo "ERROR: Missing VCF: ${INVCF}" >&2
        ((MISSING += 1))
    fi
done < "${CHROM_LIST}"


if (( ${#VCF_LIST[@]} == 0 )); then
    echo "ERROR: No input VCFs were found." >&2
    exit 1
fi

echo "[$(timestamp)] Found ${#VCF_LIST[@]} scaffold VCFs."

###############################################################################
# Step 2: Concatenate scaffold VCFs
###############################################################################

bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${VCF_LIST[@]}" \
    -Oz \
    -o "${OUTVCF}"

###############################################################################
# Step 3: Index and calculate summary statistics
###############################################################################

bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    "${OUTVCF}"

bcftools stats \
    -s - \
    "${OUTVCF}" \
    > "${OUTVCF%.vcf.gz}.stats"

echo "[$(timestamp)] Finished concatenating biallelic SNP VCFs."
echo "Output: ${OUTVCF}"
```


### 2.7 Summary stats

After running the pipeline, we can extract the summary stats from `bcftools stats`. 
```bash
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

BASE="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/Down6x"
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
collect_stage "2_Filtered_gVCF"          "${BASE}/2_FilterGVCFs/*/*.stats"
echo "Finished 2_Filtered_gVCF"
collect_stage "3_Merged"                 "${BASE}/3_MergeGVCFs_perChrm/*.filt.g.vcf.stats"
echo "Finished 3_Merged"
collect_stage "4_GroupMiss80"            "${BASE}/4_GroupMissFiltered/*.groupmiss80.vcf.stats"
echo "Finished 4_GroupMiss80"
collect_stage "5_Bisnp_with_singletons"  "${BASE}/5_bisnps/miss80/*.groupmiss80.bisnps.vcf.stats"
echo "Finished 5_Bisnp_with_singletons"
collect_stage "6_Bisnp_no_singletons"    "${BASE}/5_bisnps/miss80/*.miss80.bisnps.MAC2.vcf.stats"
echo "Finished 6_Bisnp_no_singletons"
collect_stage "7_Concat_genome_wide"     "${BASE}/6_concat/*.concat.groupmiss80.bisnps.stats"
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

```

Calculating depth percentiles from VCFs and BAMs

```bash
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

```

### Pipeline Results
After collecting the summary stats, I put the data into excel and R to visulaise the main trends or unusual results:
> 

