#!/bin/bash
#SBATCH --job-name=tests_scaf1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%10
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

set -euo pipefail
timestamp() { date '+%F %T'; }

# === Load modules ===
module purge; module load angsd/0.940 winsfs/0.7.0 bcftools/1.22

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/output_angsd/hetsites_perind"
mkdir -p ${DIR}
cd ${DIR}

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_OBPs_combined.txt"
INDIR_GVCFS="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/1_VCFs"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"



# Run ANGSD per indivudal
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BAM_LIST_COMBINED})
if [ -z "${SAMPLE}" ]; then
    echo "$(timestamp) [ERROR] No sample found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
NAME=$(basename "${SAMPLE}" .rmdup.bam)

echo "$(timestamp) [*] Starting ANGSD for individual ${NAME}"


# Filtered ANGSD - # Requierements (at least >4x coverage and <50x coverage per site to be included in the analysis)
angsd -i ${SAMPLE} -ref ${REF} -r scaffold_1 -anc ${REF} \
    -sites ${LUKE_LIST_SCAF1} -out ${NAME}_nosex_filtlist \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 \
    -baq 0 -minMapQ 30 -minQ 20 -setMinDepth 4 -setMaxDepth 50 -doCounts 1 -dumpCounts 3 \
    -nThreads ${SLURM_CPUS_PER_TASK} -GL 2 -doSaf 1 \
    -doGlf 2 -doMajorMinor 1 -doMaf 2 \
    -dopost 1 -dobcf 1 --ignore-RG 0 -dogeno 1


echo "$(timestamp) [OK] ANGSD completed for ${NAME}"



exit 0

##################################################################################################################

# PART 2:
#!/bin/bash
#SBATCH --job-name=compare_angsd_snpad_hets
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --qos=normal
#SBATCH --partition=cpuqueue
#SBATCH --array=1-35%10
#SBATCH --output=jobs/%x_%A_%a.out
#SBATCH --error=jobs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amb230@kent.ac.uk

# I want to create a table that compares het sites called by ANGSD to the genotypes at those sites in the SNPAD VCF. This will help me understand what happened to those sites in the SNPAD calling process - were they retained as hets, called as hom ref, hom alt, missing, or absent from the VCF?

set -euo pipefail

timestamp() { date '+%F %T'; }

echo "$(timestamp) Starting array task ${SLURM_ARRAY_TASK_ID}"

module purge
module load bcftools/1.22

# -----------------------------
# Input/output directories
# -----------------------------
BEAGLE_DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/output_angsd/hetsites_perind"
INDIR_GVCFS="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/1_VCFs/scaffold_1"
OUTDIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/output_angsd/hetsites_perind"

mkdir -p "${OUTDIR}/jobs"
cd "${OUTDIR}"

# -----------------------------
# Select one beagle file for this array task
# -----------------------------
BEAGLE_FILE=$(find "${BEAGLE_DIR}" -maxdepth 1 -name "*beagle.gz" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "${BEAGLE_FILE}" ]]; then
    echo "$(timestamp) ERROR: No beagle file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

name=$(basename "${BEAGLE_FILE}" "_nosex_filtlist.beagle.gz")

echo "$(timestamp) Processing sample: ${name}"
echo "$(timestamp) Beagle file: ${BEAGLE_FILE}"

VCF="${INDIR_GVCFS}/${name}_scaffold_1.vcf.gz"

if [[ ! -f "${VCF}" ]]; then
    echo "$(timestamp) ERROR: VCF not found for ${name}: ${VCF}"
    exit 1
fi

# -----------------------------
# Step 1. Extract ANGSD het sites from beagle
# -----------------------------
echo "$(timestamp) Extracting ANGSD het sites"

zcat "${BEAGLE_FILE}" | \
awk -v sample="${name}" 'BEGIN{OFS="\t"}
NR==1 {next}
{
    n = split($1, a, "_")
    pos = a[n]
    chrom = a[1]
    for (i = 2; i < n; i++) {
        chrom = chrom "_" a[i]
    }
    ANGSD_A1 = $2
    ANGSD_A2 = $3
    pAA = $4
    pAB = $5
    pBB = $6

    if (pAB > 0.95) {
        key = chrom ":" pos
        print key, sample, chrom, pos, ANGSD_A1, ANGSD_A2, pAA, pAB, pBB, "het"
    }
}' > "${name}_angsd_het_sites.noheader.tsv"

{
    echo -e "key\tsample\tchrom\tpos\tANGSD_A1\tANGSD_A2\tANGSD_pAA\tANGSD_pAB\tANGSD_pBB\tANGSD_call"
    cat "${name}_angsd_het_sites.noheader.tsv"
} > "${name}_angsd_het_sites.tsv"

awk 'BEGIN{OFS="\t"} NR>1 {print $3, $4}' \
    "${name}_angsd_het_sites.tsv" \
    > "${name}_angsd_het_sites.positions.tsv"


# -----------------------------
# Step 2. Extract information from read counts from pos and counts file 
# -----------------------------
echo "$(timestamp) Extracting SNPAD information"

bcftools view -R "${name}_angsd_het_sites.positions.tsv" "${VCF}" -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%A\t%C\t%G\t%T\t%GQ]\n' | \
awk 'BEGIN{OFS="\t"}
{
    key=$1 ":" $2
    print key, $1, $2, $3, $4, "yes", $5, $6, $7, $8, $9, $10, $11
}' > "${name}_snpad_at_angsd_hets.noheader.tsv"

{
    echo -e "key\tchrom\tpos\tREF\tALT\tSNPAD_present\tSNPAD_GT\tSNPAD_DP\tSNPAD_A_reads\tSNPAD_C_reads\tSNPAD_G_reads\tSNPAD_T_reads\tSNPAD_GQ"
    cat "${name}_snpad_at_angsd_hets.noheader.tsv"
} > "${name}_snpad_at_angsd_hets.tsv"

# STEP 3: Extract the REF and ALT allele from the MAFs file to get the correct reference and alternative allele for each site. This is necessary because the SNPAD VCF may have different REF and ALT alleles than the ANGSD output, and we need to know which one is which to compare read counts correctly.
# Match the 3rd and 4th col of the _angsd_het_sites file to the 1st and 2nd col of the mafs file to find the correct site and then extract the major, minor and ref alele from the maf file and add it to the _angsd_het_sites file as additional columns
echo "$(timestamp) Adding MAF allele information"

awk 'BEGIN{OFS="\t"}
NR==FNR {
    if (FNR==1) {
        header = $0
        next
    }
    key = $3 FS $4
    line[key] = $0
    order[++n] = key
    next
}
FNR==1 { next }
{
    key = $1 FS $2
    if (key in line) {
        maf_ref[key]   = $5
        maf_major[key] = $3
        maf_minor[key] = $4
    }
}
END {
    print header, "maf_ref", "maf_major", "maf_minor"
    for (i = 1; i <= n; i++) {
        key = order[i]
        print line[key], \
              ((key in maf_ref) ? maf_ref[key] : "NA"), \
              ((key in maf_major) ? maf_major[key] : "NA"), \
              ((key in maf_minor) ? maf_minor[key] : "NA")
    }
}
' "${name}_angsd_het_sites.tsv" <(zcat "${name}_nosex_filtlist.mafs.gz") \
> "${name}_angsd_het_sites_with_maf.tsv"


# -----------------------------
# Step 4. Add ANGSD read counts
# Memory-safe: stores only ANGSD het sites, streams through pos/counts
# -----------------------------
echo "$(timestamp) Adding ANGSD read counts"

awk 'BEGIN{OFS="\t"}
NR==FNR {
    if (FNR==1) {
        header = $0
        next
    }
    key = $3 FS $4
    line[key] = $0
    order[++n] = key
    wanted[key] = 1
    next
}

FNR==1 {
    next
}

{
    key = $1 FS $2
    if (key in wanted) {
        counts[key] = $0
    }
}

END {
    print header, "ANGSD_count_chr", "ANGSD_count_pos", "ANGSD_totDepth", "ANGSD_totA", "ANGSD_totC", "ANGSD_totG", "ANGSD_totT"

    for (i = 1; i <= n; i++) {
        key = order[i]

        if (key in counts) {
            print line[key], counts[key]
        } else {
            print line[key], "NA", "NA", "NA", "NA", "NA", "NA", "NA"
        }
    }
}
' "${name}_angsd_het_sites_with_maf.tsv" \
  <(paste <(zcat "${name}_nosex_filtlist.pos.gz") <(zcat "${name}_nosex_filtlist.counts.gz")) \
> "${name}_angsd_het_sites_with_maf_and_counts.tsv"

# -----------------------------
# Step 5. Join everything and calculate final table
# -----------------------------
echo "$(timestamp) Creating final table"

ANGSD_FILE="${name}_angsd_het_sites_with_maf_and_counts.tsv"
SNPAD_FILE="${name}_snpad_at_angsd_hets.tsv"
OUT="${name}_final_angsd_snpad_het_comparison.tsv"


awk 'BEGIN{OFS="\t"}
NR==FNR {
    if (FNR==1) {
        snpad_header=$0
        next
    }

    key=$1
    snpad[key]=$0
    next
}

FNR==1 {
    print $0, snpad_header, "outcome"
    next
}

{
    key=$1

    if (key in snpad) {
        split(snpad[key], s, OFS)

        gt=s[7]

        if (gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") {
            outcome="retained_het"
        } else if (gt=="0/0" || gt=="0|0") {
            outcome="converted_to_hom_ref"
        } else if (gt=="1/1" || gt=="1|1") {
            outcome="converted_to_hom_alt"
        } else if (gt=="./." || gt==".|." || gt==".") {
            outcome="missing_genotype"
        } else {
            outcome="other"
        }

        print $0, snpad[key], outcome
    } else {
        print $0, "NA", "NA", "NA", "NA", "NA", "no", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "absent_from_vcf"
    }
}
' "${SNPAD_FILE}" "${ANGSD_FILE}" > "${OUT}"

# clean unecessary intermediate files
rm -f "${name}_angsd_het_sites.noheader.tsv" \
      "${name}_angsd_het_sites.tsv" \
      "${name}_angsd_het_sites.positions.tsv" \
      "${name}_angsd_het_sites_with_maf.tsv" \
      "${name}_snpad_at_angsd_hets.noheader.tsv" \
      ${name}_angsd_het_sites_with_maf_and_counts.tsv \
      ${name}_snpad_at_angsd_hets.tsv


echo "$(timestamp) Completed array task ${SLURM_ARRAY_TASK_ID} for sample ${name}"

exit 0

# AFter running script compare the final tables across all individuals to see if there are any patterns in which sites were retained as hets, converted to hom ref, hom alt, missing, or absent from the VCF. This will help us understand what happened to those sites in the SNPAD calling process and whether there are any biases in which sites were retained or lost. Make sure to only include one header
awk 'BEGIN{OFS="\t"} FNR==1 && NR!=1 {next} {print}' *_final_angsd_snpad_het_comparison.tsv > all_final_samples_angsd_snpad_het_comparison.tsv




exit 0

# To extract het stites from VCFs created using ANGSD, I can use the following command:

# === Inputs ====
DIR="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/output_angsd/hetsites_perind"
mkdir -p ${DIR}
cd ${DIR}

# Variables
CALLABLE_SITES_LUKE_COMBINED="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/Luke_het/entire_genome_combined_nosex.list" #total sites 717644389 
LUKE_LIST_SCAF1="/projects/erode/people/qhk515/01_OBP/popgen_analyses/het/tests_scaf1/scaf1_callable_sites.txt.tab"
BAM_LIST_COMBINED="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_OBPs_combined.txt"
INDIR_GVCFS="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/1_VCFs"

AUTOSOMAL_SCAFFOLDS="/projects/erode/people/qhk515/01_OBP/scripts/OBPs/list_scaffolds_autosomal.txt"
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa"
REF_INDEX="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"



# Run ANGSD per indivudal
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${BAM_LIST_COMBINED})
if [ -z "${SAMPLE}" ]; then
    echo "$(timestamp) [ERROR] No sample found for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
NAME=$(basename "${SAMPLE}" .rmdup.bam)



# -----------------------------
# Step 2. Extract read counts from VCF
# -----------------------------
echo "$(timestamp) Extracting SNPAD information"

bcftools view -R "${name}_angsd_het_sites.positions.tsv" "${VCF}" -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%A\t%C\t%G\t%T\t%GQ]\n' | \
awk 'BEGIN{OFS="\t"}
{
    key=$1 ":" $2
    print key, $1, $2, $3, $4, "yes", $5, $6, $7, $8, $9, $10, $11
}' > "${name}_snpad_at_angsd_hets.noheader.tsv"

{
    echo -e "key\tchrom\tpos\tREF\tALT\tSNPAD_present\tSNPAD_GT\tSNPAD_DP\tSNPAD_A_reads\tSNPAD_C_reads\tSNPAD_G_reads\tSNPAD_T_reads\tSNPAD_GQ"
    cat "${name}_snpad_at_angsd_hets.noheader.tsv"
} > "${name}_snpad_at_angsd_hets.tsv"

