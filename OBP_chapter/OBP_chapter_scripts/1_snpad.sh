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
