#!/bin/bash
#SBATCH --job-name=5_bisnps
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:006
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
INDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/4_GroupMissFiltered"
OUTDIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/5_bisnps/miss80"
mkdir -p "${OUTDIR}"

LIST_SELECTED_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
CHROM_LIST="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

cat ${CHROM_LIST} | while read chrom; do 
    [[ -n "${chrom}" ]] || continue
 
    # Check if the group-missingness-filtered VCF for this scaffold exists
    INVCF="${INDIR}/OBP_merged.26.${chrom}.groupmiss80.vcf.gz"
 
    if [[ -s "${INVCF}" ]]; then
        OUTVCF="${OUTDIR}/OBP_merged.26.${chrom}.groupmiss80.bisnps.vcf.gz"
   OUTVCF_NOSINGLETONS="${OUTDIR}/OBP_merged.26.${chrom}.groupmiss80.bisnps.MAC2.vcf.gz"
 
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
