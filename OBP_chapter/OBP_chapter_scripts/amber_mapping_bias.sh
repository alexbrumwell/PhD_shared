# Create BAM list file has two columns: first sample ID, second file path
module load python/3.12.8

DIR="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP"
DIR_OUTPUT="/projects/erode/people/qhk515/01_OBP/data/raw_data/bams_OBP/AMBER_QC"
DIR_AMBER="/projects/erode/people/qhk515/software/AMBER"
mkdir -p $DIR_OUTPUT

cd $DIR

# Make a list of 6 samples
paste <(basename -a *bam ) <(realpath *bam) | tail -n 6 > ${DIR_OUTPUT}/BAM_AMBER.list

# Make a list of scaffolds to exclude
seq 2 2002 | sed 's/^/scaffold_/' > ${DIR_OUTPUT}/scaffolds.list

# RUN AMBER 
cd $DIR_AMBER
python -m amber --bamfiles ${DIR_OUTPUT}/BAM_AMBER.list --output ${DIR_OUTPUT}/AMBER_QC --exclude ${DIR_OUTPUT}/scaffolds.list
