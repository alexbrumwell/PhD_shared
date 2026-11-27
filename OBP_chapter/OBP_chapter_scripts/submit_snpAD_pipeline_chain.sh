# SCRIPT TO CHAIN SCRIPTS FOR PIPELINE 
DIR="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output"

chunk="/projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/ToRun_chunk_00.txt"

# Submit script 1
n=$(wc -l < "$chunk")
cd ${DIR}/2_FiltGVCFs
job1=$(sbatch --parsable --array=1-${n}%100 2_Filter.sh $chunk )
echo "Submitted script1.sh as job ${job1}"

# Submit script 2 after script 1 finishes successfully
cd ${DIR}/3_MergeGVCFs_perChrm
job2=$(sbatch --parsable --dependency=afterok:${job1} --array=1-357%100 3_Merge.sh)
echo "Submitted script2.sh as job ${job2}"

# Submit script 3 after script 2 finishes successfully
cd ${DIR}/4_stricter_tests
job3=$(sbatch --parsable --dependency=afterok:${job2} --array=1-357%100 4_Missingness_filter.sh)
echo "Submitted script3.sh as job ${job3}"

# Submit script 4 after script 3 finishes successfully
cd ${DIR}/5_bisnps
job4=$(sbatch --parsable --dependency=afterok:${job3} --array=1-357%100 5_bisnps_nosingletons.sh)
echo "Submitted script4.sh as job ${job4}"

# Submit script 5 after script 4 finishes successfully
cd ${DIR}/6_concat
job5=$(sbatch --parsable --dependency=afterok:${job4} 6_concat.sh)
echo "Submitted script5.sh as job ${job5}"
