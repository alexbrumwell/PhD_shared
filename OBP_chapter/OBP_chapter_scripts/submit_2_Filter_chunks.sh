cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/2_FiltGVCFs/

for chunk in /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/ToRun_chunk_*.txt
do
    n=$(wc -l < "$chunk")

    chunk_base=$(basename "$chunk")

    echo "$chunk_base"

    if [[ "$chunk_base" == "ToRun_chunk_00.txt" ]]; then
        sbatch --array=1-${n}%100 2_Filter.sh "$chunk"
    else
        sbatch --time=01:00:00 --array=1-${n}%100 2_Filter.sh "$chunk"
    fi
done
