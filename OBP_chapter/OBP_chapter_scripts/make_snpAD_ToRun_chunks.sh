LIST_SAMPLES="/projects/erode/people/qhk515/01_OBP/data/OBP_samples_selected_6x.list"
LIST_SCAF="/projects/erode/people/qhk515/01_OBP/data/OBP_autosomes_ge100kb.list"

cd /projects/erode/people/qhk515/01_OBP/data/Calling_snpAD/output/

: > ToRun.txt

while IFS= read -r scaffold; do
  while IFS= read -r sample; do
    printf '%s\t%s\n' "$scaffold" "$sample" >> ToRun.txt
  done < $LIST_SAMPLES
done < $LIST_SCAF

# Filter for whichever samples and scaffolds to downsample
#grep -w "scaffold_1" ToRun.txt > tmp.txt && mv tmp.txt ToRun.txt

# Now split the text file into max. 9999 lines for SLURM submission
split -l 9999 -d --additional-suffix=.txt ToRun.txt ToRun_chunk_
