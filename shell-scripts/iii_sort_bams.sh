#!/bin/bash
SAMPLES_LIST=$(<$"samples_list.txt")
OUTPUT_FILE="log-sort_bams.txt"

#input_dir
BAM_FILES="bam_files/"
OUTPUT_FILES="bam_files_sorted/"

for SAMPLE_PATH in $SAMPLES_LIST
do
    echo "[LOG] `date`: [$SAMPLE_PATH]: ALIGNING ..." &>> $OUTPUT_FILE
    time (samtools sort\
    -@ 8\
    $BAM_FILES$SAMPLE_PATH.bam\
    -o $OUTPUT_FILES/$SAMPLE_PATH.bam) &>> $OUTPUT_FILE
    echo "[LOG] `date`: [$SAMPLE_PATH]: COMPLETED. " &>> $OUTPUT_FILE
done

echo "[LOG] `date`: COMPLETED ALL SAMPLES. " &>> $OUTPUT_FILE