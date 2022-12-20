#!/bin/bash
SAMPLES_LIST=$(<$"samples_list.txt")
OUTPUT_FILE="log-generate_bams.txt"

#input_dir
OUTPUT_DIR="bam_files/"
FIXED_FASTQ="fixed_fastq/"
REFERENCE='clonostachys_rosea'

for SAMPLE_PATH in $SAMPLES_LIST
do
    echo "[LOG] `date`: [$SAMPLE_PATH]: ALIGNING ..." &>> $OUTPUT_FILE
    time (bowtie2 \
    --no-unal\
    -p 8\
    -x $REFERENCE\
    -1 $FIXED_FASTQ$SAMPLE_PATH"_1.fastq.gz"\
    -2 $FIXED_FASTQ$SAMPLE_PATH"_2.fastq.gz"\
    | samtools view -bS - > "$OUTPUT_DIR$SAMPLE_PATH.bam") &>> $OUTPUT_FILE
    echo "[LOG] `date`: [$SAMPLE_PATH]: COMPLETED. " &>> $OUTPUT_FILE
done

    echo "[LOG] `date`: COMPLETED ALL SAMPLES. " &>> $OUTPUT_FILE
