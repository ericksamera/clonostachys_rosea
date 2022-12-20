#!/bin/bash
SAMPLES_LIST=$(<$"samples_list.txt")
OUTPUT_FILE="log-fix_reads.txt"

#input_dir
INPUT_DIR="ena_files/"
FIXED_FASTQ="fixed_fastq/"
REFERENCE='clonostachys_rosea'

for SAMPLE_PATH in $SAMPLES_LIST
do
    echo "[LOG] ==================[$SAMPLE_PATH]=================" &>> $OUTPUT_FILE
    echo "[LOG] `date`: [$SAMPLE_PATH]: REPAIRING ..." &>> $OUTPUT_FILE
    (repair.sh\
    in1=$INPUT_DIR$SAMPLE_PATH"_1.fastq.gz"\
    in2=$INPUT_DIR$SAMPLE_PATH"_2.fastq.gz"\
    out1=$FIXED_FASTQ$SAMPLE_PATH"_1.fastq.gz"\
    out2=$FIXED_FASTQ$SAMPLE_PATH"_2.fastq.gz"\
    outs=$FIXED_FASTQ$SAMPLE_PATH"_singletons.fastq.gz"\
    repair=t\
    nullifybrokenquality=t) &>> $OUTPUT_FILE
done

    echo "[LOG] `date`: COMPLETED ALL SAMPLES. " &>> $OUTPUT_FILE