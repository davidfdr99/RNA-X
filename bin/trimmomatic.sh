#!/bin/bash

# Trims raw fastq files using trimmomatic
# Options include the adapter and the directory of output

sample_id=${1}
reads=${2}
trimdir=${3}/trimmed_reads
adapter=${4}

mkdir "$trimdir"

trimmomatic PE -threads 4 "$reads"\
     "$trimdir"/paired_"$sample_id"_1.fastqc.gz "$trimdir"/unpaired_"$sample_id"_1.fastq.gz\
     "$trimdir"/paired_"$sample_id"_2.fastqc.gz "$trimdir"/unpaired_"$sample_id"_2.fastq.gz\
     ILLUMINACLIP:"$adapter":2:30:10\
     SLIDINGWINDOW:4:15 MINLEN:36
