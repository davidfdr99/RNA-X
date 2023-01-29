#!/bin/bash

# Trims raw fastq files using trimmomatic
# Options include the adapter and the directory of output

sample_id=${1}
reads=${2}
adapter=${3}

trimmomatic PE -threads 4 ${reads}\
     paired_"$sample_id"_1.fastq.gz unpaired_"$sample_id"_1.fastq.gz\
     paired_"$sample_id"_2.fastq.gz unpaired_"$sample_id"_2.fastq.gz\
     ILLUMINACLIP:"$adapter":2:30:10\
     SLIDINGWINDOW:4:15 MINLEN:36
