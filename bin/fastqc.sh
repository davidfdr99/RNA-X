#!/bin/bash

# Performs fastqc on input reads 

set -u

outdir=${1}
reads=${2}

mkdir ${outdir}
fastqc -o ${outdir} -f fastq -q ${reads}
