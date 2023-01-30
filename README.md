# RNA-X

RNA-X is a nextflow implementation of a simple bulk RNA-seq pipeline using fastqc, trimmomatic and Kallisto.

## Quick Start

If not done already, install the Nextflow runtime:

```bash
$ curl -fsSL get.nextflow.io | bash
```

It is recommended to install all the required dependencies with conda or to use Docker. The former will install all the packages locally:

```bash
$ conda env create -f environment.yml
$ conda activate rnax
$ nextflow run davidfdr/rnax -profile test
```

You can also run the pipeline with docker:

```bash
$ docker pull davidfdr/rnax
$ docker run davidfdr/rnax -profile test -with-docker
```

## Test Data

The test data set provided is from one sample of *S. cerevisiae* that can be found under the accession SRS307298 of the SRA database. There are two different biological conditions (Batch and CENPK) and three technical replicates respectively. 

The data can be downloaded locally using samtools. A script is provided that does that job:

```bash
$ cd data/test_reads
$ bash fetch.sh SRR_Acc_list.txt
```


