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
```

You can also run the pipeline with docker:

```bash
$ docker pull davidfdr/rnax
```

