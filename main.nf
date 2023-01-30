/*
 * Copyright (c) 2023, David Fandrei.
 *
 *   This file is part of 'RNA-X', a simple pipeline for bulk RNA-seq implementing Kallisto.
 *
 *   RNA-X is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNA-X is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 */

params.reads = "$projectDir/data/test_reads/*_{1,2}.fastq.gz"
params.outdir = "$projectDir/data/test_results"
params.adapter="$projectDir/assets/trimmomatic/adapters/TruSeq3-PE.fa"
params.transcriptome_fasta="$projectDir/assets/transcriptomes/gencode.v42.transcripts.fa.gz"
params.annotation= null
params.chromosome_file="$projectDir/assets/transcriptomes/chromInfo.hg19.tsv"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    transcriptome: ${params.transcriptome_fasta}
    outdir       : ${params.outdir}
    adapter	 : ${params.adapter}
    annotation   : ${params.annotation)
    """
    .stripIndent()

/* Check if annotation file exists
 */

annotationFile = file(params.annotation)

/* fastqc process
 */
process FASTQC {
    tag "FASTQC on $sample_id"

    publishDir "${params.outdir}/$sample_id/qc_results", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    tuple val(sample_id), path(trimmed_read_1), path(trimmed_read_2)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    fastqc.sh "$sample_id" "$reads" "$trimmed_read_1" "$trimmed_read_2"
    """
}

/* MultiQC upon completion
 */
process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

/* Trimming process
 */
process TRIMMING {
    tag "Trimming on $sample_id"

    publishDir "${params.outdir}/$sample_id/trim_results"

    input: 
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("paired_${sample_id}_1.fastq.gz"), path("paired_${sample_id}_2.fastq.gz")

    script:
    """
    trimmomatic.sh "$sample_id" "$reads" $params.adapter
    """
}

/* Kallisto
 */
process INDEX {
    tag "Indexing using Kallisto"

    input:
    path transcriptome

    output:
    path 'kallisto_index'

    script:
    """
    kallisto index -i kallisto_index $transcriptome
    """
}

process QUANTIFICATION {
    tag "Kallisto on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    path kallisto_index
    tuple val(sample_id), path(read_1), path(read_2)

    output:
    path "kallisto_${sample_id}"

    script:
    if ( !annotationFile.exists() )
        """
        kallisto quant -i $kallisto_index -o "kallisto_${sample_id}" -t 4 --pseudobam ${read_1} ${read_2} 
        """
    if ( annotationFile.exists() )
        """
        kallisto quant -i $kallisto_index -o "kallisto_${sample_id}" -t 4 --genomebam --gtf $annotationFile --chromosomes $params.chromosome_file ${read_1} ${read_2}
        """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    trimming_ch = TRIMMING(read_pairs_ch)
    index_ch = INDEX(params.transcriptome_fasta)
    quant_ch = QUANTIFICATION(index_ch, trimming_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration: ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
    
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
