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

log.info """\
       R N A - X    P I P E L I N E
    ===================================
    reads        : ${params.reads}
    transcriptome: ${params.transcriptome_fasta}
    outdir       : ${params.outdir}
    adapter	 : ${params.adapter}
    annotation   : ${params.annotation}
    """
    .stripIndent()

/* fastqc process
 */
process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/$sample_id", mode: 'copy'
    container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"


    input:
    tuple val(sample_id), path(reads)
    tuple val(sample_id), path(trimmed_read_1), path(trimmed_read_2)

    output:
    path "fastqc_${sample_id}_untrimmed"
    path "fastqc_${sample_id}_trimmed"

    script:
    """
    fastqc.sh "fastqc_${sample_id}_untrimmed" "${reads}"
    fastqc.sh "fastqc_${sample_id}_trimmed" "${trimmed_read_1} ${trimmed_read_2}"
    """
}

/* MultiQC upon completion
 */
process MULTIQC {
    publishDir params.outdir, mode:'copy'
    container "https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0"

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
    label 'process_medium'

    tag "Trimming on $sample_id"
    container 'https://depot.galaxyproject.org/singularity/trimmomatic:0.35--hdfd78af_7'

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
    container 'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2'

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
    label 'process_high'

    tag "Kallisto on $sample_id"
    publishDir "$params.outdir/kallisto", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2'

    input:
    path kallisto_index
    tuple val(sample_id), path(read_1), path(read_2)

    output:
    path "kallisto_${sample_id}"

    script:
    """
    kallisto quant -i $kallisto_index -o "kallisto_${sample_id}" -t 4 --genomebam --gtf ${params.annotation} --chromosomes ${params.chromosome_file} ${read_1} ${read_2}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    trimming_ch = TRIMMING(read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch, trimming_ch)
    index_ch = INDEX(params.transcriptome_fasta)
    quant_ch = QUANTIFICATION(index_ch, trimming_ch)
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
    
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html : "Oops .. something went wrong" )
}
