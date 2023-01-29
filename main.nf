/*
 * pipeline input parameters
 */

params.reads = "$projectDir/data/test/*_{1,2}.fastq.gz"
params.multiqc = "./multiqc"
params.outdir = "$projectDir/data/test_results"
params.adapter="$projectDir/assets/trimmomatic/adapters/TruSeq3-PE.fa"
params.transcriptome_fasta="$projectDir/assets/transcriptomes/gencode.v42.transcripts.fa.gz"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    transcriptome: ${params.transcriptome_fasta}
    outdir       : ${params.outdir}
    adapter	 : ${params.adapter}
    """
    .stripIndent()

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
    path "$sample_id"

    script:
    """
    kallisto quant -i $kallisto_index -o "kallisto_${sample_id}" -t 4 ${read_1} ${read_2} 
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
