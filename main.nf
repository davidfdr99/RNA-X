/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/test/*_{1,2}.fastq.gz"

params.multiqc = "./multiqc"
params.outdir = "./results"
params.adapter="$projectDir/assets/trimmomatic/adapter/TruSeq2-PE.fa"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/* fasted process
 */

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    fastqc.sh "$sample_id" "$reads"
    """
}


process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path '"$outdir"/multiqc_report.html'

    script:
    """
    multiqc .
    """
}


process TRIMMING {
    tag "Trimming on $sample_id"

    input: 
    tuple val(sample_id), path(reads)

    output:
    path '"$outdir"/trimmed_reads'

    script:
    """
    trimmomatic.sh "$sample_id" "$reads" "$outdir" "$adapter"
    """
}

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file

process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}
*/

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

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
