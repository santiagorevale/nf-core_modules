// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_GENOTYPEGVCFS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--0"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }

    input:
    tuple val(meta), path(gvcf), path(gvcfIndex)
    path fasta
    path fastaIndex
    path fastaDict
    path dbsnp
    path dbsnpIndex
    path intervalsBed

    output:
    tuple val(meta), path("*.genotyped.vcf.gz"), emit: vcf
    path "*.version.txt"                       , emit: version

    script:
    def software        = getSoftwareName(task.process)
    def prefix          = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def dbsnpOption     = dbsnp ? "-D ${dbsnp}" : ""
    def intervalsOption = intervalsBed ? "-L ${intervalsBed}" : ""
    def avail_mem = 6
    if (!task.memory) {
        log.info "[GATK GenotypeGVCFs] Available memory not known - defaulting to ${avail_mem}GB. Specify process memory requirements to change this."
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk \\
        --java-options -Xmx${avail_mem}g \\
        GenotypeGVCFs \\
        $options.args \\
        $intervalsOption \\
        $dbsnpOption \\
        -R $fasta \\
        -V $gvcf \\
        -O ${prefix}.genotyped.vcf.gz

    gatk GenotypeGVCFs --version 2> /dev/null | grep "GATK" | sed -e 's/The Genome Analysis Toolkit (GATK) v//' > ${software}.version.txt
    """
}
