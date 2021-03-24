#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_GENOTYPEGVCFS } from '../../../../software/gatk4/genotypegvcfs/main.nf' addParams( options: [:] )

// Basic parameters with uncompressed VCF input test
workflow test_gatk4_genotypegvcfs_vcf_input {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.idx", checkIfExists: true) ]

    fasta        = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta", checkIfExists: true)
    fastaIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta.fai", checkIfExists: true)
    fastaDict    = file("${launchDir}/tests/data/genomics/hg38/genome/genome.dict", checkIfExists: true)

    GATK4_GENOTYPEGVCFS ( input, fasta, fastaIndex, fastaDict, [], [], [] )
}

// Basic parameters test (with compressed VCF input)
workflow test_gatk4_genotypegvcfs_gz_input {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.gz", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.gz.tbi", checkIfExists: true) ]

    fasta        = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta", checkIfExists: true)
    fastaIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta.fai", checkIfExists: true)
    fastaDict    = file("${launchDir}/tests/data/genomics/hg38/genome/genome.dict", checkIfExists: true)

    GATK4_GENOTYPEGVCFS ( input, fasta, fastaIndex, fastaDict, [], [], [] )
}

// Basic parameters + optional dbSNP test
workflow test_gatk4_genotypegvcfs_vcf_input_dbsnp {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.idx", checkIfExists: true) ]

    fasta        = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta", checkIfExists: true)
    fastaIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta.fai", checkIfExists: true)
    fastaDict    = file("${launchDir}/tests/data/genomics/hg38/genome/genome.dict", checkIfExists: true)

    dbsnp        = file("${launchDir}/tests/data/genomics/hg38/genome/dbsnp.vcf.gz", checkIfExists: true)
    dbsnpIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/dbsnp.vcf.gz.tbi", checkIfExists: true)

    GATK4_GENOTYPEGVCFS ( input, fasta, fastaIndex, fastaDict, dbsnp, dbsnpIndex, [] )
}

// Basic parameters + optional intervals test
workflow test_gatk4_genotypegvcfs_vcf_input_intervals {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.idx", checkIfExists: true) ]

    fasta        = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta", checkIfExists: true)
    fastaIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta.fai", checkIfExists: true)
    fastaDict    = file("${launchDir}/tests/data/genomics/hg38/genome/genome.dict", checkIfExists: true)

    intervalsBed = file("${launchDir}/tests/data/genomics/hg38/genome/bed/intervals1.bed", checkIfExists: true)

    GATK4_GENOTYPEGVCFS ( input, fasta, fastaIndex, fastaDict, [], [], intervalsBed )
}

// Basic parameters + optional dbSNP + optional intervals test
workflow test_gatk4_genotypegvcfs_vcf_input_dbsnp_intervals {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf", checkIfExists: true),
              file("${launchDir}/tests/data/genomics/hg38/illumina/vcf/test1.g.vcf.idx", checkIfExists: true) ]

    fasta        = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta", checkIfExists: true)
    fastaIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/genome.fasta.fai", checkIfExists: true)
    fastaDict    = file("${launchDir}/tests/data/genomics/hg38/genome/genome.dict", checkIfExists: true)

    dbsnp        = file("${launchDir}/tests/data/genomics/hg38/genome/dbsnp.vcf.gz", checkIfExists: true)
    dbsnpIndex   = file("${launchDir}/tests/data/genomics/hg38/genome/dbsnp.vcf.gz.tbi", checkIfExists: true)

    intervalsBed = file("${launchDir}/tests/data/genomics/hg38/genome/bed/intervals1.bed", checkIfExists: true)

    GATK4_GENOTYPEGVCFS ( input, fasta, fastaIndex, fastaDict, dbsnp, dbsnpIndex, intervalsBed )
}
