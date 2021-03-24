#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_INDEXFEATUREFILE } from '../../../../software/gatk4/indexfeaturefile/main.nf' addParams( options: [:] )

workflow test_gatk4_indexfeaturefile_bed {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed", checkIfExists: true) ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_bed_gz {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/genome/bed/test.bed.gz", checkIfExists: true) ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_vcf {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf", checkIfExists: true) ]

    GATK4_INDEXFEATUREFILE ( input )
}

workflow test_gatk4_indexfeaturefile_vcf_gz {
    
    def input = []
    input = [ [ id:'test' ], // meta map
              file("${launchDir}/tests/data/genomics/sarscov2/illumina/vcf/test.vcf.gz", checkIfExists: true) ]

    GATK4_INDEXFEATUREFILE ( input )
}
