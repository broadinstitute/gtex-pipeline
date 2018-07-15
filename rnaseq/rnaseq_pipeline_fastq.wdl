import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:star_v1-0_BETA/versions/3/plain-WDL/descriptor" as star_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:markduplicates_v1-0_BETA/versions/1/plain-WDL/descriptor" as markduplicates_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rsem_v1-0_BETA/versions/3/plain-WDL/descriptor" as rsem_wdl
import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rnaseqc_counts_v1-0_BETA/versions/2/plain-WDL/descriptor" as rnaseqc_wdl

workflow rnaseq_pipeline_fastq_workflow {

    File fastq1
    File? fastq2
    String prefix

    call star_wdl.star {
        input: fastq1=fastq1, fastq2=fastq2, prefix=prefix
    }

    call markduplicates_wdl.markduplicates {
        input: input_bam=star.bam_file, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rnaseqc_wdl.rnaseqc_counts {
        input: bam_file=markduplicates.bam_file, bam_index=markduplicates.bam_index, prefix=prefix
    }
}
