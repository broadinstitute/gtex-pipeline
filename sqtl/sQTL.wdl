version 1.0

import "../genotype/participant_vcfs.wdl" as participant_vcfs
import "../rnaseq/star.wdl" as rnaseq
import "../utils/Fingerprinting.wdl" as fp

workflow sQTLAnalysis{
	input {

	}	

	call participant_vcfs.participant_vcfs as get_het_vcfs{} 

	call rnaseq.star as star{
		input: 
		varVCFfile=get_het_vcfs.snps_vcf
	}

	call fp.CrossCheckSample as fingerprint{
		input:
		first=star.bam_file,
		first_index=star.bam_index,
		second=get_het_vcfs.snps_vcf,
		second_index=get_het_vcfs.snps_vcf_index
	}


	output {
		File bam_file = star.bam_file 
		File bam_index = star.bam_index 
		File transcriptome_bam = star.transcriptome_bam 
		File chimeric_junctions = star.chimeric_junctions 
		File chimeric_bam_file = star.chimeric_bam_file 
		File chimeric_bam_index = star.chimeric_bam_index 
		File read_counts = star.read_counts 
		File junctions = star.junctions 
		File junctions_pass1 = star.junctions_pass1 
		Array[File] star_logs = star.logs

		File fingerprint_metrics=fingerprint.metrics
	}
}