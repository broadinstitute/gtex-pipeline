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

	call fp.CrossCheckSample{
		input:
		first=star.bam_file,
		first_index=star.bam_file_index,
		second=get_het_vcfs.snps_vcf,
		second_index=get_het_vcfs.snps_vcf_index
	}
}