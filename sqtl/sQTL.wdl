version 1.0

import "genotype/participant_vcfs.wdl" as participant_vcfs
import "rnaseq/star.wdl" as rnaseq

workflow sQTLAnalysis{
	input {

	}	

	call participant_vcfs.participant_vcfs as get_het_vcfs{} 

	call rnaseq.star{
		input: 
		varVCFfile=get_het_vcfs.snps_vcf
	}


}