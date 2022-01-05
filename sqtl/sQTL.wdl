version 1.0

import "../genotype/participant_vcfs.wdl" as participant_vcfs
import "../rnaseq/star.wdl" as rnaseq
import "../utils/Fingerprinting.wdl" as fp
import "../utils/IdentifySample.wdl" as id

workflow sQTLAnalysis{
	input {
		File hapMap
		Boolean crosscheck=true
	}	

	call participant_vcfs.participant_vcfs as get_het_vcfs{} 

	call rnaseq.star as star{
		input: 
		varVCFfile=get_het_vcfs.snps_vcf
	}
	if (crosscheck) {
		call fp.CrossCheckSample as fingerprint {
			input:
			first=star.bam_file,
			first_index=star.bam_index,
			second=get_het_vcfs.snps_vcf,
			second_index=get_het_vcfs.snps_vcf_index,
			hapMap = hapMap	
		}

		call id.IdentifySampleWF as identifySample{
			input:
			sample=star.bam_file,
	        sample_index=star.bam_index,
	        hapMap=hapMap
		}

	}
	output {
		File snps_vcf = get_het_vcfs.snps_vcf
		File snps_vcf_index = get_het_vcfs.snps_vcf_index
		File snps_het_vcf = get_het_vcfs.snps_het_vcf
		File snps_het_vcf_index = get_het_vcfs.snps_het_vcf_index

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

		File? fingerprint_metrics=fingerprint.metrics

		File? clustered_metrics=identifySample.fp_clustered
		File? fingerprint_matrix=identifySample.fp_metrics
		String? fingerprint_match=identifySample.match_group

	}
}