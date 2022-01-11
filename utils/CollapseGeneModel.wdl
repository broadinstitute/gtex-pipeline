version 1.0


task CollapseGeneModel {
	input {
		File gene_annotation_file
		String output_prefix
		File collapse_annotation_script
	}
	command<<<
		set -euo pipefail

		pip3 install bx-python
		python3 ~{collapse_annotation_script} ~{gene_annotation_file} ~{output_prefix}.genes.gtf
	>>>
	runtime {
		docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
		memory: "4GB"
		disks: "local-disk 10 HDD"
		cpu: 1
		preemptible: 1
	}

	output {
		File genes = output_prefix +".genes.gtf"
	}

	meta {
		author: "Yossi Farjoun"
	}
}

task ExtractExonList {
	input {
		File collapsed_annotation_file
		String output_prefix
	}
	command <<<
	pip3 install qtl
	python3 <<EOF

	import pandas as pd
	import qtl.annotation

	annot = qtl.annotation.Annotation('~{collapsed_annotation_file}')
	exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
							for g in annot.genes for e in g.transcripts[0].exons],
								columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
	exon_df.to_csv('~{output_prefix}.exons.txt', sep='\t', index=False)
		
	EOF
	>>>
	output {
		File exons = output_prefix+".exons.txt"
	}

	runtime {
		docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
		memory: "5 GB"
		disks: "local-disk 50 HDD"
		cpu: 1
		preemptible: 1
	}

	meta {
		author: "Yossi Farjoun"
	}
}


workflow CollapseAnnotationFile{
	input {
		File gene_annotation_file
		String output_prefix
		File collapse_annotation_script
	}

	call CollapseGeneModel {
		input:
			gene_annotation_file = gene_annotation_file,
			output_prefix = output_prefix,
			collapse_annotation_script = collapse_annotation_script
	}

	call ExtractExonList {
		input:
			collapsed_annotation_file=CollapseGeneModel.genes,
			output_prefix=output_prefix
	}

	output {
		File genes=CollapseGeneModel.genes
		File exons=ExtractExonList.exons
	}
}
