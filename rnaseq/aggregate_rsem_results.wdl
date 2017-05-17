task aggregate_rsem_results {

    Array[File] rsem_isoform_outputs
    Array[File] rsem_gene_outputs
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        # workaround for broken 'write_lines'
        echo $(date +"[%b %d %H:%M:%S] Writing inputs to file")
        python3 <<CODE
with open("rsem_isoform_paths.tsv", "w") as f:
    for p in '${sep="," rsem_isoform_outputs}'.split(","):
        f.write(p+'\n')
with open("rsem_gene_paths.tsv", "w") as f:
    for p in '${sep="," rsem_gene_outputs}'.split(","):
        f.write(p+'\n')
CODE
        echo $(date +"[%b %d %H:%M:%S] Combining transcript-level output")
        python3 /src/aggregate_rsem_results.py rsem_isoform_paths.tsv TPM IsoPct expected_count ${prefix}
        echo $(date +"[%b %d %H:%M:%S] Combining gene-level output")
        python3 /src/aggregate_rsem_results.py rsem_gene_paths.tsv TPM expected_count ${prefix}
    }

    output {
        File rsem_transcripts_tpm="${prefix}.rsem_transcripts_tpm.txt.gz"
        File rsem_transcripts_isopct="${prefix}.rsem_transcripts_isopct.txt.gz"
        File rsem_transcripts_expected_count="${prefix}.rsem_transcripts_expected_count.txt.gz"
        File rsem_genes_tpm="${prefix}.rsem_genes_tpm.txt.gz"
        File rsem_genes_expected_count="${prefix}.rsem_genes_expected_count.txt.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow aggregate_rsem_results_workflow {
    call aggregate_rsem_results
}
