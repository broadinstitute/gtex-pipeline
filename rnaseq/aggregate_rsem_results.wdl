task aggregate_rsem_results {

    Array[File] rsem_isoform_results
    Array[File] rsem_gene_results
    String prefix

    Int memory
    Int disk_space

    command {
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_isoform_results} TPM "${prefix}_rsem_isoforms_tpm"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_gene_results} TPM "${prefix}_rsem_genes_tpm"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_isoform_results} expected_count "${prefix}_rsem_isoforms_counts"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_gene_results} expected_count "${prefix}_rsem_genes_counts"
    }

    output {
        File RSEM_isoforms_TPM="${prefix}_rsem_isoforms_tpm.txt.gz"
        File RSEM_genes_TPM="${prefix}_rsem_genes_tpm.txt.gz"
        File RSEM_isoforms_counts="${prefix}_rsem_isoforms_counts.txt.gz"
        File RSEM_genes_counts="${prefix}_rsem_genes_counts.txt.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow aggregate_rsem_results_workflow {
    call aggregate_rsem_results
}
