task combine_gcts {

    Array[File] rsem_isoform_results
    Array[File] rsem_gene_results
    String prefix

    Int memory
    Int disk_space

    command {
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_isoform_results} TPM "${prefix}_RSEM_isoforms_TPM"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_gene_results} TPM "${prefix}_RSEM_genes_TPM"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_isoform_results} expected_count "${prefix}_RSEM_isoforms_counts"
        python3 /src/aggregate_rsem_results.py ${sep=' ' rsem_gene_results} expected_count "${prefix}_RSEM_genes_counts"
    }

    output {
        File RSEM_isoforms_TPM="${prefix}_RSEM_isoforms_TPM.txt.gz"
        File RSEM_genes_TPM="${prefix}_RSEM_genes_TPM.txt.gz"
        File RSEM_isoforms_counts="${prefix}_RSEM_isoforms_counts.txt.gz"
        File RSEM_genes_counts="${prefix}_RSEM_genes_counts.txt.gz"
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


workflow combine_gcts_workflow {
    call combine_gcts
}
