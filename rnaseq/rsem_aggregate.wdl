task rsem_aggregate_results {

    Array[File] rsem_isoforms
    Array[File] rsem_genes
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        echo $(date +"[%b %d %H:%M:%S] Combining transcript-level output")
        python3 /src/aggregate_rsem_results.py ${write_lines(rsem_isoforms)} TPM IsoPct expected_count ${prefix}
        echo $(date +"[%b %d %H:%M:%S] Combining gene-level output")
        python3 /src/aggregate_rsem_results.py ${write_lines(rsem_genes)} TPM expected_count ${prefix}
    }

    output {
        File transcripts_tpm="${prefix}.rsem_transcripts_tpm.txt.gz"
        File transcripts_isopct="${prefix}.rsem_transcripts_isopct.txt.gz"
        File transcripts_expected_count="${prefix}.rsem_transcripts_expected_count.txt.gz"
        File genes_tpm="${prefix}.rsem_genes_tpm.txt.gz"
        File genes_expected_count="${prefix}.rsem_genes_expected_count.txt.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rsem_aggregate_results_workflow {
    call rsem_aggregate_results
}
