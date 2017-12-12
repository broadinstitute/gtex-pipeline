task eqtl_prepare_expression {

    File tpm_gct
    File counts_gct
    File annotation_gtf
    File sample_participant_ids
    File vcf_chr_list
    String prefix

    Float? tpm_threshold
    Int? count_threshold
    Float? sample_frac_threshold
    String? normalization_method
    String? flags  # --convert_tpm, --legacy_mode

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} \
        ${annotation_gtf} ${sample_participant_ids} ${vcf_chr_list} ${prefix} \
        ${"--tpm_threshold " + tpm_threshold} \
        ${"--count_threshold " + count_threshold} \
        ${"--sample_frac_threshold " + sample_frac_threshold} \
        ${"--normalization_method " + normalization_method} \
        ${flags}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File expression_bed="${prefix}.expression.bed.gz"
        File expression_bed_index="${prefix}.expression.bed.gz.tbi"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow eqtl_prepare_expression_workflow {
    call eqtl_prepare_expression
}
