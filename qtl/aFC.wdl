task aFC {

    File vcf_file
    File vcf_index
    File expression_bed
    File expression_bed_index
    File covariates_file
    File afc_qtl_file
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        python3 /opt/aFC/aFC.py \
            --vcf ${vcf_file} \
            --pheno ${expression_bed} \
            --qtl ${afc_qtl_file} \
            --cov ${covariates_file} \
            --log_xform 1 \
            --log_base 2 \
            --o ${prefix}.aFC.txt
        gzip ${prefix}.aFC.txt
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File afc_file="${prefix}.aFC.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow aFC_workflow {
    call aFC
}
