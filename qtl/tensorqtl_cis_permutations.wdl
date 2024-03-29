task tensorqtl_cis_permutations {

    File plink_pgen
    File plink_pvar
    File plink_psam

    File phenotype_bed
    File covariates
    String prefix

    File? phenotype_groups
    Float? fdr
    Float? qvalue_lambda
    Int? seed
    String? flags

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        plink_base=$(echo "${plink_pgen}" | rev | cut -f 2- -d '.' | rev)
        python3 -m tensorqtl \
            $plink_base ${phenotype_bed} ${prefix} \
            --mode cis \
            --covariates ${covariates} \
            ${"--phenotype_groups " + phenotype_groups} \
            ${"--fdr " + fdr} \
            ${"--qvalue_lambda " + qvalue_lambda} \
            ${"--seed " + seed} \
            ${flags}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-central1-c"]
    }

    output {
        File cis_qtl="${prefix}.cis_qtl.txt.gz"
        File log="${prefix}.tensorQTL.cis.log"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_cis_permutations_workflow {
    call tensorqtl_cis_permutations
}
