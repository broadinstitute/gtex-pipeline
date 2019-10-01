task tensorqtl_cis_nominal {

    File plink_bed
    File plink_bim
    File plink_fam

    File phenotype_bed
    File covariates
    String prefix

    File? interaction
    File? phenotype_groups

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        plink_base=$(echo "${plink_bed}" | rev | cut -f 2- -d '.' | rev)
        python3 -m tensorqtl \
            $plink_base ${phenotype_bed} ${prefix} \
            --mode cis_nominal \
            --covariates ${covariates} \
            ${"--interaction " + interaction} \
            ${"--phenotype_groups " + phenotype_groups}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: ["us-central1-c"]
    }

    output {
        Array[File] chr_parquet=glob("${prefix}*.parquet")
        File log=glob("${prefix}*.log")[0]
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_cis_nominal_workflow {
    call tensorqtl_cis_nominal
}
