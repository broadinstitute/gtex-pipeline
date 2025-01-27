task dapars {

    Array[File] bigwig_files
    Array[String] sample_ids
    String prefix
    File utr_annotation
    File size_factors
    File? sample_participant_lookup
    File? expression_bed
    Int? coverage_threshold = 10

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        python3 <<CODE
with open("bigwig_list.txt", "w") as f:
    f.write('sample_id\tbigwig\n')
    for sample_id, bigwig in zip('${sep="," sample_ids}'.split(","), '${sep="," bigwig_files}'.split(",")):
        f.write('\t'.join([sample_id, bigwig])+'\n')
CODE

    python3 /src/dapars.py \
        --bigwig_files bigwig_list.txt \
        --utr ${utr_annotation} \
        --size_factors ${size_factors} \
        --prefix ${prefix} \
        ${"--sample_participant_lookup " + sample_participant_lookup} \
        ${"--expression_bed " + expression_bed} \
        --coverage_threshold ${coverage_threshold} \
        --threads ${num_threads} \
        --parquet
    }

    output {
        File phenotype_groups = "${prefix}.dapars2_phenotype_groups.txt"
        File dapars_bed = "${prefix}.dapars2_phenotypes.bed.parquet"
        File dapars_ratios = "${prefix}.dapars2_3p_UTR_ratios.parquet"
        File bigwig_list = "bigwig_list.txt"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow dapars_workflow {
    call dapars
}
