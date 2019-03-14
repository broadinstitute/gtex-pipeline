task leafcutter_cluster {

    Array[File] junc_files
    File exon_list
    File genes_gtf
    String prefix

    Int? min_clu_reads
    Float? min_clu_ratio
    Int? max_intron_len
    Int? num_pcs

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        python /src/cluster_prepare_fastqtl.py \
            ${write_lines(junc_files)} \
            ${exon_list} \
            ${genes_gtf} \
            ${prefix} \
            ${"--min_clu_reads " + min_clu_reads} \
            ${"--min_clu_ratio " + min_clu_ratio} \
            ${"--max_intron_len " + max_intron_len} \
            ${"--num_pcs " + num_pcs} 
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File counts="${prefix}_perind.counts.gz"
        File counts_numers="${prefix}_perind_numers.counts.gz"
        File clusters_pooled="${prefix}_pooled.gz"
        File clusters_refined="${prefix}_refined.gz"
        File clusters_to_genes="${prefix}.leafcutter.clusters_to_genes.txt"
        File phenotype_groups="${prefix}.leafcutter.phenotype_groups.txt"
        File leafcutter_bed="${prefix}.leafcutter.bed.gz"
        File leafcutter_bed_index="${prefix}.leafcutter.bed.gz.tbi"
        File leafcutter_pcs="${prefix}.leafcutter.PCs.txt"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_cluster_workflow {
    call leafcutter_cluster
}
