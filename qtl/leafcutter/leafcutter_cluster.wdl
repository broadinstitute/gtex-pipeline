version 1.0

task leafcutter_cluster {
    input {
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
    }

    command <<<
        set -euo pipefail
        pip3 install qtl # TODO: add this to the docker image

        R -e 'install.packages(c("dplyr","foreach"))' #TODO: add this to docker image

        echo << EOF > temp_map.tsv

        one	one
        two	two
        three	three
        EOF

        python3 /src/cluster_prepare_fastqtl.py \
            "~{write_lines(junc_files)}" \
            "~{exon_list}" \
            "~{genes_gtf}" \
            "~{prefix}" \
            "temp_map.tsv" \
            ~{"--min_clu_reads " + min_clu_reads} \
            ~{"--min_clu_ratio " + min_clu_ratio} \
            ~{"--max_intron_len " + max_intron_len} \
            ~{"--num_pcs " + num_pcs} 
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_threads
        preemptible: num_preempt
    }

    output {
        File counts="~{prefix}_per_ind.counts.gz"
        File counts_numbers="~{prefix}_per_ind_numbers.counts.gz"
        File clusters_pooled="~{prefix}_pooled.gz"
        File clusters_refined="~{prefix}_refined.gz"
        File clusters_to_genes="~{prefix}.leafcutter.clusters_to_genes.txt"
        File phenotype_groups="~{prefix}.leafcutter.phenotype_groups.txt"
        File leafcutter_bed="~{prefix}.leafcutter.bed.gz"
        File leafcutter_bed_index="~{prefix}.leafcutter.bed.gz.tbi"
        File leafcutter_pcs="~{prefix}.leafcutter.PCs.txt"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_cluster_workflow {
    call leafcutter_cluster

    output {
        File leafcutter_counts=leafcutter_cluster.counts
        File leafcutter_counts_numbers=leafcutter_cluster.counts_numbers
        File leafcutter_clusters_pooled=leafcutter_cluster.clusters_pooled
        File leafcutter_clusters_refined=leafcutter_cluster.clusters_refined
        File leafcutter_clusters_to_genes=leafcutter_cluster.clusters_to_genes
        File leafcutter_phenotype_groups=leafcutter_cluster.phenotype_groups
        File leafcutter_bed=leafcutter_cluster.leafcutter_bed
        File leafcutter_bed_index=leafcutter_cluster.leafcutter_bed_index
        File leafcutter_pcs=leafcutter_cluster.leafcutter_pcs
    }
}
