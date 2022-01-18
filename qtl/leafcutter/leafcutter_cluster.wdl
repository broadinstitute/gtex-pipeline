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
        set -exuo pipefail
        pip3 install qtl # TODO: add this to the docker image

        R -e 'install.packages(c("dplyr","foreach"))' #TODO: add this to docker image

        # the list of files -> take their basename, and replicate the part prior to the regtools ending
        # with a \t separating, and put the results into temp_map.tsv
        cat ~{write_lines(junc)} | xargs -n1 basename | sed 's/\(.*\).regtools_junc.txt.gz/\1\t&/' > temp_map.tsv

        touch "~{prefix}_perind.counts.gz"
        touch "~{prefix}_perind_numbers.counts.gz"
        touch "~{prefix}_perind.counts.filtered.gz"
        touch "~{prefix}_pooled.gz"
        touch "~{prefix}_refined.gz"
        touch "~{prefix}.leafcutter.clusters_to_genes.txt"
        touch "~{prefix}.leafcutter.phenotype_groups.txt"
        touch "~{prefix}.leafcutter.bed.gz"
        touch "~{prefix}.leafcutter.bed.gz.tbi"
        touch "~{prefix}.leafcutter.PCs.txt"

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
        File counts="~{prefix}_perind.counts.gz"
        File counts_numbers="~{prefix}_perind_numbers.counts.gz"
        File counts_numbers_filtered="~{prefix}_perind.counts.filtered.gz"
        File clusters_pooled="~{prefix}_pooled.gz"
        File clusters_refined="~{prefix}_refined.gz"
        File clusters_to_genes="~{prefix}.leafcutter.clusters_to_genes.txt"
        File phenotype_groups="~{prefix}.leafcutter.phenotype_groups.txt"
        File leafcutter_bed="~{prefix}.leafcutter.bed.gz"
        File leafcutter_bed_index="~{prefix}.leafcutter.bed.gz.tbi"
        File leafcutter_pcs="~{prefix}.leafcutter.PCs.txt"
        File map="temp_map.tsv"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_cluster_workflow {
    input {
        Array[File] junc_files
    }

    call leafcutter_cluster{
    input:
        disk_space=20+ceil(size(junc_files, "GB")),
        junc_files=junc_files
    }

    output {
        File leafcutter_counts=leafcutter_cluster.counts
        File leafcutter_counts_numbers=leafcutter_cluster.counts_numbers
        File counts_numbers_filtered=leafcutter_cluster.counts_numbers_filtered
        File leafcutter_clusters_pooled=leafcutter_cluster.clusters_pooled
        File leafcutter_clusters_refined=leafcutter_cluster.clusters_refined
        File leafcutter_clusters_to_genes=leafcutter_cluster.clusters_to_genes
        File leafcutter_phenotype_groups=leafcutter_cluster.phenotype_groups
        File leafcutter_bed=leafcutter_cluster.leafcutter_bed
        File leafcutter_bed_index=leafcutter_cluster.leafcutter_bed_index
        File leafcutter_pcs=leafcutter_cluster.leafcutter_pcs
    }
}
