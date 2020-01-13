task phaser {

    String individual_id
    String chromosome
    Array[File]+ bam_files
    Array[File]+ bam_indices
    File genotype_vcf
    File genotype_index
    File gene_model_bed

    File? haplotype_blacklist_bed
    File? phaser_blacklist_bed

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        export TMP=$HOME/
        set -euo pipefail

        echo $(date +"[%b %d %H:%M:%S] Preparing indices")
        for index in ${sep=" " bam_indices}; do
            touch $index
        done

        echo $(date +"[%b %d %H:%M:%S] Preparing bam files")
        mkdir ./bam_staging
        for bam_file in ${sep=" " bam_files}; do
            samtools view -h $bam_file ${chromosome} | \
            grep -v "vW:i:[2-7]" | \
            samtools view -h1 | samtools sort > ./bam_staging/$(basename $bam_file)
            samtools index -@ ${num_threads} ./bam_staging/$(basename $bam_file)
        done

        touch ${genotype_index}
        echo $(date +"[%b %d %H:%M:%S] Running phASER")
        python2.7 /opt/phaser/wrapper.py phase ${individual_id} ./bam_staging/*.bam \
            ${genotype_vcf} ${gene_model_bed} .  \
            ${"--haplo-count-blacklist=" + haplotype_blacklist_bed} \
            ${"--blacklist=" + phaser_blacklist_bed} \
            --chr ${chromosome}
    }

    output {
        File allele_config = "${individual_id}.${chromosome}.allele_config.txt"
        File allelic_counts = "${individual_id}.${chromosome}.allelic_counts.txt"
        File haplotypes = "${individual_id}.${chromosome}.haplotypes.txt"
        File haplo_counts = "${individual_id}.${chromosome}.haplotypic_counts.txt"
        File variant_connections = "${individual_id}.${chromosome}.variant_connections.txt"
        File phase_vcf = "${individual_id}.${chromosome}.vcf.gz"
        File vcf_index = "${individual_id}.${chromosome}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/just-episode-184015/rnaseqc:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Aaron Graubert"
    }
}

task phaser_postprocess {

    String individual_id
    Array[File]+ chromosome_vcfs
    Array[File]+ chromosome_haplotype_counts
    File gene_model_bed

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail

        echo $(date +"[%b %d %H:%M:%S] Processing phASER Files")
        python2.7 /opt/phaser/wrapper.py postprocess ${individual_id} \
          ${write_lines(chromosome_vcfs)} \
          ${write_lines(chromosome_haplotype_counts)} \
          ${gene_model_bed} \
          .
    }

    output {
        File haplo_counts = "${individual_id}.haplotypic_counts.txt.gz"
        File phase_vcf = "${individual_id}.vcf.gz"
        File vcf_index = "${individual_id}.vcf.gz.tbi"
        File expression_counts = "${individual_id}.gene_ae.txt.gz"
    }

    runtime {
        docker: "us.gcr.io/just-episode-184015/rnaseqc:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Aaron Graubert"
    }
}


workflow phaser_workflow {
    String individual_id
    File gene_model_bed
    File contig_list

    scatter(chr in read_lines(contig_list))
    {
      call phaser { input: individual_id=individual_id, gene_model_bed=gene_model_bed, chromosome=chr}
    }

    call phaser_postprocess {
      input: individual_id=individual_id, gene_model_bed=gene_model_bed, chromosome_vcfs=phaser.phase_vcf, chromosome_haplotype_counts=phaser.haplo_counts
    }

    output {
        File haplo_counts = phaser_postprocess.haplo_counts
        File phase_vcf = phaser_postprocess.phase_vcf
        File vcf_index = phaser_postprocess.vcf_index
        File expression_counts = phaser_postprocess.expression_counts
    }
}
