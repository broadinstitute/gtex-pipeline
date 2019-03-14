    task ase_gatk_readcounter {

    File gatk_jar
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dict
    File het_vcf
    File het_vcf_index
    File bam_file
    File bam_index
    String prefix
    Boolean? filter_wasp = false

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command <<<
        set -euo pipefail
        if [[ ${filter_wasp} = "true" ]]
        then
            echo $(date +"[%b %d %H:%M:%S] Filtering out reads with allelic mapping bias")
            samtools view -h ${bam_file} | grep -v "vW:i:[2-7]" | samtools view -1 > filtered.bam
            samtools index filtered.bam
            python3 /src/run_GATK_ASEReadCounter.py ${gatk_jar} ${genome_fasta} ${het_vcf} filtered.bam ${prefix}
        else
            python3 /src/run_GATK_ASEReadCounter.py ${gatk_jar} ${genome_fasta} ${het_vcf} ${bam_file} ${prefix}
        fi

        # filter out chrX
        mv ${prefix}.readcounts.txt.gz ${prefix}.readcounts.all.txt.gz
        zcat ${prefix}.readcounts.all.txt.gz | awk '$1!="chrX" && $1!="X" {print $0}' | gzip -c > ${prefix}.readcounts.txt.gz
        zcat ${prefix}.readcounts.all.txt.gz | awk '$1=="contig" || $1=="chrX" || $1=="X" {print $0}' | gzip -c > ${prefix}.readcounts.chrX.txt.gz
    >>>

    output {
        File ase_read_counts = "${prefix}.readcounts.txt.gz"
        File ase_read_counts_chrX = "${prefix}.readcounts.chrX.txt.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow ase_gatk_readcounter_workflow {
    call ase_gatk_readcounter
}
