task star {

    File fastq1
    File? fastq2
    String prefix
    File star_index

    # STAR options
    Int? outFilterMultimapNmax
    Int? alignSJoverhangMin
    Int? alignSJDBoverhangMin
    Int? outFilterMismatchNmax
    Float? outFilterMismatchNoverLmax
    Int? alignIntronMin
    Int? alignIntronMax
    Int? alignMatesGapMax
    String? outFilterType
    Float? outFilterScoreMinOverLread
    Float? outFilterMatchNminOverLread
    Int? limitSjdbInsertNsj
    String? outSAMstrandField
    String? outFilterIntronMotifs
    String? alignSoftClipAtReferenceEnds
    String? quantMode
    String? outSAMattrRGline
    String? outSAMattributes
    Int? chimSegmentMin
    Int? chimJunctionOverhangMin
    String? chimOutType
    Int? chimMainSegmentMultNmax
    File? sjdbFileChrStartEnd

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir star_out
        /src/run_STAR.py \
            star_index $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir star_out \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            --threads ${num_threads}
    }

    output {
        File bam_file = "star_out/${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/${prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "star_out/${prefix}.Chimeric.out.junction"
        File chimeric_bam_file = "star_out/${prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "star_out/${prefix}.ReadsPerGene.out.tab"
        File junctions = "star_out/${prefix}.SJ.out.tab"
        File junctions_pass1 = "star_out/${prefix}._STARpass1/SJ.out.tab"
        Array[File] logs = ["star_out/${prefix}.Log.final.out", "star_out/${prefix}.Log.out", "star_out/${prefix}.Log.progress.out"]
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


task rsem {

    File transcriptome_bam
    File rsem_reference
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    Int? max_frag_len
    String? estimate_rspd
    String? is_stranded
    String? paired_end

    command {
        set -euo pipefail
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            ${"--paired_end " + paired_end} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
    }

    output {
        File genes="${prefix}.rsem.genes.results"
        File isoforms="${prefix}.rsem.isoforms.results"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


task rnaseqc {
    File bam_file
    File bam_index
    File genes_gtf
    File genome_fasta
    File genome_fasta_index

    String prefix
    String? gatk_flags

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${bam_index}
        touch ${genome_fasta_index}
        python3 /src/run_rnaseqc.py ${bam_file} ${genes_gtf} ${genome_fasta} ${prefix}\
            --java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java\
            --memory ${memory}\
            --rnaseqc_flags noDoC strictMode\
            ${" --gatk_flags " + gatk_flags}
    }

    output {
        File gene_rpkm = "${prefix}.gene_rpkm.gct.gz"
        File gene_counts = "${prefix}.gene_reads.gct.gz"
        File exon_counts = "${prefix}.exon_reads.gct.gz"
        File count_metrics = "${prefix}.metrics.tsv"
        File count_outputs = "${prefix}.tar.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseq_pipeline_fastq_workflow {

    File fastq1
    File fastq2
    String prefix

    call star {
        input: fastq1=fastq1, fastq2=fastq2, prefix=prefix
    }

    call rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rnaseqc {
        input: bam_file=star.bam_file, bam_index=star.bam_index, prefix=prefix
    }
}
