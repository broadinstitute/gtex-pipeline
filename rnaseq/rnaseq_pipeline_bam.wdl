task samtofastq {

    File input_bam
    String prefix

    Int memory
    Int disk_space
    Int num_preempt

    command {
        # make sure path is absolute
        input_bam_abs=${input_bam}
        if [[ $input_bam_abs != /* ]]; then
            input_bam_abs=$PWD/$input_bam_abs
        fi

        mkdir samtofastq  # workaround for named pipes
        python3 -u /src/run_SamToFastq.py $input_bam_abs -p ${prefix} --output_dir samtofastq
        mv samtofastq/${prefix}_*.fastq.gz .
    }

    output {
        File fastq1="${prefix}_1.fastq.gz"
        File fastq2="${prefix}_2.fastq.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


task star {
    
    File fastq1
    File fastq2
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
    
    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
    
    command {
        set -euo pipefail
        
        # make sure paths are absolute
        fastq1_abs=${fastq1}
        fastq2_abs=${fastq2}
        if [[ $fastq1_abs != /* ]]; then
            fastq1_abs=$PWD/$fastq1_abs
            fastq2_abs=$PWD/$fastq2_abs
        fi
        
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
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
    
    meta {
        author: "Francois Aguet"
    }
}


task bamsync {

    File source_bam
    File target_bam
    File target_bam_index
    String prefix

    Int memory
    Int disk_space
    Int num_preempt

    command {
        echo $(date +"[%b %d %H:%M:%S] Running bamsync")
        /src/run_bamsync.sh ${source_bam} ${target_bam} ${prefix}
        echo $(date +"[%b %d %H:%M:%S] Running samtools flagstat")
        samtools flagstat ${prefix}.Aligned.sortedByCoord.out.patched.bam > ${prefix}.Aligned.sortedByCoord.out.patched.bam.flagstat.txt
    }

    output {
        File patched_bam_file="${prefix}.Aligned.sortedByCoord.out.patched.bam"
        File patched_bam_index="${prefix}.Aligned.sortedByCoord.out.patched.bam.bai"
        File patched_bam_flagstat="${prefix}.Aligned.sortedByCoord.out.patched.bam.flagstat.txt"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
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

    command {
        mkdir rsem_reference
        tar -xvvf ${rsem_reference} -C rsem_reference --strip-components=1

        /src/run_RSEM.py \
            ${"--max_frag_len " + max_frag_len} \
            ${"--estimate_rspd " + estimate_rspd} \
            ${"--is_stranded " + is_stranded} \
            --threads ${num_threads} \
            rsem_reference ${transcriptome_bam} ${prefix}
    }

    output {
        File genes="${prefix}.rsem.genes.results"
        File isoforms="${prefix}.rsem.isoforms.results"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
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
    String notes
    String? gatk_flags

    Int memory
    Int disk_space
    Int num_preempt

    command {
        set -e
        touch ${bam_index}
        touch ${genome_fasta_index}

        /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx${memory}g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
        -s ${prefix},${bam_file},${notes} -t ${genes_gtf} -r ${genome_fasta} -o . -noDoC -strictMode ${" -gatkFlags " + gatk_flags}
        
        # remove tmp files
        rm genes.rpkm.gct
        rm ${prefix}/${prefix}.metrics.tmp.txt
        rm ${prefix}/${prefix}.metrics.txt
        
        mv ${prefix}/${prefix}.transcripts.rpkm.gct ${prefix}.gene_rpkm.gct
        python3 /src/convert_rnaseqc_counts.py ${prefix}.gene_rpkm.gct ${prefix}/${prefix}.exon_intron_report.txt --exon_report ${prefix}/${prefix}.exon_report.txt ${genes_gtf} -o .
        mv metrics.tsv ${prefix}.metrics.tsv
        gzip ${prefix}.gene_rpkm.gct
        tar -cvzf ${prefix}.tar.gz ${prefix}/*
    }

    output {
        File gene_rpkm = "${prefix}.gene_rpkm.gct.gz"
        File gene_counts = "${prefix}.gene_reads.gct.gz"
        File exon_counts = "${prefix}.exon_reads.gct.gz"
        File count_metrics = "${prefix}.metrics.tsv"
        File count_outputs = "${prefix}.tar.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseq_pipeline_bam_workflow {

    File input_bam
    String prefix

    File star_index
    File rsem_reference

    File genes_gtf
    File genome_fasta
    File genome_fasta_index
    String rnaseqc_notes

    Int samtofastq_memory
    Int samtofastq_disk
    Int samtofastq_preempt

    Int star_memory
    Int star_disk
    Int star_threads
    Int star_preempt

    Int bamsync_memory
    Int bamsync_disk
    Int bamsync_preempt

    Int rsem_memory
    Int rsem_disk
    Int rsem_threads
    Int rsem_preempt

    Int rnaseqc_memory
    Int rnaseqc_disk
    Int rnaseqc_preempt

    call samtofastq {
        input: input_bam=input_bam, prefix=prefix, memory=samtofastq_memory, disk_space=samtofastq_disk, num_preempt=samtofastq_preempt
    }

    call star {
        input: fastq1=samtofastq.fastq1, fastq2=samtofastq.fastq2, prefix=prefix, star_index=star_index, memory=star_memory, disk_space=star_disk, num_threads=star_threads, num_preempt=star_preempt
    }

    call bamsync {
        input: source_bam=input_bam, target_bam=star.bam_file, target_bam_index=star.bam_index, prefix=prefix, memory=bamsync_memory, disk_space=bamsync_disk, num_preempt=bamsync_preempt
    }

    call rsem {
        input: transcriptome_bam=star.transcriptome_bam, rsem_reference=rsem_reference, prefix=prefix, memory=rsem_memory, disk_space=rsem_disk, num_threads=rsem_threads, num_preempt=rsem_preempt
    }
    
    call rnaseqc {
        input: bam_file=bamsync.patched_bam_file, bam_index=bamsync.patched_bam_index, genes_gtf=genes_gtf, genome_fasta=genome_fasta, genome_fasta_index=genome_fasta_index, prefix=prefix, notes=rnaseqc_notes, memory=rnaseqc_memory, disk_space=rnaseqc_disk, num_preempt=rnaseqc_preempt
    }
}
