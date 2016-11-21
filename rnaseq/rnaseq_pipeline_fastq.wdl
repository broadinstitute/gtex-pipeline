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
    
    Int memory
    Int disk_space
    Int num_threads
    
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
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
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
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
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

    Int memory
    Int disk_space

    command {
        touch ${bam_index}
        /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx${memory}g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
        -s ${prefix},${bam_file},${notes} -t ${genes_gtf} -r ${genome_fasta} -o . -noDoC -strictMode 
        mv ${prefix}/${prefix}.transcripts.rpkm.gct ${prefix}/${prefix}.gene_rpkm.gct
        python3 /src/convert_rnaseqc_counts.py ${prefix}/${prefix}.gene_rpkm.gct ${prefix}/${prefix}.exon_intron_report.txt
    }

    output {
        File rnaseqc_gene_rpkm = "${prefix}/${prefix}.gene_rpkm.gct"
        File rnaseqc_gene_counts = "${prefix}.gene_reads.gct"
        File rnaseqc_exon_counts = "${prefix}/${prefix}.exon_report.txt"
        File rnaseqc_metrics = "metrics.tsv"
        # Array[File] rnaseqc_report_files = ["report.html", "meanCoverageNorm_high.png", "meanCoverageNorm_high.txt", "meanCoverageNorm_low.png", "meanCoverageNorm_low.txt", "meanCoverageNorm_medium.png", "meanCoverageNorm_medium.txt", "meanCoverage_high.png", "meanCoverage_high.txt", "meanCoverage_low.png", "meanCoverage_low.txt", "meanCoverage_medium.png", "meanCoverage_medium.txt"]
        Array[File] rnaseqc_outputs = ["${prefix}/${prefix}.chimericPairs.txt", "${prefix}/${prefix}.exon_intron_report.txt", "${prefix}/${prefix}.intron_report.txt", "${prefix}/${prefix}.introns.rpkm.gct", "${prefix}/${prefix}.libraryComplexity.txt", "${prefix}/${prefix}.metrics.tmp.txt", "${prefix}/${prefix}.metrics.txt", "${prefix}/${prefix}.rRNA_counts.txt", "${prefix}/${prefix}.totalExonQuantifiedReads.txt"]
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseq_pipeline_fastq_workflow {

    File fastq1
    File fastq2
    String prefix

    File star_index
    File rsem_reference

    File genes_gtf
    File genome_fasta
    File genome_fasta_index
    String rnaseqc_notes

    Int star_memory
    Int star_disk
    Int star_threads

    Int rsem_memory
    Int rsem_disk
    Int rsem_threads

    Int rnaseqc_memory
    Int rnaseqc_disk


    call star {
        input: fastq1=fastq1, fastq2=fastq2, prefix=prefix, star_index=star_index, memory=star_memory, disk_space=star_disk, num_threads=star_threads
    }

    call rsem {
        input: transcriptome_bam=star.transcriptome_bam, rsem_reference=rsem_reference, prefix=prefix, memory=rsem_memory, disk_space=rsem_disk, num_threads=rsem_threads
    }
    
    call rnaseqc {
        input: bam_file=star.bam_file, bam_index=star.bam_index, genes_gtf=genes_gtf, genome_fasta=genome_fasta, genome_fasta_index=genome_fasta_index, prefix=prefix, notes=rnaseqc_notes, memory=rnaseqc_memory, disk_space=rnaseqc_disk
    }
}
