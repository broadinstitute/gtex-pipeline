task fastqtl {

    File rpkm_gct
    File reads_gct
    File annotation_gtf
    File vcf
    File vcf_index	
    String prefix

    String? expression_threshold
    String? count_threshold
    String? min_samples

    String num_peer
    File genotype_pcs
    File? add_covariates

    String chunks
    String permutations
    String? cis_window
    String? ma_sample_threshold
    String? maf_threshold

    String fdr
    File variant_lookup

    Int memory
    Int disk_space
    Int num_threads


    command {
        set -euo pipefail
        touch ${vcf_index}  # avoid tabix "index older than vcf" error
        # pre-processing
        /src/normalize_expression.py ${rpkm_gct} ${reads_gct} ${annotation_gtf} ${vcf} ${prefix} ${"--expression_threshold " + expression_threshold} ${"--count_threshold " + count_threshold} ${"--min_samples " + min_samples}
        Rscript /src/run_PEER.R ${prefix}.expression.txt ${prefix} ${num_peer}		
        /src/combine_covariates.py ${genotype_pcs} ${prefix}_PEER_covariates.txt ${prefix} ${"--add_covariates " + add_covariates}
        # nominal pass
        /opt/fastqtl/python/run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} --covariates ${prefix}.combined_covariates.txt ${"--window " + cis_window} ${"--ma_sample_threshold " + ma_sample_threshold} ${"--maf_threshold " + maf_threshold} --chunks ${chunks} --threads ${num_threads}
        # 2nd nominal pass (FPKM)
        /opt/fastqtl/python/run_FastQTL_threaded.py ${vcf} ${prefix}.expression.fpkm.bed.gz ${prefix}.fpkm --covariates ${prefix}.combined_covariates.txt ${"--window " + cis_window} ${"--ma_sample_threshold " + ma_sample_threshold} ${"--maf_threshold " + maf_threshold} --chunks ${chunks} --threads ${num_threads}
        # permutation pass
        /opt/fastqtl/python/run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} --covariates ${prefix}.combined_covariates.txt --permute ${permutations} ${"--window " + cis_window} ${"--ma_sample_threshold " + ma_sample_threshold} ${"--maf_threshold " + maf_threshold} --chunks ${chunks} --threads ${num_threads}
        # post-processing
        /opt/fastqtl/python/annotate_outputs.py ${prefix}.egenes.txt.gz ${fdr} ${annotation_gtf} ${variant_lookup} --nominal_results ${prefix}.allpairs.txt.gz --nominal_results_unnormalized ${prefix}.fpkm.allpairs.txt.gz fpkm
    }

    runtime {
        docker: "broadinstitute/gtex_eqtl"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File covariates="${prefix}.combined_covariates.txt"
        File expression_bed="${prefix}.expression.bed.gz"
        File expression_bed_index="${prefix}.expression.bed.gz.tbi"
        File expression_fpkm_bed="${prefix}.expression.fpkm.bed.gz"
        File expression_fpkm_bed_index="${prefix}.expression.fpkm.bed.gz.tbi"
        File fastqtl_egenes="${prefix}.egenes.annotated.txt.gz"
        File fastqtl_egenes_log="${prefix}.egenes.log"
        File fastqtl_signifpairs="${prefix}.signifpairs.txt.gz"		
        File fastqtl_allpairs="${prefix}.allpairs.txt.gz"
        File fastqtl_allpairs_log="${prefix}.allpairs.log"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow fastqtl_workflow {
    call fastqtl
}
