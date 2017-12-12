task ase_aggregate_by_individual {

    Array[File] ase_readcount_files
    Array[String] sample_ids
    Array[String] tissue_site_details
    File het_vcf
    File vep_dict
    File simulation_bias
    File mappability_bigwig
    File tissue_abbreviations
    File lamp_values
    String individual_id

    Int? coverage_cutoff
    Float? other_ratio_cutoff
    Float? mono_cutoff

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        # workaround for broken 'write_lines'
        python3 <<CODE
with open("ase_readcount_file_paths.tsv", "w") as f:
    f.write('sample_id\ttissue_site_detail\tase_readcount_file\n')
    for i,t,p in zip('${sep="," sample_ids}'.split(","), '${sep="," tissue_site_details}'.split(","), '${sep="," ase_readcount_files}'.split(",")):
        f.write('\t'.join([i,t,p])+'\n')
CODE
        python3 -u /src/ase_aggregate_by_individual.py ase_readcount_file_paths.tsv ${het_vcf} ${vep_dict} ${simulation_bias} ${mappability_bigwig} ${tissue_abbreviations} ${lamp_values} ${individual_id}
    }

    output {
        File ase_table = "${individual_id}.ase_table.tsv.gz"
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


workflow ase_aggregate_by_individual_workflow {
    call ase_aggregate_by_individual
}
