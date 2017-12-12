task combine_signif_pairs {

    Array[File] signifpairs
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/combine_signif_pairs.py ${write_lines(signifpairs)} ${prefix}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File combined_signifpairs="${prefix}.combined_signifpairs.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}


task extract_pairs {

    File input_pairs
    File extract_pairs
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/extract_pairs.py ${input_pairs} ${extract_pairs} ${prefix} --feather
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File extracted_pairs="${prefix}.extracted_pairs.ft"
    }

    meta {
        author: "Francois Aguet"
    }
}


task metasoft_prepare_input {

    Array[File] pair_files
    String prefix
    Int? chunks

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/metasoft_prepare_input.py ${write_lines(pair_files)} ${prefix} ${"--chunks " + chunks} --write_full
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File metasoft_input = "${prefix}.metasoft_input.txt.gz"
        Array[File] metasoft_input_chunks = glob("${prefix}.metasoft_input.chunk*.txt.gz")
    }

    meta {
        author: "Francois Aguet"
    }
}



task metasoft_scatter {

    File metasoft_input
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/run_metasoft.py /opt/metasoft/Metasoft.jar ${metasoft_input} ${prefix}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File metasoft_output="${prefix}.metasoft.txt.gz"
        File metasoft_log="${prefix}.metasoft.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


task metasoft_postprocess {

    Array[File] metasoft_chunks
    Array[File] pair_files
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        /src/metasoft_postprocess.py ${write_lines(metasoft_chunks)} ${write_lines(pair_files)} ${prefix}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File metasoft_output="${prefix}.metasoft.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow metasoft_workflow {

    Array[File] signifpairs
    Array[File] allpairs
    Array[String] sample_ids
    String prefix

    call combine_signif_pairs { input: signifpairs=signifpairs, prefix=prefix }

    # for each sample set, extract same set of pairs
    scatter(i in range(length(allpairs))) {
        call extract_pairs { input: input_pairs=allpairs[i], extract_pairs=combine_signif_pairs.combined_signifpairs, prefix=sample_ids[i]}
    }

    call metasoft_prepare_input { input: pair_files=extract_pairs.extracted_pairs, prefix=prefix}

    scatter(chunk in metasoft_prepare_input.metasoft_input_chunks) {
        call metasoft_scatter {
            input: metasoft_input=chunk, prefix=prefix
        }
    }

    call metasoft_postprocess {
        input: metasoft_chunks=metasoft_scatter.metasoft_output, pair_files=extract_pairs.extracted_pairs, prefix=prefix
    }

}
