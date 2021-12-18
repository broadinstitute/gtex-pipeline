version 1.0


task IdentifySample {
    input {
        File sample
        File sample_index
        File vcf
        File vcf_index

        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String? gatkTag
    }
    String gatkTag_final = select_first([gatkTag, "4.2.4.0"])


    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size([hapMap, vcf], "GB"))

    parameter_meta {
        sample: {
            localization_optional: true
        }
        sample_index: {
            localization_optional: true
        }
        vcf: {
            localization_optional: true
        }
        vcf_index: {
            localization_optional: true
        }
    }

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{sample} \
            -SI ~{vcf} \
            -H ~{hapMap} \
            --CALCULATE_TUMOR_AWARE_RESULTS false \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --CROSSCHECK_BY SAMPLE \
            --OUTPUT sample.crosscheck_metrics \
            --EXIT_CODE_WHEN_MISMATCH 3
    >>>
    output {
        File metrics="sample.crosscheck_metrics"
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag_final
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
            continueOnReturnCode: [0,3]

    }
}


workflow IdentifySampleWF{
    call IdentifySample{}

    output {
        File fp_metrics=IdentifySample.metrics
    }
}