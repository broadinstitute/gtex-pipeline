version 1.0

task CrossCheckSample {
    input {
        File first
        File first_index
        File second
        File second_index
        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String gatkTag="4.2.4.0"
    }

    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size(hapMap, "GB"))

    parameter_meta {
        first: {
            localization_optional: true
        }
        first_index: {
            localization_optional: true
        }
        second: {
            localization_optional: true
        }
        second_index: {
            localization_optional: true
        }
    }

    command <<<
        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{first} \
            -SI ~{second} \
            -H ~{hapMap} \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --CROSSCHECK_BY FILE \
            --EXPECT_ALL_GROUPS_TO_MATCH \
            --OUTPUT sample.crosscheck_metrics 
    >>>
    output {
        File metrics="sample.crosscheck_metrics"
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }
}

task IdentifySample {
    input {
        File bam
        File bam_index
        File vcf

        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String gatkTag="4.2.4.0"
    }

    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size([hapMap, vcf], "GB"))

    parameter_meta {
        bam: {
            localization_optional: true
        }
        bam_index: {
            localization_optional: true
        }
    }

    command <<<
        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{bam} \
            -SI ~{vcf} \
            -H ~{hapMap} \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --CROSSCHECK_BY SAMPLE \
            --OUTPUT sample.crosscheck_metrics 
    >>>
    output {
        File metrics="sample.crosscheck_metrics"
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }
}



workflow CrosscheckFingerprints {
    call CrossCheckSample


    output {
        File fingerprinting_metrics=CrossCheckSample.metrics

    }
}