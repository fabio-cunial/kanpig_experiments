version 1.0


#
workflow TruvariIntersample {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }

    call TruvariIntersampleImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi
    }
    
    output {
        File bcftools_merged_vcf_gz = TruvariIntersampleImpl.bcftools_merged_vcf_gz
        File bcftools_merged_tbi = TruvariIntersampleImpl.bcftools_merged_tbi
        File truvari_collapsed_vcf_gz = TruvariIntersampleImpl.truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi = TruvariIntersampleImpl.truvari_collapsed_tbi
    }
}


task TruvariIntersampleImpl {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }
    
    Int ram_size_gb = 64  # Arbitrary
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))        
        
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        
        # BCFTOOLS
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z tmp1.vcf.gz > bcftools_merged.vcf.gz
        tabix -f bcftools_merged.vcf.gz
        rm -f tmp1.vcf*
        
        # TRUVARI
        ${TIME_COMMAND} truvari collapse --input bcftools_merged.vcf.gz --sizemin 0 --sizemax 1000000 --keep common --gt all > tmp1.vcf
        ${TIME_COMMAND} bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z tmp1.vcf > truvari_collapsed.vcf.gz
        tabix -f truvari_collapsed.vcf.gz
        rm -f tmp1.vcf*
    >>>
    
    output {
        File bcftools_merged_vcf_gz = work_dir + "/bcftools_merged.vcf.gz"
        File bcftools_merged_tbi = work_dir + "/bcftools_merged.vcf.gz.tbi"
        File truvari_collapsed_vcf_gz = work_dir + "/truvari_collapsed.vcf.gz"
        File truvari_collapsed_tbi = work_dir + "/truvari_collapsed.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}