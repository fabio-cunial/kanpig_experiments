version 1.0


# 
workflow SnifflesSingleSample {
    input {
        String sample_id
        File input_bam
        File input_bai
        File reference_fa
        File tandems_bed
    }

    call SingleSampleCalling {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            reference_fa = reference_fa,
            tandems_bed = tandems_bed
    }

    output {
         File output_vcf_gz = SingleSampleCalling.vcf_gz
         File output_tbi = SingleSampleCalling.tbi
         File output_snf = SingleSampleCalling.snf
    }
}


task SingleSampleCalling {
    input {
        String sample_id
        File input_bam
        File input_bai
        File reference_fa
        File tandems_bed
    }
    
    Int disk_size_gb = 2*ceil(size(input_bam, "GB")) + ceil(size(reference_fa, "GB")) + 50
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"

        ${TIME_COMMAND} sniffles --threads ${N_THREADS} \
                 --input ~{input_bam} \
                 --reference ~{reference_fa} \
                 --tandem-repeats ~{tandems_bed} \
                 --sample-id ~{sample_id} \
                 --vcf ~{sample_id}.sniffles.vcf \
                 --snf ~{sample_id}.sniffles.snf
        bgzip ~{sample_id}.sniffles.vcf
        tabix -f ~{sample_id}.sniffles.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/" + sample_id + ".sniffles.vcf.gz"
        File tbi = work_dir + "/" + sample_id + ".sniffles.vcf.gz.tbi"
        File snf = work_dir + "/" + sample_id + ".sniffles.snf"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: 32
        memory: "64GB"  # Arbitrary
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
