version 1.0


# Resources used when run on a 16-virtual-core, 32-GB machine:
# [input, time, cost, CPU, RAM]
#
# single sample 32x CCS     1h      $ 1.26      27 %      2 GB
# single sample 32x ONT     3h      $ 1.45      98 %      2 GB
# truvari collapse 32x CCS  1h      $ 1.50      67 %     16 GB
# truvari collapse 32x ONT  1h50m   $ 2        200 %     17 GB
# bcftools merge 32x CCS    40m     $ 1.50     120 %     20 GB
# bcftools merge 32x ONT    1h      $ 1.50     300 %     20 GB
#
workflow SnifflesGenotyper {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call SnifflesGenotyperImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
        File regenotyped_sniffles = SnifflesGenotyperImpl.regenotyped_sniffles
        File regenotyped_sniffles_tbi = SnifflesGenotyperImpl.regenotyped_sniffles_tbi
    }
}


task SnifflesGenotyperImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "sniffles_regenotyped"
    Int disk_size_gb = 2*( 100 + ceil(size(reference_fa,"GB")) + 10*ceil(size(input_vcf_gz,"GB")) + ceil(size(alignments_bam,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --reference ~{reference_fa} --input ~{alignments_bam} --genotype-vcf ~{input_vcf_gz} --vcf ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_sniffles = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_sniffles_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: 2
        memory: "64GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
