version 1.0


# Resources used when run on a 16-virtual-core, 32-GB machine:
# [input, time, cost, CPU, RAM]
#
# single sample 32x CCS      7h      $ 5.95      1500 %    22 GB
# single sample 32x ONT     10h      $ 9.12       900 %    22 GB
# truvari collapse 32x CCS  38h      $ 28        1400 %    28 GB
# truvari collapse 32x ONT  
# bcftools merge 32x CCS    73h      $ 60        1400 %    38 GB
# bcftools merge 32x ONT
#
workflow SvjediGenotyper {
    input {
        File input_vcf_gz
        File input_tbi
        File reads_fastq
        File reference_fa
        File reference_fai
        Int min_support = 3
        Int n_cpu = 16
        Int ram_size_gb = 64
    }
    parameter_meta {
        min_support: "Minimum number of total reads in the region. Default=3."
    }
    
    call SvjediGenotyperImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            reads_fastq = reads_fastq,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            min_support = min_support,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb
    }
    
    output {
        File regenotyped_jedi = SvjediGenotyperImpl.regenotyped_jedi
        File regenotyped_jedi_tbi = SvjediGenotyperImpl.regenotyped_jedi_tbi
    }
}


task SvjediGenotyperImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File reads_fastq
        File reference_fa
        File reference_fai
        Int min_support
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "svjedi_regenotyped"
    Int disk_size_gb = 100 + ceil(size(reference_fa,"GB")) + 10*ceil(size(input_vcf_gz,"GB")) + ceil(size(reads_fastq,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        gunzip -c ~{input_vcf_gz} > input.vcf
        ${TIME_COMMAND} python ~{docker_dir}/svjedigraph/svjedi-graph.py --threads ${N_THREADS} --minsupport ~{min_support} --vcf input.vcf --ref ~{reference_fa} --reads ~{reads_fastq} --prefix ~{output_prefix}
        bgzip -c ~{output_prefix}_genotype.vcf > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
        
    >>>

    output {
        File regenotyped_jedi = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_jedi_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
