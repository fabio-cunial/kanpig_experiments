version 1.0


#
workflow Resolve {
    input {
        String sample_id
        File vcf_gz
        File tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call ResolveImpl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            tbi = tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File resolved_vcf_gz = ResolveImpl.resolved_vcf_gz
    	File resolved_tbi = ResolveImpl.resolved_tbi
    }
}


#
task ResolveImpl {
    input {
        String sample_id
        File vcf_gz
        File tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools norm --multiallelics -any --check-ref s --do-not-normalize --fasta-ref ~{reference_fa} ~{vcf_gz} > tmp1.vcf
        
        bcftools view --include "SVTYPE != 'BND'" --output-type z tmp1.vcf > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf
        
        ${TIME_COMMAND} python ~{docker_dir}/resolve.py tmp2.vcf.gz ~{reference_fa} > tmp3.vcf
        rm -f tmp2.vcf.gz*
        
        #${TIME_COMMAND} python ~{docker_dir}/inversion_guesser.py -i tmp3.vcf -o tmp4.vcf
        #rm -f tmp3.vcf
        
        bgzip -c tmp3.vcf > ~{sample_id}_resolved.vcf.gz
        tabix -f ~{sample_id}_resolved.vcf.gz
    >>>
    
    output {
    	File resolved_vcf_gz = "~{work_dir}/~{sample_id}_resolved.vcf.gz"
    	File resolved_tbi = "~{work_dir}/~{sample_id}_resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
