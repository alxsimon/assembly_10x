
rule rm_dup_pcr:
    input:
        unpack(get_fastq)
    output:
        temp(expand("results/preprocessing/{{sample}}/{{sample}}_dedup_{R}.fastq", R=["R1", "R2"]))
    log:
        "logs/nubeam-dedup.{sample}.log"
    shell:
        """
        set +e
        nubeam-dedup -i1 {input.fq1} -i2 {input.fq2} \
        -o1 {output[0]} -o2 {output[1]} \
        > {log} 2>&1
        
        exitcode=$?
	echo $exitcode
        if [ $exitcode -eq 1 ]
        then
            exit 0
        else
            exit 0
        fi
        """


rule proc10x_process:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_{R}.fastq", R=["R1", "R2"])
    output:
        protected("results/preprocessing/{sample}/{sample}_dedup_proc_barcodes.txt"),
        protected(expand("results/preprocessing/{{sample}}/{{sample}}_dedup_proc_{R}_001.fastq.gz", R=["R1", "R2"]))
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_proc"
    log:
        "logs/process_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/process_10xReads.py \
        -o {params.out_prefix} \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """

rule filter_barcodes:
    input:
        in_barcodes = "results/preprocessing/{sample}/{sample}_dedup_proc_barcodes.txt"
    output:
        out_barcodes = "results/preprocessing/{sample}/{sample}_filt_barcodes.txt",
        figure = "results/preprocessing/{sample}/{sample}_barcode_plot.pdf"
    log:
        "logs/filter_barcodes.{sample}.log"
    script:
        "../scripts/filter_barcodes.R"


rule proc10x_filter:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_proc_{R}_001.fastq.gz", R=["R1", "R2"]),
        barcodes = "results/preprocessing/{sample}/{sample}_filt_barcodes.txt"
    output:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_filt_{R}_001.fastq.gz", R=["R1", "R2"])
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_filt"
    log:
        "logs/filter_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/filter_10xReads.py \
        -o {params.out_prefix} \
        -L {input.barcodes} \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """


rule proc10x_regen:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_filt_{R}_001.fastq.gz", R=["R1", "R2"])
    output:
        protected(expand("results/preprocessing/{{sample}}/{{sample}}_dedup_regen_{R}_001.fastq.gz", R=["R1", "R2"]))
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_regen"
    log:
        "logs/regen_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/regen_10xReads.py \
        -o {params.out_prefix} \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """
    
