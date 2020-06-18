
rule rm_dup_pcr:
    input:
        unpack(get_fastq)
    output:
        temp(unpack(rmdup_expand("results/preprocessing/{sample}/{sample}_dedup")))
    log:
        "logs/nubeam-dedup.{sample}.log"
    shell:
        """
        nubeam-dedup -i1 {input.fq1} -i2 {input.fq2} \
        -o1 {output.fq1} -o2 {output.fq2} \
        > {log} 2>&1
        """


rule proc10x_process:
    input:
        unpack(rmdup_expand("results/preprocessing/{sample}/{sample}_dedup"))
    output:
        protected(unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_proc"))),
        protected("results/preprocessing/{sample}/{sample}_dedup_proc_barcodes.txt")
    params:
        out_prefix = lambda w, output: output[0].strip("_R1_001.fastq.gz")
    log:
        "logs/process_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/process_10xReads.py \
        -o {params.out_prefix} \
        -1 {input.fq1} -2 {input.fq2} \
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
        unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_proc")),
        barcodes = "results/preprocessing/{sample}/{sample}_filt_barcodes.txt"
    output:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_filt_{R}_001.fastq.gz", R=["R1", "R2"])
    params:
        out_prefix = lambda w, output: output[0].strip("_R1_001.fastq.gz")
    log:
        "logs/filter_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/filter_10xReads.py \
        -o {params.out_prefix} \
        -L {input.barcodes} \
        -1 {input.fq1} -2 {input.fq2} \
        > {log} 2>&1
        """


rule proc10x_regen:
    input:
        unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_filt"))
    output:
        protected(unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_regen")))
    params:
        out_prefix = lambda w, output: output[0].strip("_R1_001.fastq.gz")
    log:
        "logs/regen_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/regen_10xReads.py \
        -o {params.out_prefix} \
        -1 {input.fq1} -2 {input.fq2} \
        > {log} 2>&1
        """
    
