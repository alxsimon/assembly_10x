
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


rule fastp:
    input:
        unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_proc"))
    output:
        unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_proc_fastp"))
    params:
        report = lambda w, output: os.path.dirname(output[0]) + f'/{w.sample}_fastp'
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp.{sample}.log"
    threads:
        8
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} \
        -o {output.fq1} -O {output.fq2} \
        --disable_length_filtering \
        --correction \
        --trim_poly_g \
        --json {params.report}.json \
        --html {params.report}.html \
        -w {threads}
        """


rule proc10x_filter_regen:
    input:
        unpack(proc10x_expand("results/preprocessing/{sample}/{sample}_dedup_proc_fastp")),
        barcodes = "results/preprocessing/{sample}/{sample}_filt_barcodes.txt"
    output:
        protected(expand("results/preprocessing/{{sample}}/{{sample}}_regen_{R}_001.fastq.gz",
            R=["R1", "R2"]))
    params:
        out_prefix = lambda w, output: output[0].strip("_R1_001.fastq.gz")
    log:
        "logs/filter_regen_10xReads.{sample}.log"
    conda:
        "../envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/filter_10xReads.py \
        -L {input.barcodes} \
        -1 {input.fq1} -2 {input.fq2} |
        /opt/proc10xG/regen_10xReads.py \
        -o {params.out_prefix} \
        > {log} 2>&1
        """    
