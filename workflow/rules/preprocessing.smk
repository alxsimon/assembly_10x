
rule rm_dup_pcr:
    input:
        unpack(get_fastq)
    output:
        temp(multiext("results/preprocessing/{sample}/{sample}_dedup",
            "_R1.fastq.gz", "_R2.fastq.gz"))
    log:
        "logs/nubeam-dedup.{sample}.log"
    shell:
        """
        nubeam-dedup -i1 {input.fq1} -i2 {input.fq2} \
        -o1 {output[0]} -o2 {output[1]} \
        > {log} 2>&1
        """


rule proc10x_process:
    input:
        multiext("results/preprocessing/{sample}/{sample}_dedup",
            "_R1.fastq.gz", "_R2.fastq.gz")
    output:
        protected(multiext("results/preprocessing/{sample}/{sample}_dedup_proc",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")),
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


rule fastp:
    input:
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
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
        fastp -i {input[0]} -I {input[1]} \
        -o {output[0]} -O {output[1]} \
        --disable_length_filtering \
        --correction \
        --trim_poly_g \
        --json {params.report}.json \
        --html {params.report}.html \
        -w {threads}
        """


rule proc10x_filter_regen:
    input:
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz"),
        barcodes = "results/preprocessing/{sample}/{sample}_filt_barcodes.txt"
    output:
        protected(multiext("results/preprocessing/{sample}/{sample}_regen",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz"))
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
        -1 {input[0]} -2 {input[1]} |
        /opt/proc10xG/regen_10xReads.py \
        -o {params.out_prefix} \
        > {log} 2>&1
        """
