def get_fastq(wildcards):
    prefix = config["raw_names"][wildcards.sample]
    fq1 = glob.glob("resources/10x_reads/" +
        prefix + '*_R1*.fastq.gz')
    fq2 = glob.glob("resources/10x_reads/" +
        prefix + '*_R2*.fastq.gz')
    return {'fq1': fq1, 'fq2': fq2}


rule rm_dup_pcr:
    input:
        unpack(get_fastq)
    output:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_{R}.fastq.gz", R=["R1", "R2"])
    log:
        "results/preprocessing/{sample}/nubeam-dedup.{sample}.log"
    shell:
        """
        nubeam-dedup -i1 {input.fq1} -i2 {input.fq2} \
        -o1 {output[0]} -o2 {output[1]} \
        > {log} 2>&1
        """


rule proc10x_process:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_{R}.fastq.gz", R=["R1", "R2"])
    output:
        "results/preprocessing/{sample}/{sample}_dedup_proc_barcodes.txt",
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_proc_{R}.fastq.gz", R=["R1", "R2"])
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_proc"
    log:
        "results/preprocessing/{sample}/process_10xReads.{sample}.log"
    conda:
        "envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/process_10xReads.py \
        -o {params.out_prefix} \
        -t 0 \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """

rule filter_barcodes:
    Here goes the R script


rule proc10x_filter:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_process_{R}.fastq.gz", R=["R1", "R2"]),
        barcodes = ""
    output:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_filt_{R}.fastq.gz", R=["R1", "R2"])
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_filt"
    log:
        "results/preprocessing/{sample}/filter_10xReads.{sample}.log"
    conda:
        "envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/filter_10xReads.py \
        -o {output} \
        -L {input.barcodes} \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """


rule proc10x_regen:
    input:
        barcodes = ""
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_filt_{R}.fastq.gz", R=["R1", "R2"])
    output:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_regen_{R}.fastq.gz", R=["R1", "R2"])
    params:
        out_prefix = "results/preprocessing/{sample}/{sample}_dedup_regen"
    log:
        "results/preprocessing/{sample}/regen_10xReads.{sample}.log"
    conda:
        "envs/py2.yaml"
    shell:
        """
        /opt/proc10xG/regen_10xReads.py \
        -o {params.out_prefix} \
        -1 {input[0]} -2 {input[1]} \
        > {log} 2>&1
        """
    