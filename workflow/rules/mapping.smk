rule index_refs:
    input:
        "results/fasta/{sample}_{version}.pseudohap.fasta.gz"
    output:
        multiext("results/fasta/{sample}_{version}.pseudohap.fasta.gz",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        "logs/bwa-mem2_index_{sample}_{version}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa-mem2 index {input} > {log} 2>&1
        """

rule map_final:
    input:
        fa = "results/fasta/{sample}_{version}.pseudohap.fasta.gz",
        index = multiext("results/fasta/{sample}_{version}.pseudohap.fasta.gz",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        fastqs = multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        "results/mapping/{sample}_{version}.bam",
        "results/mapping/{sample}_{version}.bam.bai"
    log:
        "logs/bwa-mem2_{sample}_{version}.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        config['mapping']['threads']
    shell:
        """
        bwa-mem2 mem -t {threads} {input.fa} \
        {input.fastqs} 2> {log} | \
        samtools sort -@ {threads} -o {output[0]} -
        samtools index {output[0]}
        """

rule mapping_stats:
    input:
        "results/mapping/{sample}_{version}.bam"
    output:
        "results/mapping/{sample}_{version}.stats",
    threads:
        config['mapping']['threads']
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools stats -@ {threads} {input} > {output[0]}
        """

rule mosdepth:
    input:
        "results/mapping/{sample}_{version}.bam"
    output:
        "results/mapping/{sample}_{version}.mosdepth.global.dist.txt",
        "results/mapping/{sample}_{version}.per-base.d4",
    log:
        "logs/mosdepth/mosdepth_{exp}_{sample}.log"
    params:
        prefix = lambda w, output: output[1].replace('.per-base.d4', ''),
        min_mapq = 10,
    threads:
        config['mapping']['threads']
    shell:
        """
        mosdepth \
        -t {threads} \
        --d4 \
        --mapq {params.min_mapq} \
        {params.prefix} \
        {input} \
        > {log} 2>&1
        """