rule map_reads:
    input:
        "results/supernova_assemblies/{sample}_{version}/fasta/{sample}_v2.pseudohap.fasta.gz",
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp_filt",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.fasta",
        "results/purge_dups/{sample}/{sample}.bam"
    params:
        ref = lambda w, input: input[0].replace(".gz", "")
    log:
        "logs/mapping_purge.{sample}.log"
    conda:
        "../envs/purge_dups.yaml"
    threads:
        16
    shell:
        """
        zcat {input[0]} > {output[0]}
        bwa index {output[0]}
        bwa mem -t {threads} {params.ref} {input[1]} {input[2]} | \
        samtools view -b -@ {threads} - > {output[1]} \
        2> {log}
        """

rule purge_stats:
    input:
        "results/purge_dups/{sample}/{sample}.bam"
    output:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov")
    params:
        workdir = lambda w, input: os.path.dirname(input[0]),
        input = lambda w, input: os.path.basename(input[0])
    log:
        "logs/purge_stats.{sample}.log"
    conda:
        "../envs/purge_dups.yaml"
    threads:
        16
    shell:
        """
        cd {params.workdir}
        (ngscstat {params.input} && \
        calcuts TX.stats > cutoffs) 2> {log}
        """

rule split_fa:
    input: 
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.fasta"
    output:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fasta"
    log:
        "logs/split_fa.{sample}.log"
    shell:
        """
        split_fa {input} > {output} 2> {log}
        """

rule self_map:
    input:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fasta"
    output:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.self.paf.gz"
    log:
        "logs/self_map.{sample}/log"
    threads:
        ...
    shell:
        """
        minimap2 -xasm5 -DP {input} {input} | gzip -c - > {output} 2> {log}
        """

