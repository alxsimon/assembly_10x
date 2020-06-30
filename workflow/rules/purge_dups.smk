rule map_reads:
    input:
        "results/supernova_assemblies/{sample}_{version}/fasta/{sample}_v2.pseudohap.fasta.gz",
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp_filt",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
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
        zcat {input[0]} > {params.ref}
        bwa index {params.ref}
        bwa mem -t {threads} {params.ref} {input[1]} {input[2]} | \
        samtools view -b -@ {threads} - > {output} \
        2> {log}
        """

rule ngscstat:
    input:
        "results/purge_dups/{sample}/{sample}.bam"
    output:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov")
    params:
        workdir = lambda w, input: os.path.dirname(input[0]),
        input = lambda w, input: os.path.basename(input[0])
    log:
        "logs/ngscstat.{sample}.log"
    conda:
        "../envs/purge_dups.yaml"
    threads:
        16
    shell:
        """
        cd {params.workdir}
        ngscstat {params.input}
        """