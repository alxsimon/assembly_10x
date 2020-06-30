rule map_reads:
    input:
        "results/supernova_assemblies/{sample}_{version}/fasta/{sample}_v2.pseudohap.fasta.gz",
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp_filt",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz") # for now, would be better to use output of filter
    output:
        "results/purge_dups/{sample}.bam"
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