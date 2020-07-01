rule map_reads:
    input:
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp_filt",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz"),
        fa = "results/supernova_assemblies/{sample}_v2/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        "results/purge_dups/{sample}/{sample}.bam"
    log:
        "logs/mapping_purge.{sample}.log"
    conda:
        "../envs/purge_dups.yaml"
    threads:
        16
    shell:
        """
        bwa index {input.fa}
        bwa mem -t {threads} {input.fa} {input[0]} {input[1]} | \
        samtools view -b -@ {threads} - > {output[1]} \
        2> {log}
        """

rule purge_stats:
    input:
        "results/purge_dups/{sample}/{sample}.bam"
    output:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov"),
        "results/purge_dups/{sample}/cutoffs"
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
        "results/supernova_assemblies/{sample}_v2/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        temp("results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fasta")
    log:
        "logs/split_fa.{sample}.log"
    shell:
        "split_fa {input} > {output} 2> {log}"

rule self_map:
    input:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fasta"
    output:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.self.paf.gz"
    log:
        "logs/self_map.{sample}/log"
    threads:
        16
    shell:
        "minimap2 -t {threads} "
        "-xasm5 -DP {input} {input} "
        "| gzip -c - > {output} 2> {log}"

rule purge_dups:
    input:
        selfmap = "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.self.paf.gz",
        basecov = "results/purge_dups/{sample}/TX.base.cov",
        cutoffs = "results/purge_dups/{sample}/cutoffs"
    output:
        "results/purge_dups/{sample}/{sample}.dups.bed"
    log:
        "logs/purge_dups.{sample}.log"
    shell:
        "purge_dups -2 -T {input.cutoffs} "
        "-c {input.basecov} {input.selfmap} > {output} "
        "2> {log}"

rule get_sequences:
    input:
        bed = "results/purge_dups/{sample}/{sample}.dups.bed",
        fa = "results/supernova_assemblies/{sample}_v2/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        #"results/purge_dups/{sample}/{sample}_v2.pdups.fasta.gz"
        "TEST_purge_dups_{sample}"
    shell:
        #"get_seqs {input.bed} {input.fa};"
        "touch {output}"

