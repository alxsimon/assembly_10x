rule lr_mkref:
    input:
        fa = "results/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        directory("results/fasta/refdata-{sample}_v2.pseudohap")
    log:
        "logs/lr_mkref.{sample}.log",
    container:
        "containers/supernova.sif"
    shell:
        """
        longranger mkref {input.fa} > {log} 2>&1
        mv refdata-{wildcards.sample}_v2.pseudohap {output}
        """

rule lr_align:
    input:
        unpack(get_fastq),
        directory("results/fasta/refdata-{sample}_v2.pseudohap")
    output:
        directory("results/purge_dups/lr_align_{sample}_v2"),
        "results/purge_dups/lr_align_{sample}_v2/outs/possorted_bam.bam"
    params:
        input_dir = lambda w, input: os.path.dirname(input[0]),
        run_id = lambda w: f'{w.sample}_v2',
        sample = lambda w, input: re.sub("_S.+_L.+_R1_001.fastq.gz", "", os.path.basename(input[0])),
        mem = config['supernova_mem']
    threads: 
        workflow.cores
    log:
        "logs/lr_align.{sample}.log"
    container:
        "containers/supernova.sif"
    shell:
        """
        longranger align \
        --id={params.run_id} \
        --fastqs={params.input_dir} \
        --sample={params.sample} \
        --reference={input[2]} \
        --localcores={threads} \
        --localmem={params.mem} \
        > {log} 2>&1
        mv {params.run_id} {output}
        """

rule ngscstat:
    input:
        "results/purge_dups/lr_align_{sample}_v2/outs/possorted_bam.bam"
    output:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov")
    params:
        workdir = lambda w, input: os.path.dirname(input[0]),
        input = lambda w, input: os.path.basename(input[0])
    log:
        "logs/purge_stats_ngsstat.{sample}.log"
    shell:
        """
        cd {params.workdir}
        ngscstat {params.input} 2> ../../../{log}
        """

rule calcuts:
    input:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov")
    output:
        "results/purge_dups/{sample}/cutoffs"
    log:
        "logs/purge_stats_calcuts.{sample}.log"
    shell:
        """
        calcuts {input[0]} > {output} 2> {log}
        """

rule split_fa:
    input: 
        "results/fasta/{sample}_v2.pseudohap.fasta.gz"
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
        "logs/self_map.{sample}.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        16
    shell:
        "(minimap2 -t {threads} "
        "-xasm5 -DP {input} {input} "
        "| gzip -c - > {output}) 2> {log}"

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
        fa = "results/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        purged = "results/purge_dups/{sample}/{sample}_v2.pseudohap.purged.fa.gz",
        haps = "results/purge_dups/{sample}/{sample}_v2.pseudohap.hap.fa.gz"
    params:
        prefix = lambda w, output: output[0].replace(".purged.fa.gz", "")
    shell:
        "get_seqs {input.bed} {input.fa} "
        "-p {params.prefix}; "
        "gzip {params.prefix}.*.fa"

