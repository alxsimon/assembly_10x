rule lr_mkref:
    input:
        fa = "results/fasta/{sample}_v2.pseudohap.fasta.gz"
    output:
        directory("results/purge_dups/{sample}/refdata-{sample}_v2.pseudohap"),
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.fa"
    params:
        tmp_fa = lambda w: f'results/purge_dups/{w.sample}/{w.sample}_v2.pseudohap.fa'
    log:
        "logs/lr_mkref.{sample}.log",
    container:
        "containers/supernova.sif"
    shell:
        """
        zcat {input.fa} > {params.tmp_fa}
        longranger mkref {params.tmp_fa} > {log} 2>&1
        mv refdata-{wildcards.sample}_v2.pseudohap {output[0]}
        """

rule lr_align:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_S1_L001_{R}_001.fastq.gz", R=["R1", "R2"]),
        rules.lr_mkref.output
    output:
        directory("results/purge_dups/{sample}/lr_align_{sample}_v2"),
        "results/purge_dups/{sample}/lr_align_{sample}_v2/outs/possorted_bam.bam"
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
        
        rm -r {output[0]} \
        && cp -al {params.run_id} {output[0]} \
        && rm -r {params.run_id}
        """

rule sort_by_name:
    input:
        "results/purge_dups/{sample}/lr_align_{sample}_v2/outs/possorted_bam.bam"
    output:
        temp("results/purge_dups/{sample}/namesorted_bam.bam")
    threads:
        16
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/purge_dups_samtools_namesort.{sample}.log"
    shell:
        "samtools sort -n -O BAM -o {output} {input} > {log} 2>&1"

rule ngscstat:
    input:
        "results/purge_dups/{sample}/namesorted_bam.bam"
    output:
        multiext("results/purge_dups/{sample}/TX", ".stat", ".base.cov")
    params:
        workdir = lambda w, output: os.path.dirname(output[0]),
        input = "namesorted_bam.bam"
    log:
        "logs/purge_stats_ngscstat.{sample}.log"
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

rule make_hist:
    input:
        "results/purge_dups/{sample}/cutoffs",
        "results/purge_dups/{sample}/TX.stat"
    output:
        "results/purge_dups/{sample}/hist_cutoffs.png"
    shell:
        "/opt/purge_dups/scripts/hist_plot.py "
        "-c {input[0]} "
        "{input[1]} "
        "{output}"

rule split_fa:
    input: 
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.fa"
    output:
        temp("results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fa")
    log:
        "logs/split_fa.{sample}.log"
    shell:
        "split_fa {input} > {output} 2> {log}"

rule self_map:
    input:
        "results/purge_dups/{sample}/{sample}_v2.pseudohap.split.fa"
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
    params:
        M = config['purge_dups']['M'],
        E = config['purge_dups']['E']
    log:
        "logs/purge_dups.{sample}.log"
    shell:
        "purge_dups -2 -M{params.M} -E{params.E} -T {input.cutoffs} "
        "-c {input.basecov} {input.selfmap} > {output} "
        "2> {log}"

rule get_sequences:
    input:
        bed = "results/purge_dups/{sample}/{sample}.dups.bed",
        fa = "results/purge_dups/{sample}/{sample}_v2.pseudohap.fa",
        hist = "results/purge_dups/{sample}/hist_cutoffs.png"
    output:
        purged = "results/purge_dups/{sample}/{sample}_v2.pseudohap.purged.fa.gz",
        haps = "results/purge_dups/{sample}/{sample}_v2.pseudohap.hap.fa.gz"
    params:
        prefix = lambda w, output: output[0].replace(".purged.fa.gz", "")
    shell:
        """
        get_seqs {input.bed} {input.fa} \
        -p {params.prefix}
        gzip {params.prefix}.*.fa
        """

