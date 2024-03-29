rule lr_mkref:
    input:
        fa = "results/fasta/{sample}_{version}.pseudohap.fasta.gz"
    output:
        directory("results/purge_dups/{sample}_{version}/refdata-{sample}_{version}.pseudohap"),
        "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.fa"
    params:
        tmp_fa = lambda w: f'results/purge_dups/{w.sample}_{w.version}/{w.sample}_{w.version}.pseudohap.fa'
    log:
        "logs/lr_mkref.{sample}_{version}.log",
    container:
        "containers/supernova.sif"
    shell:
        """
        zcat {input.fa} > {params.tmp_fa}
        longranger mkref {params.tmp_fa} > {log} 2>&1
        mv refdata-{wildcards.sample}_v1.pseudohap {output[0]}
        """

rule lr_align:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_S1_L001_{R}_001.fastq.gz", R=["R1", "R2"]),
        rules.lr_mkref.output
    output:
        directory("results/purge_dups/{sample}_{version}/lr_align_{sample}_{version}"),
        "results/purge_dups/{sample}_{version}/lr_align_{sample}_{version}/outs/possorted_bam.bam"
    params:
        input_dir = lambda w, input: os.path.dirname(input[0]),
        run_id = lambda w: f'{w.sample}_{w.version}',
        sample = lambda w, input: re.sub("_S.+_L.+_R1_001.fastq.gz", "", os.path.basename(input[0])),
        mem = config['supernova_mem']
    threads: 
        workflow.cores
    log:
        "logs/lr_align.{sample}_{version}.log"
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
        "results/purge_dups/{sample}_{version}/lr_align_{sample}_{version}/outs/possorted_bam.bam"
    output:
        temp("results/purge_dups/{sample}_{version}/namesorted_bam.bam")
    threads:
        16
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/purge_dups_samtools_namesort.{sample}_{version}.log"
    threads:
        config['purge_dups']['threads']
    shell:
        "samtools sort -@ {threads} -n -O BAM -o {output} {input} > {log} 2>&1"

rule ngscstat:
    input:
        "results/purge_dups/{sample}_{version}/namesorted_bam.bam"
    output:
        multiext("results/purge_dups/{sample}_{version}/TX", ".stat", ".base.cov")
    params:
        workdir = lambda w, output: os.path.dirname(output[0]),
        input = "namesorted_bam.bam"
    log:
        "logs/purge_stats_ngscstat.{sample}_{version}.log"
    shell:
        """
        cd {params.workdir}
        ngscstat {params.input} 2> ../../../{log}
        """

rule calcuts:
    input:
        multiext("results/purge_dups/{sample}_{version}/TX", ".stat", ".base.cov")
    output:
        "results/purge_dups/{sample}_{version}/cutoffs"
    log:
        "logs/purge_stats_calcuts.{sample}_{version}.log"
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

rule split_fa_purge_dups:
    input: 
        "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.fa"
    output:
        temp("results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.split.fa")
    log:
        "logs/split_fa.{sample}_{version}.log"
    shell:
        "split_fa {input} > {output} 2> {log}"

rule self_map:
    input:
        "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.split.fa"
    output:
        "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.split.self.paf.gz"
    log:
        "logs/self_map.{sample}_{version}.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        config['purge_dups']['threads']
    shell:
        "(minimap2 -t {threads} "
        "-xasm5 -DP {input} {input} "
        "| gzip -c - > {output}) 2> {log}"

rule purge_dups:
    input:
        selfmap = "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.split.self.paf.gz",
        basecov = "results/purge_dups/{sample}_{version}/TX.base.cov",
        cutoffs = "results/purge_dups/{sample}_{version}/cutoffs"
    output:
        "results/purge_dups/{sample}_{version}/{sample}.dups.bed"
    params:
        M = config['purge_dups']['M'],
        E = config['purge_dups']['E']
    log:
        "logs/purge_dups.{sample}_{version}.log"
    shell:
        "purge_dups -2 -M{params.M} -E{params.E} -T {input.cutoffs} "
        "-c {input.basecov} {input.selfmap} > {output} "
        "2> {log}"

rule get_sequences:
    input:
        bed = "results/purge_dups/{sample}_{version}/{sample}.dups.bed",
        fa = "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.fa",
        hist = "results/purge_dups/{sample}_{version}/hist_cutoffs.png"
    output:
        purged = "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.purged.fa.gz",
        haps = "results/purge_dups/{sample}_{version}/{sample}_{version}.pseudohap.hap.fa.gz"
    params:
        prefix = lambda w, output: output[0].replace(".purged.fa.gz", "")
    shell:
        """
        get_seqs {input.bed} {input.fa} \
        -p {params.prefix}
        gzip {params.prefix}.*.fa
        """

rule copy_output_v3:
    input:
        "results/purge_dups/{sample}_v2/{sample}_v2.pseudohap.purged.fa.gz"
    output:
        protected("results/fasta/{sample}_v3.pseudohap.fasta.gz")
    shell:
        "cp {input} {output}"

rule copy_output_v4:
    input:
        "results/purge_dups/{sample}_v1/{sample}_v1.pseudohap.purged.fa.gz"
    output:
        protected("results/fasta/{sample}_v4.pseudohap.fasta.gz")
    shell:
        "cp {input} {output}"
