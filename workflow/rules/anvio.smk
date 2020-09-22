rule reformat_fasta:
    input:
        "results/fasta/{sample}_v5.pseudohap.fasta.gz"
    output:
        "results/anvio/{sample}/{sample}_fixed.fa"
    container: 
        "docker://meren/anvio:6.2"
    log:
        "logs/anvio_simplify_names_{sample}.log"
    shell: 
        """
        anvi-script-reformat-fasta \
        -l 1000 --simplify-names -o {output} <(zcat {input}) \
        > {log} 2>&1
        """

rule gen_db:
    input:
        "results/anvio/{sample}/{sample}_fixed.fa"
    output:
        "results/anvio/{sample}/{sample}.db"
    container: 
        "docker://meren/anvio:6.2"
    log:
        "logs/anvio_gen_db_{sample}.log"
    shell: 
        """
        anvi-gen-contigs-database -f {input} \
        -o {output} --skip-mindful-splitting \
        -n '{wildcards.sample}' \
        > {log} 2>&1
        """

rule run_hmms:
    input:
        "results/anvio/{sample}/{sample}.db"
    output:
        multiext("results/anvio/{sample}/{sample}.db", ".hits", ".genes")
    container: 
        "docker://meren/anvio:6.2"
    log:
        "logs/anvio_run_hmms_{sample}.log"
    threads: 
        config['anvio']['threads']
    shell: 
        """
        anvi-run-hmms -T {threads} -c {input} \
        > {log} 2>&1
        anvi-script-gen_stats_for_single_copy_genes.py {input}
        """

rule filter_fasta:
    input:
        "results/anvio/{sample}/{sample}.db.hits",
        "results/anvio/{sample}/{sample}_fixed.fa"
    output:
        "results/fasta/{sample}_v6.pseudohap.fasta.gz",
        "results/anvio/{sample}/{sample}_trashed_scaff.fa"
    conda: 
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -f <(tail -n +2 {input[0]} | cut -f 2) {input[1]} > {output[1]}
        seqkit grep -v -f <(tail -n +2 {input[0]} | cut -f 2) {input[1]} | \
        gzip > {output[0]}
        """