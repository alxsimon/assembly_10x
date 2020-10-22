rule split_on_btk_info:
    input:
        "results/blobtoolkit/blobdirs/{sample}_v5"
    output:
        kept = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_kept.csv",
        bact = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_bacteria.csv",
        virus = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_viruses.csv",
        euka = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_eukaryota.csv"
    script:
        "../scripts/btk_conta_extraction.py"


# Remove leading and trailing Ns, 
# then remove sequences < 1000bp and with > 90% N.
rule clean_fasta:
    input:
        "results/fasta/{sample}_v5.pseudohap.fasta.gz"
    output:
        "results/fasta/{sample}_v5.cleaned.fasta.gz"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit replace -is -p "^N+|N+$" -r "" {input} > {input}_tmp
        seqkit fx2tab -n -i --gc --length -B N {input}_tmp \
        | awk '($2 > 1000 && $4 < 90) {{print $1}}' > {input}_filt_list
        seqkit grep -f {input}_filt_list {input}_tmp | gzip -c > {output}
        rm {input}_tmp {input}_filt_list
        """

#===============================
# Most obvious contamination removal

rule get_potential_conta:
    input:
        expand("results/blobtoolkit/blobdirs/{sample}_v5/bestsumorder_phylum.json",
            sample=config['samples'])
    output:
        "results/phyloligo/potential_conta_btk_phylum.json"
    run:
        phylums = []
        import json
        for file in input:
            with open(file) as fr:
                phylums += json.load(fr)['keys']
        phylums = [x for x in phylums if x not in ['Mollusca', 'no-hit']]
        phylums = list(set(phylums))
        with open(output[0], 'w') as fw:
            json.dump(phylums, fw, indent=2)

rule extract_potential_conta:
    input:
        "results/phyloligo/potential_conta_btk_phylum.json",
        "results/blobtoolkit/blobdirs/{sample}_v5",
        "results/fasta/{sample}_v5.cleaned.fasta.gz"
    output:
        "results/phyloligo/{sample}/{sample}.btk_conta.fa",
        "results/phyloligo/{sample}/{sample}.btk_mollusca.fa"
    conda:
        "../envs/seqkit.yaml"
    script:
        "../scripts/btk_extraction.py"

rule prepare_fasta:
    input: 
        "results/fasta/{sample}_v5.cleaned.fasta.gz",
        "results/phyloligo/{sample}/{sample}.btk_conta.fa",
        "results/phyloligo/{sample}/{sample}.btk_mollusca.fa"
    output: 
        "results/phyloligo/{sample}/{sample}_v5.cleaned.fa",
        "results/phyloligo/{sample}/{sample}.subsample.fa",
        "results/phyloligo/{sample}/{sample}.for_selection.fa"
    params: 
        perc = config['phyloligo']['perc_sampling']
    conda: 
        "../envs/phyloligo.yaml"
    shell:
        """
        zcat {input[0]} > {output[0]}
        phylopreprocess.py -i {output[0]} -g {params.perc} -r -o {output[1]}
        cat {output[1]} {input[1]} {input[2]} > {output[2]}
        """

rule distance_matrix:
    input: 
        "results/phyloligo/{sample}/{sample}.for_selection.fa"
    output: 
        "results/phyloligo/{sample}/{sample}.distmat"
    params:
        dist = config['phyloligo']['dist_matrix'],
        pattern = config['phyloligo']['pattern']
    conda: 
        "../envs/phyloligo.yaml"
    threads: 
        workflow.cores
    shell:
        """
        phyloligo.py -i {input} -d {params.dist} --method joblib \
        -c {threads} -o {output}  \
        -p {params.pattern}
        """

checkpoint clustering:
    input: 
        fa = "results/phyloligo/{sample}/{sample}.for_selection.fa",
        mat = "results/phyloligo/{sample}/{sample}.distmat"
    output: 
        directory("results/phyloligo/{sample}/{sample}_clust")
    conda: 
        "../envs/phyloligo.yaml"
    shell:
        """
        phyloselect.py -i {input.mat} -m hdbscan \
        -f {input.fa} -t -o {output}
        """

def get_clusters(wildcards):
    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
    return expand("results/phyloligo/{sample}/{sample}_clust/data_fasta_cl{num}.fa",
        sample=wildcards.sample,
        num=glob_wildcards(os.path.join(checkpoint_output, "data_fasta_cl{num}.fa")).num)

rule recursive_decontamination:
    input: 
        "results/phyloligo/{sample}/{sample}_v5.cleaned.fa",
        get_clusters
    output:
        "results/phyloligo/{sample}/contalocate/DONE",
        "results/fasta/{sample}_v6.pseudohap.fasta.gz",
        "results/fasta/{sample}_v6.contaminants.fasta.gz",
        "results/fasta/{sample}_v6.contaminants.gff"
    params:
        dist = config['phyloligo']['dist_contalocate'],
        wd = lambda w: f'results/phyloligo/{w.sample}/contalocate/',
        win_step = config['phyloligo']['win_step'],
        win_size = config['phyloligo']['win_size'],
        pattern = config['phyloligo']['pattern']
    conda: 
        "../envs/phyloligo.yaml"
    threads:
        workflow.cores
    script:
        "../scripts/recursive_decontamination.py"

rule btk_filter:
    input:
        "results/fasta/{sample}_v6.contaminants.gff",
        "results/blobtoolkit/blobdirs/{sample}_v5"
    output:
        directory("results/blobtoolkit/blobdirs/{sample}_v6"),
        directory("results/blobtoolkit/blobdirs/{sample}_v6_contam")
    params:
        blobtools_bin = config['btk']['blobtools_path']
    conda:
        "../envs/btk_env.yaml"
    shell:
        """
        cut -f 1 {input[0]} > tmp_filt_list_{wildcards.sample}

        {params.blobtools_bin} filter \
        --list tmp_filt_list_{wildcards.sample} \
        --invert --output {output[0]} {input[1]}

        {params.blobtools_bin} filter \
        --list tmp_filt_list_{wildcards.sample} \
        --output {output[1]} {input[1]}

        rm tmp_filt_list_{wildcards.sample}
        """