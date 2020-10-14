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
# Contamination removal with phyloligo
rule prepare_fasta:
    input: 
        "results/fasta/{sample}_v5.cleaned.fasta.gz"
    output: 
        "results/phyloligo/{sample}/{sample}_v5.cleaned.fa",
        "results/phyloligo/{sample}/{sample}.25p.fa"
    conda: 
        "../envs/phyloligo.yaml"
    shell:
        """
        zcat {input} > {output[0]}
        phylopreprocess.py -i {output[0]} -g 25 -r -o {output[1]}
        """

rule distance_matrix:
    input: 
        "results/phyloligo/{sample}/{sample}.25p.fa"
    output: 
        "results/phyloligo/{sample}/{sample}.JSD.mat"
    params:
        dist = config['phyloligo']['dist'],
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
        fa = "results/phyloligo/{sample}/{sample}.25p.fa",
        mat = "results/phyloligo/{sample}/{sample}.JSD.mat"
    output: 
        directory("results/phyloligo/{sample}/{sample}_clust")
    conda: 
        "../envs/phyloligo.yaml"
    shell:
        """
        phyloselect.py -i {input.mat} -m hdbscan -f {input.fa} -t -o {output}
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
        dist = config['phyloligo']['dist'],
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
