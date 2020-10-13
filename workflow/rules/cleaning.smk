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
    conda: 
        "../envs/phyloligo.yaml"
    threads: 
        workflow.cores
    shell:
        """
        phyloligo.py -i {input} -d JSD --method joblib -c {threads} -o {output}
        """

checkpoint clustering:
    input: 
        fa = "results/phyloligo/{sample}/{sample}.20p.fa",
        mat = "results/phyloligo/{sample}/{sample}.JSD.mat"
    output: 
        directory("results/phyloligo/{sample}/{sample}_clust")
    conda: 
        "../envs/phyloligo.yaml"
    shell:
        """
        phyloselect.py -i {input.mat} -m hdbscan -f {input.fa} -t -o {output}
        """

# rule blast_makedb:
#     input:
#         "resources/Fraisse2016_contigs.fasta"
#     output:
#         "resources/Fraisse2016_contigs.fasta.ndb"
#     conda:
#         "../envs/blast.yaml"
#     shell:
#         """
#         makeblastdb -in {input} -dbtype nucl -title 'Fraisse2016_contigs'
#         """

# rule blast_clust:
#     input:
#         db = "resources/Fraisse2016_contigs.fasta.ndb",
#         fa = "results/phyloligo/{sample}/{sample}_clust/data_fasta_{cl}.fa"
#     output:
#         "results/phyloligo/{sample}/{sample}_clust/data_blastn_{cl}.txt"
#     params:
#         db = lambda w, input: re.sub('\.ndb', '', input.db)
#     conda:
#         "../envs/blast.yaml"
#     shell:
#         """
#         blastn -query {input.fa} -db {params.db} -outfmt 6 > {output}
#         """


# rule Kount:
#     input: 
#         asm = "results/phyloligo/{sample}/{sample}_v5.cleaned.fa",
#         host = "...",
#         conta = "..."
#     output:
#         "{sample}_v5.cleaned.fa.mcp_hostwindows_vs_host_{sample}_host.fa_KL.dist",
#         "{sample}_v5.cleaned.fa.mcp_hostwindows_vs_conta_{sample}_{cl}.fa_KL.dist"
#     params:
#         dist = "JSD",
#         wd = lambda w: f'results/phyloligo/{w.sample}'
#     conda: 
#         "../envs/phyloligo.yaml"
#     threads:
#         config['phyloligo']['threads']
#     shell:
#         """
#         Kount.py -i {input.asm} -r {input.host} -c {input.conta} \
#         -u {threads} -d {params.dist} -W {params.wd}
#         """

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
        dist = "JSD",
        wd = lambda w: f'results/phyloligo/{w.sample}/contalocate/'
    conda: 
        "../envs/phyloligo.yaml"
    threads:
        workflow.cores
    script:
        "../scripts/recursive_decontamination.py"
