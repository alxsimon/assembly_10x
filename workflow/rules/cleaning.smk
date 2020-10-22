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