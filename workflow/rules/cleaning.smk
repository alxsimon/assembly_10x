rule split_on_btk_info:
    input:
        "results/blobtoolkit/blobdirs/{sample}_v5"
    output:
        kept = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_kept.csv",
        bact = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_bacteria.csv",
        virus = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_viruses.csv",
        euka = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_eukaryota.csv"
    params:
        mollusca_taxids = "resources/mollusca_taxids.txt"
    script:
        "../scripts/btk_conta_extraction.py"

rule btk_filter_conta:
    input:
        "results/blobtoolkit/blobdirs/{sample}_v5",
        "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_kept.csv"
    output:
        directory("results/blobtoolkit/blobdirs/{sample}_v5_kept"),
        directory("results/blobtoolkit/blobdirs/{sample}_v5_conta")
    params:
        blobtools_bin = config['btk']['blobtools_path'],
        tmp_kept_list = lambda w, input: f'{input[0]}/tmp_kept_list.ids'
    conda:
        "../envs/btk_env.yaml"
    shell:
        """
        tail -n +2 {input[1]} | cut -d ',' -f 1 > {params.tmp_kept_list}

        {params.blobtools_bin} filter \
        --list {params.tmp_kept_list} \
        --output {output[0]} {input[0]}

        {params.blobtools_bin} filter \
        --list {params.tmp_kept_list} \
        --invert \
        --output {output[1]} {input[0]}

        rm {params.tmp_kept_list}
        """

# Remove contamination scaffolds
# Remove leading and trailing Ns, 
# then remove sequences < 1000bp and with > 90% N.
rule v6_clean_rename:
    input:
        fa = "results/fasta/{sample}_v5.pseudohap.fasta.gz",
        kept = "results/blobtoolkit/blobdirs/{sample}_v5/{sample}_kept.csv"
    output:
        "results/fasta/{sample}_v6.pseudohap.fasta.gz"
    params: 
        scaff_prefix = lambda w: config['scaff_prefix'][w.sample],
        nr_width = 5
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -f <(tail -n +2 {input.kept} | cut -d ',' -f 1) {input[0]} \
        | seqkit replace -is -p "^N+|N+$" -r "" > {input.fa}_tmp

        seqkit fx2tab -n -i --gc --length -B N {input.fa}_tmp \
        | awk '($2 >= 1000 && $4 < 90) {{print $1}}' > {input.fa}_filt_list
        seqkit grep -f {input.fa}_filt_list {input.fa}_tmp > {input.fa}_tmp2

        seqkit sort --by-length -2 --reverse {input.fa}_tmp2 \
        | seqkit replace -p .+ -r "{params.scaff_prefix}_{{nr}}" --nr-width {params.nr_width} \
        | gzip -c > {output}

        rm {input.fa}_tmp* {input.fa}_filt_list
        """