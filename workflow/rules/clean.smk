rule clean_fasta:
    input:
        "results/fasta/{sample}_v6.pseudohap.fasta.gz"
    output:
        "results/fasta/{sample}_v7.pseudohap.fasta.gz"
    params:
        prefix = lambda w: config['scaff_prefix'][w.sample]
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -s -r -p '^N+$' -v {input} | \
        seqkit sort --by-length --two-pass | \
        seqkit replace -p .+ -r "{params.prefix}_{nr}" --nr-width 5 | \
        gzip > {output}
        """