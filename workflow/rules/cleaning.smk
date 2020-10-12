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
        seqkit grep -f {input}_filt_list | gzip -c > {output}
        rm {input}_tmp {input}_filt_list
        """