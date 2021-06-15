
rule freebayes:
    input:
        bams = "results/mapping/{sample}_v7.bam",
        ref = "results/fasta/{sample}_v7.pseudohap.fasta.gz",
    output:
        "results/calling/{sample}_v7.bcf"
    conda:
        "../envs/calling.yaml"
    shell:
        """
        freebayes \
        -f {input.ref} \
        --standard-filters \
        {input.bams} | \
        bcftools view -Ob -o {output}
        """