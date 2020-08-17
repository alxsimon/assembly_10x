def get_fasta(w):
    if w.sample == "gallo" and w.version == 0:
        return "resources/GCA_001676915.1_ASM167691v1/GCA_001676915.1_ASM167691v1_genomic.fna.gz"
    else:
        return f"results/fasta/{w.sample}_{w.version}.pseudohap.fasta.gz"

rule assembly_stats:
    input:
        get_fasta
    output:
        "results/stats/{sample}_{version}.stats.json"
    conda:
        "../envs/asm_stats.yaml"
    shell:
        """
        zcat {input} | \
        assembly_stats /dev/stdin \
        > {output}
        """

rule merge_stats:
    input:
        expand("results/stats/{sample}_{version}.stats.json",
            sample=config['samples'], version=["v1", "v2", "v3"]),
        "results/stats/gallo_v0.stats.json"
    output:
        "results/stats/assembly_stats.csv"
    shell:
        "touch {output}"