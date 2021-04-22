def get_fasta(w):
    if (w.version in ["GCA900618805", "GCA017311375"]):
        return ancient(f"resources/{w.version}.fasta.gz")
    else:
        return ancient(f"results/fasta/{w.sample}_{w.version}.pseudohap.fasta.gz")

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
        "results/stats/coruscus_GCA017311375.stats.json",
        "results/stats/gallo_GCA900618805.stats.json",
        expand("results/stats/{sample}_{version}.stats.json",
            sample=config['samples'], 
            version=["v1", "v2", "v3", "v4", "v5", "v6", "v7"]),
    output:
        "results/stats/assembly_stats.csv"
    conda:
        "../envs/asm_stats.yaml"
    script:
        "../scripts/merge_stats.py"

