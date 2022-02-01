def get_fasta(w):
    if (w.version in ["GCA900618805", "GCA017311375", "GCA905397895", "GCA019925415"]):
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
        "results/stats/edu_GCA905397895.stats.json",
        "results/stats/eduam_GCA019925415.stats.json",
        expand("results/stats/{sample}_{version}.stats.json",
            sample=config['samples'], 
            version=["v1", "v2", "v3", "v4", "v5", "v6", "v7"]),
    output:
        "results/stats/assembly_stats.csv"
    conda:
        "../envs/asm_stats.yaml"
    script:
        "../scripts/merge_stats.py"

#========================================
# Coverage histograms

rule d4_hist:
    input:
        "results/mapping/{sample}_{version}.per-base.d4"
    output:
        "results/stats/{sample}_{version}.cov.hist"
    shell:
        """
        d4tools stat -s hist {input} > {output}
        """
        # for now d4tools have a maxbin at 1000X

rule plot_coverage_hist:
    input:
        hist = "results/stats/{sample}_{version}.cov.hist"
    output:
        pdf = "results/stats/{sample}_{version}.cov.hist.pdf"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_coverage_hist.R"