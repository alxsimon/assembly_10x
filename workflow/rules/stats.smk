def get_fasta(w):
    if w.sample == "gallo" and w.version == "v0":
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
            sample=config['samples'], version=["v1", "v2"]), # , "v3"]),
        "results/stats/gallo_v0.stats.json"
    output:
        "results/stats/assembly_stats.csv"
    conda:
        "../envs/asm_stats.yaml"
    run:
        import pandas as pd
        from pathlib import Path
        import json
        out_df = pd.DataFrame()
        for stat_file in input:
            with open(stat_file, 'r') as f:
                data = json.loads(f.read())
            data['C'] = data.pop('Contig Stats')
            data['S'] = data.pop('Scaffold Stats')
            tmpdf = pd.json_normalize(data)
            df.insert(0, "asm", [f"{wildcards.sample}_{wildcards.version}"])
            out_df = out_df.append(df)
        out_df.to_csv(output[0], sep=',', index=False)

