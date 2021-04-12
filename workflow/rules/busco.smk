rule dwld_busco_databases:
    output:
        directory('resources/busco_databases/metazoa_odb10'),
        directory('resources/busco_databases/mollusca_odb10')
    params:
        metazoa = config['busco']['metazoa'],
        mollusca = config['busco']['mollusca'], 
        outdir = 'resources/busco_databases'
    shell:
        """
        cd {params.outdir}
        wget -c {params.metazoa} -O - | tar -xz
        wget -c {params.mollusca} -O - | tar -xz
        """

def get_busco_input(w):
    if (w.version in ["GCA900618805", "GCA017311375"]):
        return ancient(f"resources/{w.version}.fasta.gz")
    else:
        return ancient(f"results/fasta/{w.sample}_{w.version}.pseudohap.fasta.gz")

rule unzip_fasta:
    input:
        get_busco_input
    output:
        temp("results/fasta/{sample}_{version}.fa")
    shell:
        "zcat {input} > {output}"

rule busco:
    input:
        "results/fasta/{sample}_{version}.fa",
        rules.dwld_busco_databases.output
    output:
        "results/busco/{sample}_{version}_{db}/run_{db}/short_summary.txt",
        "results/busco/{sample}_{version}_{db}/run_{db}/full_table.tsv"
    params:
        db = lambda w: f'resources/busco_databases/{w.db}',
        fa = lambda w, input: input[0].replace(".fasta.gz", ".fa"),
        outdir = lambda w: f'{w.sample}_{w.version}_{w.db}'
    wildcard_constraints:
        db = '\w+_\w+'
    log:
        "logs/busco.{sample}_{version}_{db}.log"
    threads:
        config['busco']['threads']
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -f -m genome -i {params.fa} -o {params.outdir} \
        -q -c {threads} \
        -l {params.db} > {log} 2>&1
        cp -r {params.outdir} results/busco/ && rm -r {params.outdir}
        """


rule summarize_busco:
    input:
        expand("results/busco/{sample}_{version}_{db}/short_summary.specific.{db}.{sample}_{version}_{db}.txt",
            sample=config['samples'], 
            version=["v1", "v2", "v3", "v4", "v5", "v6"], # "v7"
            db=["metazoa_odb10", "mollusca_odb10"]),
        expand("results/busco/tros_v7_{db}/short_summary.specific.{db}.tros_v7_{db}.txt",
            db=["metazoa_odb10", "mollusca_odb10"]),
        expand("results/busco/{pub}_{db}/short_summary.specific.{db}.{pub}_{db}.txt",
            pub=['coruscus_GCA017311375', 'gallo_GCA900618805'], 
            db=["metazoa_odb10", "mollusca_odb10"])
    output:
        "results/stats/busco_summary.tsv"
    conda:
        "../envs/busco.yaml"
    shell:
        """
        python workflow/scripts/summarize_buscos.py \
        -o {output} \
        --files {input}
        """