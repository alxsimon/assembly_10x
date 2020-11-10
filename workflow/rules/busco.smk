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
    if (w.version == "LNJA01" or w.version == "UYJE01") and w.sample == "gallo":
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