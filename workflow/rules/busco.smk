rule dwld_busco_databases:
    output:
        directory('resources/busco_databases/metazoa_odb10'),
        directory('resources/busco_databases/mollusca_odb10')
    params:
        metazoa = config['busco']['metazoa'],
        mollusca = config['busco']['mollusca']
    shell:
        "wget -P resources/busco_databases/ {params.metazoa} {params.mollusca} && "
        "tar -xf resources/busco_databases/*.tar.gz && "
        "rm resources/busco_databases/*.tar.gz"

rule unzip_fasta:
    input:
        "results/fasta/{sample}_{version}.pseudohap.fasta.gz"
    output:
        temp("results/fasta/{sample}_{version}.fa")
    wildcard_constraints:
        version = 'v[0-9]+'
    shell:
        "zcat {input} > {output}"

rule busco:
    input:
        "results/fasta/{sample}_{version}.pseudohap.fa",
        rules.dwld_busco_databases.output
    output:
        "results/busco/{sample}_{version}_{db}/run_{db}/short_summary.txt"
    params:
        db = lambda w: f'resources/busco_databases/{w.db}',
        fa = lambda w, input: input[0].replace(".fasta.gz", ".fa"),
        outdir = lambda w: f'{w.sample}_{w.version}_{w.db}'
    wildcard_constraints:
        db = '\w+_\w+',
        version = 'v[0-9]+'
    log:
        "logs/busco.{sample}_{version}_{db}.log"
    threads:
        16
    conda:
        "../envs/busco.yaml" 
    shell:
        """
        busco -f -m genome -i {params.fa} -o {params.outdir} \
        -q -c {threads} \
        -l {params.db} > {log} 2>&1
        cp -r {params.outdir} results/busco/ && rm -r {params.outdir}
        """