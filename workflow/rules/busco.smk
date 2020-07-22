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

rule busco:
    input:
        glob.glob("results/fasta/{sample}_{version}*"),
        rules.dwld_busco_databases.output
    output:
        directory("results/busco/{sample}_{version}_{db}")
    params:
        db = lambda w: f'resources/busco_databases/{w.db}_odb10',
        fa = lambda w, input: input[0].replace(".fasta.gz", ".fa")
    wildcard_constraints:
        db = '\w_\w',
        version = 'v[0-9]+'
    log:
        "logs/busco.{sample}_{version}_{db}.log"
    threads:
        16
    conda:
        "../envs/busco.yaml" 
    shell:
        """
        zcat {input[0]} > {params.fa}
        
        busco -m genome -i {params.fa} -o {output} \
        -q -c {threads} \
        -l {params.db} > {log} 2>&1

        rm {params.fa}
        """