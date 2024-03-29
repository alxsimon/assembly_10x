def get_asm(w):
    if w.asm in ["gallo_v7", "tros_v7", "edu_v7"]:
        return f"results/fasta/{w.asm}.pseudohap.fasta.gz"
    elif w.asm in ["GCA017311375", "GCA900618805"]:
        return f"resources/{w.asm}.fasta.gz"

rule unzip_fasta_repeats:
    input:
        get_asm
    output:
        temp("results/repeats/{asm}.fa")
    wildcard_constraints:
        asm = '\w+_v7|\w+'
    shell:
        "zcat {input} > {output}"

rule build_db:
    input:
        "results/repeats/{asm}.fa"
    output:
        "results/repeats/{asm}_db/{asm}_db.nhr"
    params:
        db_name = lambda w, output: output[0].replace(".nhr", ""),
    container:
        config['dfam_container']
    log:
        "logs/repeats/build_db_{asm}.log"
    shell:
        """
        export LC_ALL=C
        BuildDatabase -name {params.db_name} {input} \
        > {log} 2>&1
        """

rule repeat_modeler:
    input:
        "results/repeats/{asm}_db/{asm}_db.nhr"
    output:
        "results/repeats/{asm}_db/{asm}_db-families.fa"
    params:
        pa = int(math.floor(config['repeats']['threads']/4)),
        db_name = lambda w, input: input[0].replace(".nhr", ""),
    threads:
        config['repeats']['threads']
    container:
        config['dfam_container']
    log:
        "logs/repeats/repeat_modeler_{asm}.log"
    shell:
        """
        export LC_ALL=C
        RepeatModeler -database {params.db_name} -LTRStruct -pa {params.pa} \
        > {log} 2>&1
        """

rule cdhit_families_merging:
    input:
        expand("results/repeats/{asm}_db/{asm}_db-families.fa",
            asm=config['repeats']['asm_db'])
    output:
        "results/repeats/Mytilus_sp_repeats-families.fasta"
    params:
        tmp_merge = "results/repeats/combined_db_families.fa",
    log:
        "logs/repeats/cdhit_merge.log"
    conda:
        "../envs/repeats.yaml"
    threads:
        config['repeats']['threads']
    shell:
        """
        cat {input} > {params.tmp_merge}

        cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 \
        -T {threads} \
        -i {params.tmp_merge} \
        -o {output} \
        > {log} 2>&1

        rm {params.tmp_merge}
        """

rule repeat_masker:
    input:
        families = "results/repeats/Mytilus_sp_repeats-families.fasta",
        fa = "results/repeats/{asm}.fa"
    output:
        "results/repeats/{asm}.fa.masked"
    params:
        pa = int(math.floor(config['repeats']['threads']/4)),
        out_dir = lambda w, output: os.path.dirname(output[0])
    threads:
        config['repeats']['threads']
    container:
        config['dfam_container']
    log:
        "logs/repeats/repeat_masker_{asm}.log"
    shell:
        """
        export LC_ALL=C
        RepeatMasker -dir {params.out_dir} \
        -lib {input.families} \
        -xsmall -gff -pa {params.pa} \
        {input.fa} \
        > {log} 2>&1
        """
