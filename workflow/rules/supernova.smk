
#==================================================
# Version 1 is with raw input fastq

rule supernova_v1:
    input:
        unpack(get_fastq)
    output:
        "results/supernova_assemblies/{sample}_v1/outs/report.txt"
    params:
        mem = config['supernova_mem'],
        input_dir = lambda w, input: os.path.dirname(input[0]),
        output_dir = "results/supernova_assemblies",
        run_id = lambda w: f'{w.sample}_v1',
        sample = lambda w: config['raw_names'][w.sample]
    threads: 
        workflow.cores
    log: 
        "logs/supernova_v1.{sample}.log"
    container:
        "containers/supernova.sif"
    shell:
        """
        cd tmp
        supernova run \
        --id={params.run_id} \
        --fastqs=../{params.input_dir} \
        --sample={params.sample} \
        --maxreads='all' \
        --localcores={threads} \
        --localmem={params.mem} \
        --accept-extreme-coverage \
        > {log} 2>&1;
        cp -r {params.run_id} ../{params.output_dir}/ && \
        rm -r {params.run_id}
        """


#==================================================
# Version 2 is with filtered input fastq (preprocessing)
rule supernova_v2:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_regen_{R}_001.fastq.gz", R=["R1", "R2"])
    output:
        "results/supernova_assemblies/{sample}_v2/outs/report.txt"
    params:
        mem = config['supernova_mem'],
        input_dir = lambda w, input: os.path.dirname(input[0]),
        output_dir = "results/supernova_assemblies",
        run_id = lambda w: f'{w.sample}_v2',
        sample = lambda w: f'{w.sample}_dedup_regen'
    threads: 
        workflow.cores
    log: 
        "logs/supernova_v2.{sample}.log"
    container:
        "containers/supernova.sif"
    shell:
        """
        cd tmp
        supernova run \
        --id={params.run_id} \
        --fastqs=../{params.input_dir} \
        --sample={params.sample} \
        --maxreads='all' \
        --localcores={threads} \
        --localmem={params.mem} \
        --accept-extreme-coverage \
        > {log} 2>&1;
        cp -r {params.run_id} ../{params.output_dir}/ && \
        rm -r {params.run_id}
        """


#==================================================
# common for both versions

rule supernova_fasta:
    input:
        "results/supernova_assemblies/{sample}_{version}/outs/report.txt"
    output:
        multiext("results/supernova_assemblies/{sample}_{version}/fasta/{sample}_{version}",
            ".raw.fasta.gz", ".megabubbles.fasta.gz", ".pseudohap.fasta.gz",
            ".pseudohap2.1.fasta.gz", ".pseudohap2.2.fasta.gz")
    params:
        fasta_dir = lambda w, output: os.path.dirname(output[0]),
        asm_dir = lambda w, input: os.path.dirname(input[0]) + "/assembly",
        outprefix = lambda w, output: os.path.dirname(output[0]) + f"/{w.sample}_{w.version}"
    log:
        "logs/supernova_fasta.{sample}_{version}"
    container:
        "containers/supernova.sif"
    shell:
        """
        [ ! -d {params.fasta_dir} ] && mkdir {params.fasta_dir};
        for style in ( 'raw' 'megabubbles' 'pseudohap' pseudohap2 ); do
            supernova mkoutput \
            --style = $style \
            --asmdir = {params.asm_dir} \
            --outprefix = "{params.outprefix}.$style" \
            --headers = full \
            > "{log}.$style.log" 2>&1
        done
        """


rule supernova_compress:
    input:
        multiext("results/supernova_assemblies/{sample}_{version}/fasta/{sample}_{version}",
            ".raw.fasta.gz", ".megabubbles.fasta.gz", ".pseudohap.fasta.gz",
            ".pseudohap2.1.fasta.gz", ".pseudohap2.2.fasta.gz")
    output:
        archive = "results/supernova_assemblies/{sample}_{version}/outs/assembly.tar.zst",
        donefile = "results/supernova_assemblies/{sample}_{version}/DONE"
    params:
        input_dir = lambda w: f'results/supernova_assemblies/{w.sample}_{w.version}/outs/assembly'
    threads: 
        8
    log:
        "logs/supernova_compress.{sample}_{version}.log"
    shell:
        """
        tar -cf - {params.input_dir} | zstdmt -T{threads} > {output.archive} 2> {log} \
        && touch {output.donefile}
        """