def get_supernova_input(w):
    if w.version == "v1":
        return unpack(get_fastq(w))
    else:
        return expand("results/preprocessing/{{sample}}/{{sample}}_regen_{R}_001.fastq.gz", R=["R1", "R2"])

def get_order(w):
    supernova_order = config['supernova_order']
    current_assembly = f'{w.sample}_{w.version}'
    if current_assembly != supernova_order[0]:
        previous_assembly = supernova_order[supernova_order.index(current_assembly) - 1]
        return f'results/supernova_assemblies/{previous_assembly}/DONE'
    else:
        # dummy file that already exist
        return 'workflow/rules/supernova.smk'

rule supernova_assembly:
    input:
        get_supernova_input,
	get_order
    output:
        "results/supernova_assemblies/{sample}_{version}/outs/report.txt"
    params:
        mem = config['supernova_mem'],
        input_dir = lambda w, input: os.path.dirname(input[0]),
        output_dir = "results/supernova_assemblies",
        run_id = lambda w: f'{w.sample}_{w.version}',
        sample = lambda w: config['raw_names'][w.sample]
    threads: 
        workflow.cores
    log: 
        "logs/supernova_{version}.{sample}.log"
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
        > ../{log} 2>&1;
        cp -r {params.run_id} ../{params.output_dir}/ && \
        rm -r {params.run_id}
        """

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
        workflow.cores
    log:
        "logs/supernova_compress.{sample}_{version}.log"
    shell:
        """
        tar -cf - {params.input_dir} | zstdmt -T{threads} > {output.archive} 2> {log} \
        && touch {output.donefile}
        """
