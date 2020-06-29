def get_supernova_input(w):
    if w.version == "v1":
        return unpack(get_fastq(w))
    else:
        return expand("results/preprocessing/{{sample}}/{{sample}}_S1_L001_{R}_001.fastq.gz", R=["R1", "R2"])

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
        "tmp/{sample}_{version}/outs/report.txt"
    params:
        mem = config['supernova_mem'],
        input_dir = lambda w, input: os.path.dirname(input[0]),
        run_id = lambda w: f'{w.sample}_{w.version}',
        sample = lambda w, input: re.sub("_S.+_L.+_R1_001.fastq.gz", "", os.path.basename(input[0]))
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
        """

rule supernova_fasta:
    input:
        "tmp/{sample}_{version}/outs/report.txt"
    output:
        "tmp/{sample}_{version}/fasta/{sample}_{version}.{style}.fasta.gz"
    params:
        style = lambda w: w.style.strip('.1'),
        fasta_dir = lambda w, output: os.path.dirname(output[0]),
        asm_dir = lambda w, input: os.path.dirname(input[0]) + "/assembly",
        outprefix = lambda w, output: os.path.dirname(output[0]) + f"/{w.sample}_{w.version}.{w.style.strip('.1')}"
    wildcard_constraints:
        style = '\w+|\w+.1'
    log:
        "logs/supernova_fasta.{sample}_{version}.{style}.log"
    container:
        "containers/supernova.sif"
    shell:
        """
        [ ! -d {params.fasta_dir} ] && mkdir {params.fasta_dir};
        supernova mkoutput \
        --style = {params.style} \
        --asmdir = {params.asm_dir} \
        --outprefix = {params.outprefix} \
        --headers = full \
        > {log} 2>&1
        """


rule supernova_compress_move:
    input:
        multiext("tmp/{sample}_{version}/fasta/{sample}_{version}",
            ".raw.fasta.gz", ".megabubbles.fasta.gz", ".pseudohap.fasta.gz",
            ".pseudohap2.1.fasta.gz")
    output:
        archive = "results/supernova_assemblies/{sample}_{version}/outs/assembly.tar.zst",
        donefile = "results/supernova_assemblies/{sample}_{version}/DONE"
    params:
        input_dir = lambda w: f'tmp/{w.sample}_{w.version}',
        tmp_archive = lambda w: f'tmp/{w.sample}_{w.version}/outs/assembly.tar.zst',
        output_dir = 'results/supernova_assemblies/'
    threads: 
        workflow.cores
    log:
        "logs/supernova_compress.{sample}_{version}.log"
    shell:
        """
        tar -cf - {params.input_dir}/outs/assembly | zstdmt -T{threads} > {params.tmp_archive} 2> {log}
        cp -r {params.input_dir} {params.output_dir} && rm -r {params.input_dir}
        touch {output.donefile}
        """
