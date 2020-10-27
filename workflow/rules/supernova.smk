rule supernova_assembly:
    input:
        unpack(get_supernova_input),
#	    get_order
    output:
        temp("tmp/{sample}_{version}/outs/report.txt")
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
        [ ! -e {params.run_id}/{params.run_id}.mri.tgz ] && rm -r {params.run_id}
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
        protected("results/fasta/{sample}_{version}.{style}.fasta.gz")
    params:
        style = lambda w: w.style.strip('.1'),
        fasta_dir = lambda w, output: os.path.dirname(output[0]),
        asm_dir = lambda w, input: os.path.dirname(input[0]) + "/assembly",
        outprefix = lambda w, output: os.path.dirname(output[0]) + f"/{w.sample}_{w.version}.{w.style.strip('.1')}"
    wildcard_constraints:
        style = '\w+|\w+.1',
        version = 'v1|v2'
    log:
        "logs/supernova_fasta.{sample}_{version}.{style}.log"
    container:
        "containers/supernova.sif"
    shell:
        """
        supernova mkoutput \
        --style={params.style} \
        --asmdir={params.asm_dir} \
        --outprefix={params.outprefix} \
        --headers=full \
        > {log} 2>&1
        """


rule supernova_compress_move:
    input:
        multiext("results/fasta/{sample}_{version}",
            ".raw.fasta.gz", ".megabubbles.fasta.gz", ".pseudohap.fasta.gz",
            ".pseudohap2.1.fasta.gz"),
        "tmp/{sample}_{version}/outs/report.txt"
    output:
        archive = "results/supernova_assemblies/{sample}_{version}/outs/assembly.tar.zst",
        donefile = touch("results/supernova_assemblies/{sample}_{version}/DONE")
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
        tar -cf - -C {params.input_dir}/outs assembly | zstdmt -T{threads} > {params.tmp_archive} 2> {log} && \
        rm -r {params.input_dir}/outs/assembly && \
        cp -r {params.input_dir} {params.output_dir} && rm -r {params.input_dir}
        """

rule collect_stats:
    input:
        expand("results/supernova_assemblies/{sample}_{version}/outs/summary.csv",
            sample=config['samples'], version=["v1", "v2"])
    output:
        "results/supernova_assemblies/supernova_assemblies_stats.csv"
    shell:
        "bash workflow/scripts/collect_assembly_stats.sh"