
# rule supernova_v1



#==================================================
rule supernova_v2:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_regen_{R}_001.fastq.gz", R=["R1", "R2"])
    output:
        "results/supernova_assemblies/{sample}_v2/outs/summary.txt"
    params:
        mem = config['supernova_mem'],
        input_dir = lambda w: f'results/preprocessing/{w.sample}',
        output_dir = lambda w: f'results/supernova_assemblies/{w.sample}_v2'
    threads: 
        workflow.cores
    log: 
        "logs/supernova_v2.{sample}.log"
    shell:
        """
        supernova run \
        --id {params.output_dir} \
        --fastqs {params.input_dir} \
        --sample {sample}_dedup_regen \
        --maxreads='all' \
        --localcores={threads} \
        --localmem={params.mem} \
        --allow-extreme-coverage \
        > {log} 2>&1
        """


#==================================================
def get_output_fasta(wildcards, just_prefix=True):
    prefix = "results/supernova_assemblies/{sample}_{version}/fasta/{sample}_{version}.{style}"
    if just_prefix:
        return(prefix)
    else:
        if wildcards.style is "pseudohap2":
            return([f'{prefix}.1.fasta.gz', f'{prefix}.2.fasta.gz'])
        else:
            return(f'{prefix}.1.fasta.gz')

rule supernova_fasta:
    input:
        "results/supernova_assemblies/{sample}_{version}/outs/summary.txt"
    output:
        get_output_fasta(wildcards, just_prefix=False)
    params:
        fasta_dir = lambda w, output: os.path.dirname(output[0]),
        asm_dir = lambda w, input: os.path.dirname(input) + "/assembly",
        outprefix = get_output_fasta(wildcards, just_prefix=True)
    log:
        "logs/supernova_fasta.{sample}_{version}.{style}.log"
    shell:
        """
        [ ! -d {params.fasta_dir} ] && mkdir {params.fasta_dir}
        supernova mkoutput \
        --style = {wildcards.style} \
        --asmdir = {params.assembly_dir} \
        --outprefix = {params.outprefix} \
        --headers = full \
        > {log} 2>&1
        """


rule supernova_compress:
    input:
        rule.supernova_fasta.output
    output:
        archive = "results/supernova_assemblies/{sample}_{version}/outs/assembly.tar.zst",
        donefile = "results/supernova_assemblies/{sample}_{version}/DONE"
    params:
        input_dir = lambda w: f'results/supernova_assemblies/{w.sample}_{w.version}/outs/assembly'
    threads: 
        8
    shell:
        """
        tar -cf - {params.input_dir} | zstdmt -T{threads} > {output.archive} \
        && touch {output.donefile}
        """