
# rule supernova_v1


rule supernova_v2:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_dedup_regen_{R}_001.fastq.gz", R=["R1", "R2"])
    output:
        "results/supernova_assemblies/{sample}_v2/outs/summary.txt"
    params:
        mem = config['supernova_mem']
    threads: 
        workflow.cores
    log: 
        "logs/supernova_v2.{sample}.log"
    shell:
        """
        supernova run \
        ...
        """


rule supernova_fasta:
    {type} {sample}


rule supernova_compress:
    input:

    output:

    threads: 
        8
    shell:
        """
        tar -cf - {input folder} | zstdmt -T{threads} > {output}
        """