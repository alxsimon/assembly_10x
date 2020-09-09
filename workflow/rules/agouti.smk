checkpoint split_fa_augustus:
    input:
        "results/fasta/{sample}_v4.pseudohap.fasta.gz"
    output:
        "results/agouti/{sample}/{sample}_v4.pseudohap.fa",
        directory("results/agouti/{sample}/split")
    params:
        split_size = 50000000
    conda:
        "../envs/augustus.yaml"
    shell:
        """
        zcat {input} > {output[0]}
        mkdir {output[1]}
        splitMfasta.pl {output[0]} \
        --outputpath={output[1]} --minsize={params.split_size}
        """

rule augustus:
    input:
        "results/agouti/{sample}/split/{sample}_v4.pseudohap.split.{i}.fa"
    output:
        "results/agouti/{sample}/split/pred_{i}.gff3"
    conda:
        "../envs/augustus.yaml"
    shell:
        """
        augustus --gff3=on --species=caenorhabditis {input} > {output}
        """

def aggregate_input_gff3(wildcards):
    checkpoint_output = checkpoints.split_fa_augustus.get(**wildcards).output[1]
    return expand("results/agouti/{sample}/split/pred_{i}.gff3",
           sample=wildcards.sample,
           i=glob_wildcards(os.path.join(checkpoint_output, f"{wildcards.sample}_v4.pseudohap.split." + "{i}.fa")).i)

rule aggregate_gff3:
    input:
        aggregate_input_gff3
    output:
        "results/agouti/{sample}/{sample}_v4.pseudohap.gff3"
    conda:
        "../envs/augustus.yaml"
    shell:
        "cat {input} | join_aug_pred.pl > {output}"

#rule index_ref:
#
#
#rule map_RNAseq:
#
#    conda:
#        "../envs/mapping.yaml"
#
#
#rule merge_RNA_bams:
#
#
#rule agouti_scaff:
#    input: 
#        fasta
#        bam
#        gff
#    output: 
#        ??
#    params:
#        outdir = lambda ...,
#        minMQ = 20,
#        maxFracMM = 0.05
#    conda: 
#        "../envs/py2.yaml"
#    shell:
#        """
#        python /opt/agouti/agouti.py scaffold \
#        -assembly {input.fa} \
#        -bam {input.bam} \
#        -gff {input.gff} \
#        -outdir {params.outdir} \
#        -minMQ {params.minMQ} -maxFracMM {params.maxFracMM}
#        """