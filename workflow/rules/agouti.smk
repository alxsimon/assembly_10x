# Prepare gff3 for genomes
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

#===============================
# Preprocess RNAseq data

# The RNAseq reads need to be in the folder resources/RNAseq_raw/{sample}
# Files must be named {run}_R1.fastq.gz and {run}_R2.fastq.gz for globbing to work
# globbing is done in the rule merge_RNA_bams

rule rna_rcorrector:
    input:
        expand("resources/RNAseq_raw/{{sample}}/{{run}}_{R}.fastq.gz",
            R=['R1', 'R2'])
    output:
        temp(expand("results/agouti/{{sample}}/RNA_preproc/{{run}}_{R}.cor.fq",
            R=['R1', 'R2']))
    params:
        outdir = lambda w, output: os.path.dirname(output[0])
    log:
        "logs/rcorrector_{sample}_{run}.log" 
    threads:
        config['agouti']['threads']
    conda:
        "../envs/rna_seq.yaml"
    shell:
        """
        perl rcorrector.pl \
        -1 {input[0]} -2 {input[1]} \
        -k 23 -t {threads} -maxcorK 4 -wk 0.95 -ek 10000000 \
        -od {params.outdir} \
        2>&1 {log}
        """

rule rna_trimgalore:
    input:
        expand("results/agouti/{{sample}}/RNA_preproc/{{run}}_{R}.cor.fq",
            R=['R1', 'R2'])
    output:
        temp(expand("results/agouti/{{sample}}/RNA_preproc/{{run}}_trimgal_val_{i}.fq",
            i=['1', '2']))
    params:
        outdir = lambda w, output: os.path.dirname(output[0]),
        basename = lambda w: f'{w.run}_trimgal'
    log:
        "logs/trimgalore_{sample}_{run}.log" 
    threads:
        config['agouti']['threads']
    conda:
        "../envs/rna_seq.yaml"
    shell: 
        """
        trim_galore --cores {threads} \
        --phred33 \
        --quality 20 \
        --stringency 1 \
        -e 0.1 \
        --length 70 \
        --output_dir {params.outdir} \
        --basename {params.basename} \
        --dont_gzip \
        --paired \
        {input} \
        2>&1 {log}
        """

#===============================
# Map the RNAseq reads
rule index_ref:
    input: 
        "results/agouti/{sample}/{sample}_v4.pseudohap.fa"
    output:
        multiext("results/agouti/{sample}/{sample}_v4.pseudohap.fa",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".bwt.8bit.32", ".pac")
    conda:
        "../envs/mapping.yaml"
    shell:
        "bwa-mem2 index {input}"

rule map_RNAseq:
    input: 
        expand("results/agouti/{{sample}}/RNA_preproc/{{run}}_trimgal_val_{i}.fq",
            i=['1', '2']),
        "results/agouti/{sample}/{sample}_v4.pseudohap.fa",
        multiext("results/agouti/{sample}/{sample}_v4.pseudohap.fa",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".bwt.8bit.32", ".pac")
    output:
        "results/agouti/{sample}/mapping/{run}.bam"
    log:
        "logs/bwa_rna_{sample}_{run}.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        config['agouti']['threads']
    shell:
        "bwa-mem2 mem -t {threads} {input[2]} {input[0]} {input[1]} "
        "| samtools view -b -@ {threads} -o {output} "
        "2>&1 {log}"

def get_sample_rna_runs(w):
    list_R1_files = glob.glob(f"resources/RNAseq_raw/{w.sample}/*_R1.fastq.gz")
    list_runs = [re.sub('*/', '', re.sub('_R1\.fastq\.gz$', '', f)) for f in list_R1_files]
    return [f'results/agouti/{w.sample}/mapping/{run}.bam' for run in list_runs]

rule merge_RNA_bams:
    input:
        get_sample_rna_runs
    output:
        "results/agouti/{sample}/RNAseq_mapped_merged.bam"
    params:
        tmp_merge = lambda w: f'results/agouti/{w.sample}/tmp_merge.bam'
    conda:
        "../envs/mapping.yaml"
    threads:
        config['agouti']['threads']
    shell:
        """
        samtools merge -@ {threads} {params.tmp_merge} {input}               
        samtools sort -@ {threads} -n -o {output} {params.tmp_merge}
        rm {params.tmp_merge}
        """

#===============================
# Run agouti on all that
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