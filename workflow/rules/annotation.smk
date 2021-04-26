# Annotation pipeline

rule download_orthodb_mollusca:
    output:
        clustids = "resources/annotation/orthodb_mollusca_proteins.clustids",
        fasta = "resources/annotation/orthodb_mollusca_proteins.fa",
    run:
        import subprocess
        import json
        if os.path.exists(output.fasta):
            os.remove(output.fasta)
        # Get clusters for Mollusca taxid=6447
        url = 'http://www.orthodb.org/search?level=6447&limit=50000'
        cmd = f'wget "{url}" -O {output.clustids}'
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
        # Loop for each cluster
        with open(output.clustids, 'r') as fr:
            clusters = json.load(fr)
        for C_id in clusters['data']:
            url = f'http://www.orthodb.org/fasta?id={C_id}&species=all'
            cmd = f'wget "{url}" -O - >> {output.fasta}'
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)

#=================================
# We use here preprocessed RNAseq data from the AGOUTI step of the assembly
# see workflow/rules/agouti.smk
# preprocessing: rcorrector + trim_galore
#=================================

rule hisat2_index_ref:
    input:
        "results/fasta/{asm}.pseudohap.fasta.gz"
    output:
        expand("results/annotation/RNAseq/hisat2_index_{{asm}}/{{asm}}.{n}.ht2",
            n=range(1,9))
    params:
        index_prefix = lambda w, output: output[0].replace('.1.ht2', '') 
    conda:
        "../envs/annotation_hisat2.yaml"
    log:
        "logs/annotation/hisat2_index_{asm}.log"
    shell:
        """
        zcat {input} > tmp_ref.{wildcards.asm}.fa
        hisat2-build tmp_ref.{wildcards.asm}.fa {params.index_prefix} \
        > {log} 2>&1
        rm tmp_ref.{wildcards.asm}.fa
        """

def get_fastq_rna(w):
    sample = w.asm.replace("_v7", "")
    fastq = expand("results/RNA_preproc/{sample}/{{run}}_trimgal_val_{i}.fq",
            i=['1', '2'], sample=sample)
    return fastq

rule hisat2_map:
    input: 
        reads = get_fastq_rna,
        index = "results/annotation/RNAseq/hisat2_index_{asm}/{asm}.1.ht2",
    output:
        "results/annotation/RNAseq/{asm}/{run}.bam",
        "results/annotation/RNAseq/{asm}/{run}.bam.bai",
    params:
        index_prefix = lambda w, input: input['index'].replace('.1.ht2', '') 
    log:
        "logs/annotation/hisat2_map_rna_{asm}_{run}.log"
    conda:
        "../envs/annotation_hisat2.yaml"
    threads:
        config['annotation']['hisat2_threads']
    shell:
        """
        hisat2 -p {threads} \
        -x {params.index_prefix} \
        -1 {input.reads[0]} -2 {input.reads[1]} \
        2> {log} | \
        samtools view -b -@ {threads} | \
        samtools sort -@ {threads} - > {output}
        samtools index {output}
        """

def get_sample_rna_runs_annotation(w):
    sample = w.asm.replace("_v7", "")
    list_R1_files = glob.glob(f"resources/RNAseq_raw/{sample}/*_R1.fastq.gz")
    list_runs = [re.sub('_R1\.fastq\.gz$', '', os.path.basename(f)) for f in list_R1_files]
    return [f'results/annotation/{w.asm}/{run}.bam' for run in list_runs]

rule braker:
    input:
        genome = "results/repeats/{asm}.fa.masked",
        rna_bam = get_sample_rna_runs_annotation,
        prot_db = "resources/annotation/orthodb_mollusca_proteins.fa",
    output:
        "dummy"
    conda:
        "../envs/annotation_braker.yaml"
    threads:
        config['annotation']['braker_threads']
    shell:
        """
        braker.pl --genome {input.genome} \
        --prot_seq {input.prot_db} \
        --bam {input.rna_bam} \
        --etpmode --softmasking --cores {threads} \
        > {log} 2>&1
        """
        # --bam option all bams separated by commas