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
        samtools sort -@ {threads} -O BAM -o {output[0]}
        samtools index {output}
        """

def get_sample_rna_runs_annotation(w, asm=None):
    if asm is None:
        sample = w.asm.replace("_v7", "")
    else:
        sample = asm
    list_R1_files = glob.glob(f"resources/RNAseq_raw/{sample}/*_R1.fastq.gz")
    list_runs = [re.sub('_R1\.fastq\.gz$', '', os.path.basename(f)) for f in list_R1_files]
    return [f'results/annotation/RNAseq/{sample}_v7/{run}.bam' for run in list_runs]

species_dict = {
    'edu_v7': 'Medulis',
    'gallo_v7': 'Mgalloprovincialis',
    'tros_v7': 'Mtrossulus',
}

# To run this step, you need to have a GeneMark key available in your home folder.
# Refer to Braker2 installation steps.
rule braker:
    input:
        genome = "results/repeats/{asm}.fa.masked",
        rna_bams = get_sample_rna_runs_annotation,
        prot_db = "resources/annotation/orthodb_mollusca_proteins.fa",
    output:
        "results/annotation/braker/{asm}/braker.gtf",
    params:
        out_dir = lambda w: f"results/annotation/braker/{w.asm}",
        species = lambda w: species_dict[w.asm],
        list_bams = lambda w, input: ','.join(input['rna_bams']),
        genemark_path = config['annotation']['genemark_path'],
        to_rm = 'GeneMark*',
    conda:
        "../envs/annotation_braker.yaml"
    threads:
        config['annotation']['braker_threads']
    log:
        "logs/annotation/braker2_{asm}.log"
    shell:
        """
        find {params.out_dir} -name '{params.to_rm}' -type d -delete
        cp /opt/gm_key_64 ~/.gm_key
        braker.pl \
        --genome {input.genome} \
        --prot_seq {input.prot_db} \
        --bam {params.list_bams} \
        --workingdir {params.out_dir} \
        --species {params.species} \
        --useexisting \
        --etpmode --softmasking --cores {threads} \
        --gff3 \
        --GENEMARK_PATH {params.genemark_path} \
        --PROTHINT_PATH {params.genemark_path}/ProtHint/bin \
        > {log} 2>&1
        """


#===========================================
# Mantis annotation

rule download_mantis:
    output:
        directory("results/annotation/mantis/mantis")
    shell:
        """
        git clone https://github.com/PedroMTQ/mantis.git {output}
        cd {output}
        git checkout a59937d
        """

org_detail = {
    'MgalMED': 'Mytilus galloprovincialis',
    'MeduEUS': 'Mytilus edulis',
    'MeduEUN': 'Mytilus edulis',
}

rule get_prot_seq:
    input:
        gff = "results/final/{asm}/{asm}.gff3.gz",
        fa = "results/final/{asm}/{asm}.fa.gz",
    output:
        "results/annotation/mantis/{asm}_pep.fa"
    params:
        prefix = lambda w, output: output[0].replace('_pep.fa', '')
    log:
        "logs/annotation/get_prot_seq.{asm}.log"
    conda:
        "../envs/gff3tool.yaml"
    shell:
        """
        gff3_to_fasta -g {input.gff} -f {input.fa} \
        -st pep -d simple -o {params.prefix} \
        > log 2>&1
        """

rule run_mantis:
    input:
        config = "config/mantis.config",
        mantis_dir = "results/annotation/mantis/mantis",
        pep = "results/annotation/mantis/{asm}_pep.fa",
    output:
        "results/annotation/mantis/{asm}/consensus_annotation.tsv"
    params:
        outdir = lambda w, output: output[0].replace('/consensus_annotation.tsv', ''),
        org_detail = lambda w: org_detail[w.asm],
    log:
        "logs/annotation/mantis.{asm}.log"
    conda:
        "../envs/mantis_env.yaml"
    threads:
        32
    shell:
        """
        rm -r {params.outdir}
        python {input.mantis_dir} run_mantis \
        -t {input.pep} \
        -o {params.outdir} \
        -mc {input.config} \
        -od '{params.org_detail}' \
        -km \
        -gff \
        -c {threads} \
        > {log} 2>&1
        """
        # remove folder created by Snakemake first, overwise Mantis creates another folder with date-time

