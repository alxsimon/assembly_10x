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
# InterProScan

rule clean_prot_seq:
    input:
        "results/annotation/braker/{asm}/augustus.hints.aa"
    output:
        "results/annotation/interproscan/{asm}/{asm}_augustus.hints.aa"
    conda: 
        "../envs/seqkit.yaml"
    shell:
        "seqkit replace -s -p '\*' -r '' {input} > {output}"

rule run_interproscan:
    input:
        "results/annotation/interproscan/{asm}/{asm}_augustus.hints.aa"
    output:
        multiext("results/annotation/interproscan/{asm}/interpro_{asm}",
            ".gff3", ".xml", ".json", ".tsv")
    params:
        ips_path = config['annotation']['interproscan_path'],
        prefix = lambda w, output: output[0].replace(".gff3", ""),
    log:
        "logs/annotation/interproscan_{asm}.log"
    threads:
        config['annotation']['interproscan_threads']
    shell:
        """
        {params.ips_path}/interproscan.sh \
        -b {params.prefix} \
        -f TSV,XML,JSON,GFF3 \
        -i {input} \
        --seqtype p \
        -goterms \
        -iprlookup \
        --pathways \
        --cpu {threads} \
        > {log} 2>&1
        """


#===========================================
# Comparative Annotation Toolkit

rule convert_Mcor_gff:
    input:
        'resources/Mco_ProteinCodingGenes.gff',
        'resources/GCA_017311375.1_Mcoruscus_HiC_assembly_report.txt'
    output:
        'resources/GCA017311375.gff'
    script:
        "../scripts/convert_mcor_gff.py"

# with python gff3tool 
# gff3_QC -g GCA....gff -f ....fa.corrected -o ....gff3_qc.txt -s ....gff3_qc.stat
# gff3_fix -qc_r GCA017311375.gff3_qc.txt -g GCA017311375.gff -og GCA017311375.corrected.gff


rule prepare_cat_config:
    input:
        mcor_gff = 'resources/GCA017311375.modif.gff',
        #mgal_gff = 'resources/GCA900618805.gff',
        bams_mgal = lambda w: get_sample_rna_runs_annotation(w, asm='gallo'),
        bams_mtro = lambda w: get_sample_rna_runs_annotation(w, asm='tros'),
        bams_medu = lambda w: get_sample_rna_runs_annotation(w, asm='edu'),
        prot_db = "resources/annotation/orthodb_mollusca_proteins.fa",
    output:
        "results/annotation/CAT/cat.config"
    run:
        bams_genomes = {
            'mgal_02': input['bams_mgal'],
            'mtro_02': input['bams_mtro'],
            'medu_02': input['bams_medu'],
        }
        with open(output[0], 'w') as fw:
            fw.write('[ANNOTATION]\n')
            fw.write(f"GCA017311375 = {input['mcor_gff']}\n")
            #fw.write(f"GCA900618805 = {input['mgal_gff']}\n\n")
            fw.write("[BAM]\n")
            for genome in ['mgal_02', 'mtro_02', 'medu_02']:
                fw.write(f"{genome} = {','.join(bams_genomes[genome])}\n")
            fw.write('\n')
            fw.write("[PROTEIN_FASTA]\n")
            for genome in ['mgal_02', 'mtro_02', 'medu_02']:
                fw.write(f"{genome} = {input['prot_db']}\n")


rule comparative_annotation:
    input:
        config = "results/annotation/CAT/cat.config",
        hal = "results/cactus/myt_cactus.hal",
    output:
        dir = directory("results/annotation/CAT/cat_annot")
    log:
        "logs/cat_pipeline.log"
    threads:
        10
    container:
        "containers/cat.sif"
    shell:
        """
        luigi --module cat RunCat \
        --binary-mode=local \
        --hal={input.hal} \
        --ref-genome='GCA017311375' \
        --maxCores=1 \
        --workers={threads} \
        --config={input.config} \
        --work-dir {output.dir} \
        --out-dir {output.dir} \
        --local-scheduler \
        --augustus --augustus-species 'Caenorhabditis_elegans' \
        --augustus-cgp \
        --assembly-hub \
        > {log} 2>&1
        """
