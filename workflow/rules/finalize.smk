
def get_final_fasta(w):
    alt = {
        'MeduEUS': 'edu_v7',
        'MeduEUN': 'tros_v7',
        'MgalMED': 'gallo_v7',
    }
    fasta = f"results/fasta/{alt[w.asm]}.pseudohap.fasta.gz"
    return fasta

def get_gff(w):
    alt = {
        'MeduEUS': 'edu_v7',
        'MeduEUN': 'tros_v7',
        'MgalMED': 'gallo_v7',
    }
    gff = f"results/annotation/braker/{alt[w.asm]}/braker.gff3"
    return gff

rule cp_final_fasta:
    input:
        get_final_fasta
    output:
        "results/final/{asm}.fa.gz"
    shell:
        """
        cp {input} {output}
        """

rule clean_gff_id:
    input:
        get_gff
    output:
        temp("results/final/{asm}_newid.gff3")
    params:
        prefix = lambda w: w.asm
    conda:
        "../envs/agat.yaml"
    log:
        "logs/finalize/clean_gff_id.{asm}.log"
    shell:
        """
        agat_sp_manage_IDs.pl --gff {input} \
        --ensembl --prefix {params.prefix} \
        -o {output} > {log} 2>&1
        """

rule clean_gff_format:
    input:
        "results/final/{asm}_newid.gff3"
    output:
        temp("results/final/{asm}_fix.gff3")
    conda:
        "../envs/agat.yaml"
    log:
        "logs/finalize/clean_gff_format.{asm}.log"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input} -gvi 3 -gvo 3 \
        -o {output} > {log} 2>&1
        """

rule sort_bgzip_gff:
    input:
        "results/final/{asm}_fix.gff3"
    output:
        "results/final/{asm}.gff3.gz"
    params:
        tmp = lambda w, output: output[0].replace(".gz", "")
    conda:
        "../envs/agat.yaml"
    log:
        "logs/finalize/sort_bgzip_gff.{asm}.log"
    shell:
        """
        gff3sort.pl {input} > {params.tmp}
        bgzip {params.tmp}
        tabix -p gff {output}
        """

rule cp_pep_seq:
    input:
        "results/annotation/mantis/{asm}_pep.fa"
    output:
        "results/final/{asm}_pep.fa.gz"
    shell:
        "cat {input} | gzip -c > {output}"

rule cp_consensus_annotation:
    input:
        "results/annotation/mantis/{asm}/consensus_annotation.tsv"
    output:
        "results/final/{asm}_consensus_annotation.tsv"
    shell:
        "cp {input} {output}"