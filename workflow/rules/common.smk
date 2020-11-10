
def get_fastq(wildcards):
    prefix = config["raw_names"][wildcards.sample]
    fq1 = glob.glob(f"resources/10x_reads/{prefix}*R1*.fastq.gz")[0]
    fq2 = glob.glob(f"resources/10x_reads/{prefix}*R2*.fastq.gz")[0]
    return {'fq1': fq1, 'fq2': fq2}

def get_supernova_input(w):
    if os.path.exists(f'results/supernova_assemblies/{w.sample}_{w.version}/outs/report.txt'):
        return ['']
    elif w.version == "v1":
        return get_fastq(w)
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

rule download_assemblies:
    output:
        "resources/GCA001676915.fasta.gz",
        "resources/UYJE01.fasta.gz"
    params:
        ftp = config['published_assemblies']
    shell:
        """
        wget -O {output[0]} {params.ftp[0]}
        wget -O {output[1]} {params.ftp[1]}
        """