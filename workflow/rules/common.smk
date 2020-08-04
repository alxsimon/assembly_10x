
def get_fastq(wildcards):
    prefix = config["raw_names"][wildcards.sample]
    fq1 = glob.glob(f"resources/10x_reads/{prefix}*R1*.fastq.gz")[0]
    fq2 = glob.glob(f"resources/10x_reads/{prefix}*R2*.fastq.gz")[0]
    return {'fq1': fq1, 'fq2': fq2}
