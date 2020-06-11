
def get_fastq(wildcards):
    prefix = config["raw_names"][wildcards.sample]
    fq1 = glob.glob(f"resources/10x_reads/{prefix}*_R1*.fastq.gz")
    fq2 = glob.glob(f"resources/10x_reads/{prefix}*_R2*.fastq.gz")
    return {'fq1': fq1, 'fq2': fq2}
