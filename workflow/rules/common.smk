
def get_fastq(wildcards):
    prefix = config["raw_names"][wildcards.sample]
    fq1 = glob.glob(f"resources/10x_reads/{prefix}*_R1*.fastq.gz")
    fq2 = glob.glob(f"resources/10x_reads/{prefix}*_R2*.fastq.gz")
    return {'fq1': fq1, 'fq2': fq2}

def proc10x_expand(prefix):
    fq1 = f'{prefix}_R1_001.fastq.gz'
    fq2 = f'{prefix}_R2_001.fastq.gz'
    return({'fq1': fq1, 'fq2': fq2})

def rmdup_expand(prefix):
    fq1 = f'{prefix}_R1.fastq'
    fq2 = f'{prefix}_R2.fastq'
    return({'fq1': fq1, 'fq2': fq2})