# Scaffold according to Mytilus coruscus genome
# Finally perform correction and gap filling with Pilon

def get_fasta_to_improve(w):
    if w.asm=="mgal":
        return ancient("results/fasta/gallo_v6.pseudohap.fasta.gz")
    elif w.asm=="medu":
        return ancient("results/fasta/edu_v6.pseudohap.fasta.gz")

rule cp_asm_fasta:
    input:
        get_fasta_to_improve
    output:
        "results/{asm}_02/{asm}_01.fa",
    conda:
        "../envs/asm_improvement.yaml"
    shell:
        """
        zcat {input} > {output}
        """

# ragtag
rule ragtag:
    input:
        ref_coruscus = "resources/GCA017311375.fasta.gz",
        draft_assembly = "results/{asm}_02/{asm}_01.fa"
    output:
        "results/{asm}_02/ragtag/{asm}_01.ragtag.fa"
    params:
        out_dir = lambda w, output: f"{os.path.dirname(output[0])}/res"
    conda:
        "../envs/asm_improvement.yaml"
    threads:
        workflow.cores
    log:
        "logs/{asm}_improvement/ragtag_{asm}.log"
    shell:
        """
        ragtag.py scaffold -t {threads} \
        {input.ref_coruscus} \
        {input.draft_assembly} \
        -o {params.out_dir} \
        --mm2-params '-x asm10' -u \
        > {log} 2>&1

        cp {params.out_dir}/ragtag.scaffolds.fasta {output}
        """

# Pilon

rule prepare_bwa_index:
    input:
        "results/{asm}_02/ragtag/{asm}_01.ragtag.fa"
    output:
        "results/{asm}_02/ragtag/{asm}_01.ragtag.fa.0123"
    conda:
        "../envs/asm_improvement.yaml"
    shell:
        """
        bwa-mem2 index {input}
        """

def get_fastq_pilon_map(w):
    if w.asm=="mgal":
        R1 = "results/preprocessing/gallo/gallo_dedup_proc_fastp_R1_001.fastq.gz"
        R2 = "results/preprocessing/gallo/gallo_dedup_proc_fastp_R2_001.fastq.gz"
    elif w.asm=="medu":
        R1 = "results/preprocessing/edu/edu_dedup_proc_fastp_R1_001.fastq.gz"
        R2 = "results/preprocessing/edu/edu_dedup_proc_fastp_R2_001.fastq.gz"
    return {'R1': R1, 'R2': R2}

rule pilon_map_pe:
    input:
        unpack(get_fastq_pilon_map),
        ref = "results/{asm}_02/ragtag/{asm}_01.ragtag.fa",
        index = "results/{asm}_02/ragtag/{asm}_01.ragtag.fa.0123",
    output:
        "results/{asm}_02/pilon/{asm}_01.ragtag.pe_mapped_pilon.bam"
    log:
        "logs/{asm}_improvement/map_pe_reads_pilon.log"
    conda:
        "../envs/asm_improvement.yaml"
    threads:
        workflow.cores
    shell:
        """
        [ ! -e {input.ref}.0123 ] && bwa-mem2 index {input.ref}
        bwa-mem2 mem -t {threads} {input.ref} {input.R1} {input.R2} 2> {log} | \
        samtools sort -@ {threads} | \
        samtools view -b -@ {threads} -o {output}
        samtools index {output}
        """

rule download_pilon:
    output:
        "results/{asm}_02/pilon/pilon-1.24.jar"
    params:
        link = "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar"
    shell:
        "wget -O {output} {params.link}"


rule split_in_targets:
    input:
        ref = "results/{asm}_02/ragtag/{asm}_01.ragtag.fa"
    output:
        expand("results/{asm}_02/pilon/target_{num}.txt",
            num=[f'{i+1:02}' for i in range(15)])
    params:
        out_dir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/asm_improvement.yaml"
    shell:
        """
        for i in {{1..14}}; do
            grep '>' {input.ref} | sed 's/>//' | awk -v record=$i 'NR==record {{print $0}}' \
            > {params.out_dir}/target_$(printf '%02d' $i).txt
        done
        grep '>' {input.ref} | tail -n +15 > {output[14]}
        """

rule pilon_asm_02:
    input:
        target = "results/{asm}_02/pilon/target_{num}.txt",
        ref = "results/{asm}_02/ragtag/{asm}_01.ragtag.fa",
        pilon = "results/{asm}_02/pilon/pilon-1.24.jar",
        pe_bam = "results/{asm}_02/pilon/{asm}_01.ragtag.pe_mapped_pilon.bam",
    output:
        "results/{asm}_02/pilon/res_{asm}_02_{num}/{asm}_02_{num}.fasta"
    params:
        out_dir = lambda w, output: os.path.dirname(output[0]),
        prefix = lambda w, output: os.path.basename(output[0]).replace('.fasta', ''),
        java_mem = config['asm_improvement']['java_mem'],
    log:
        "logs/{asm}_improvement/pilon_run_target_{num}.log"
    conda:
        "../envs/asm_improvement.yaml"
    threads: 
        math.floor(workflow.cores/config['asm_improvement']['max_concurrent_pilon'])
        # limit number of concurrent processes due to available RAM
    shell:
        """
        java -Xmx{params.java_mem}G -jar {input.pilon} \
        --genome {input.ref} \
        --frags {input.pe_bam} \
        --targets {input.target} \
        --output {params.prefix} \
        --outdir {params.out_dir} \
        --changes --vcf --tracks --diploid --fix all \
        > {log} 2>&1
        """

def decide_asm(w):
    if w.sample=="gallo":
        asm = "mgal"
    elif w.sample=="edu":
        asm = "medu"
    return expand(f"results/{asm}_02/pilon/res_{asm}_02_{num}/{asm}_02_{num}.fasta",
        num=[f'{i+1:02}' for i in range(15)])

rule merge_pilon_res:
    input:
        decide_asm
    output:
        "results/fasta/{sample}_v7.pseudohap.fasta.gz"
    wildcard_constraints:
        sample = '!tros'
    shell:
        "cat {input} | bgzip -c > {output}"