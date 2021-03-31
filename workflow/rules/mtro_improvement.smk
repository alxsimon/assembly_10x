# Use oxford nanopore reads to slightly improve scaffolding
# Then scaffold according to Mytilus coruscus genome
# Finally perform correction and gap filling with Pilon

# rule cp_mtro_fasta:
#     input:
#         ancient("results/fasta/tros_v6.pseudohap.fasta.gz")
#     output:
#         "results/mtro_02/mtro_01.fa",
#     conda:
#         "../envs/mtro_improvement.yaml"
#     shell:
#         """
#         zcat {input} > {output}
#         """

rule download_lrscaf:
    output:
        "results/mtro_02/lrscaf/LRScaf-1.1.10.jar"
    params:
        link = "https://github.com/shingocat/lrscaf/releases/download/1.1.10/LRScaf-1.1.10.jar"
    shell:
        "wget -O {output} {params.link}"

rule map_ont_reads:
    input:
        ref = "results/mtro_02/mtro_01.fa",
        ont_reads = config['mtro_improvement']['mtro_ont'],
    output:
        "results/mtro_02/lrscaf/mtro_01_ont_mapped.paf"
    log:
        "logs/mtro_improvement/map_ont_reads.log"
    conda:
        "../envs/mtro_improvement.yaml"
    threads:
        workflow.cores
    shell:
        """
        minimap2 -t {threads} -x map-ont \
        {input.ref} {input.ont_reads} \
        > {output} 2> {log}
        """

rule lrscaf:
    input:
        lrscaf = "results/mtro_02/lrscaf/LRScaf-1.1.10.jar",
        paf = "results/mtro_02/lrscaf/mtro_01_ont_mapped.paf",
        ref = "results/mtro_02/mtro_01.fa",
    output:
        "results/mtro_02/lrscaf/mtro_01.lrscaf.fa"
    params:
        res = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/mtro_improvement.yaml"
    threads:
        workflow.cores
    log:
        "logs/mtro_improvement/lrscaf_mtro.log"
    shell:
        """
        java -jar {input.lrscaf} \
        -c {input.ref} \
        -a {input.paf} \
        -o {params.res}/res \
        -micl 1000 \
        -t mm \
        -i 0.3 \
        -p {threads} \
        > {log} 2>&1

        cp {params.res}/res/scaffolds.fasta {output}
        """


# ragtag
rule ragtag_mtro:
    input:
        ref_coruscus = "resources/GCA017311375.fasta.gz",
        draft_assembly = "results/mtro_02/lrscaf/mtro_01.lrscaf.fa"
    output:
        "results/mtro_02/ragtag/mtro_01.ragtag.fa"
    params:
        out_dir = lambda w, output: f"{os.path.dirname(output[0])}/res"
    conda:
        "../envs/mtro_improvement.yaml"
    threads:
        workflow.cores
    log:
        "logs/mtro_improvement/ragtag_mtro.log"
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
rule pilon_map_ont:
    input:
        ref = "results/mtro_02/ragtag/mtro_01.ragtag.fa",
        ont_reads = config['mtro_improvement']['mtro_ont'],
    output:
        "results/mtro_02/pilon/mtro_01.ragtag.ont_mapped_pilon.bam"
    log:
        "logs/mtro_improvement/map_ont_reads_pilon.log"
    conda:
        "../envs/mtro_improvement.yaml"
    threads:
        workflow.cores
    shell:
        """
        minimap2 -t {threads} -ax map-ont \
        {input.ref} {input.ont_reads} 2> {log} | \
        samtools sort -@ {threads} | \
        samtools view -b -@ {threads} -o {output}
        samtools index {output}
        """

rule prepare_bwa_index:
    input:
        "results/mtro_02/ragtag/mtro_01.ragtag.fa"
    output:
        "results/mtro_02/ragtag/mtro_01.ragtag.fa.0123"
    conda:
        "../envs/mtro_improvement.yaml"
    shell:
        """
        bwa-mem2 index {input}
        """

rule pilon_map_pe:
    input:
        ref = "results/mtro_02/ragtag/mtro_01.ragtag.fa",
        R1 = "results/preprocessing/tros/tros_dedup_proc_fastp_R1_001.fastq.gz",
        R2 = "results/preprocessing/tros/tros_dedup_proc_fastp_R2_001.fastq.gz",
        index = "results/mtro_02/ragtag/mtro_01.ragtag.fa.0123"
    output:
        "results/mtro_02/pilon/mtro_01.ragtag.pe_mapped_pilon.bam"
    log:
        "logs/mtro_improvement/map_pe_reads_pilon.log"
    conda:
        "../envs/mtro_improvement.yaml"
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
        "results/mtro_02/pilon/pilon-1.24.jar"
    params:
        link = "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar"
    shell:
        "wget -O {output} {params.link}"


rule split_in_targets:
    input:
        ref = "results/mtro_02/ragtag/mtro_01.ragtag.fa"
    output:
        expand("results/mtro_02/pilon/target_{num}.txt",
            num=[f'{i+1:02}' for i in range(15)])
    params:
        out_dir = lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/mtro_improvement.yaml"
    shell:
        """
        for i in {{1..14}}; do
            grep '>' {input.ref} | sed 's/>//' | awk -v record=$i 'NR==record {{print $0}}' \
            > {params.out_dir}/target_$(printf '%02d' $i).txt
        done
        grep '>' {input.ref} | tail -n +15 > {output[14]}
        """

rule pilon_mtro_02:
    input:
        target = "results/mtro_02/pilon/target_{num}.txt",
        ref = "results/mtro_02/ragtag/mtro_01.ragtag.fa",
        pilon = "results/mtro_02/pilon/pilon-1.24.jar",
        pe_bam = "results/mtro_02/pilon/mtro_01.ragtag.pe_mapped_pilon.bam",
        ont_bam = "results/mtro_02/pilon/mtro_01.ragtag.ont_mapped_pilon.bam",
    output:
        "results/mtro_02/pilon/res_mtro_02_{num}/mtro_02_{num}.fasta"
    params:
        out_dir = lambda w, output: os.path.dirname(output[0]),
        prefix = lambda w, output: os.path.basename(output[0]).replace('.fasta', ''),
        java_mem = config['mtro_improvement']['java_mem'],
    log:
        "logs/mtro_improvement/pilon_run_target_{num}.log"
    conda:
        "../envs/mtro_improvement.yaml"
    threads: 
        math.floor(workflow.cores/config['mtro_improvement']['max_concurrent_pilon'])
        # limit number of concurrent processes due to available RAM
    shell:
        """
        java -Xmx{params.java_mem}G -jar {input.pilon} \
        --genome {input.ref} \
        --frags {input.pe_bam} \
        --nanopore {input.ont_bam} \
        --targets {input.target} \
        --output {params.prefix} \
        --outdir {params.out_dir} \
        --changes --vcf --tracks --diploid --fix all \
        > {log} 2>&1
        """

rule merge_pilon_res:
    input:
        expand("results/mtro_02/pilon/res_mtro_02_{num}/mtro_02_{num}.fasta",
            num=[f'{i+1:02}' for i in range(15)]),
    output:
        "results/fasta/tros_v7.pseudohap.fasta.gz"
    shell:
        "cat {input} | bgzip -c > {output}"