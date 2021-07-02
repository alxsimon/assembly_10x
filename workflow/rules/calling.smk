


rule angsd_call:
    input:
        bams = expand("results/mapping/{sample}_v7.bam", sample=['edu', 'gallo', 'tros']),
        ref = "results/fasta/gallo_v7.pseudohap.fasta.gz",
        targets = "resources/angsd_subset.sites",
    output:
        "results/calling/genomes_angsd.bcf",
        "results/calling/genomes_angsd.beagle.gz",
    params:
        prefix = lambda w, output: output[0].replace('.bcf', '')
    container:
        "containers/angsd.sif"
    log:
        "logs/calling/angsd_call.log"
    threads:
        16
    shell:
        """
        echo '{input.bams}' | tr ' ' '\\n' > results/calling/bam.list
        angsd sites index {input.targets} > {log} 2>&1
        angsd -P {threads} \
        -bam results/calling/bam.list \
        -ref {input.ref} \
        -sites {input.targets} \
        -remove_bads 1 -uniqueOnly 1 -minMapQ 20 -minQ 20 \
        -gl 2 -doMajorMinor 3 -doGlf 2 -doBcf 1 \
        -doPost 1 -doMaf 1 -doGeno 1 -doCounts 1 \
        -setMinChunkSize 10000 \
        -out {params.prefix} \
        >> {log} 2>&1
        """

rule merge_beagle:
    input:
        "resources/angsd_ref_subset.beagle.gz",
        "results/calling/genomes_angsd.beagle.gz",
    output:
        "results/calling/merged_angsd.beagle.gz"
    container:
        "containers/angsd.sif"
    shell:
        """
        paste <(zcat {input[0]}) <(zcat {input[1]} | cut -f 4-) | \
        bgzip -c > {output}
        """

rule pcangsd:
    input:
        "results/calling/merged_angsd.beagle.gz",
    output:
        multiext("results/calling/pcangsd",
            ".cov", ".sites", ".weights.npy")
    params:
        min_maf = 0.05,
        eigenvectors = 6,
        prefix = lambda w, output: output[0].replace('.cov', ''),
    threads:
        16
    log:
        "logs/calling/pcangsd.log"
    container:
        "containers/angsd.sif"
    shell:
        """
        python /opt/pcangsd/pcangsd.py \
        -beagle {input} \
        -minMaf {params.min_maf} \
        -iter 200 \
        -e {params.eigenvectors} \
        -o {params.prefix} \
        -sites_save \
        -snp_weights \
        -admix \
        -threads {threads} \
        > {log} 2>&1
        """

#==============================================
# try with Chris' contigs

rule index_contigs:
    input:
        "resources/Fraisse2016/final_contig_set.fasta"
    output:
        multiext("resources/Fraisse2016/final_contig_set.fasta",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    log:
        "logs/calling/index_contigs.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        bwa-mem2 index {input} > {log} 2>&1
        """

rule map_genomes_contigs:
    input:
        fa = "resources/Fraisse2016/final_contig_set.fasta",
        index = multiext("resources/Fraisse2016/final_contig_set.fasta",
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        fastqs = multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        "results/calling/{sample}_v7.bam",
        "results/calling/{sample}_v7.bam.bai"
    params:
        rg_string = lambda w: f'@RG\\tID:{w.sample}\\tSM:{w.sample}'
    log:
        "logs/calling/map_genomes_contigs_{sample}_v7.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        config['mapping']['threads']
    shell:
        """
        bwa-mem2 mem -t {threads} \
        -R \"{params.rg_string}\" \
        {input.fa} \
        {input.fastqs} 2> {log} | \
        samtools sort -@ {threads} -o {output[0]} -
        samtools index {output[0]}
        """

rule get_sites_calling:
    input:
        "resources/Fraisse2016/SNP_GQ20_0%missing.vcf"
    output:
        "results/calling/targets_SNP_contigs.pos.txt"
    conda:
        "../envs/calling.yaml"
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\n' {input} > {output}
        """

rule bcftools_call:
    input:
        bams = expand("results/calling/{sample}_v7.bam", sample=['edu', 'gallo', 'tros']),
        ref = "resources/Fraisse2016/final_contig_set.fasta",
        targets = "results/calling/targets_SNP_contigs.pos.txt",
    output:
        "results/calling/genomes_calling_contigs.bcf"
    conda:
        "../envs/calling.yaml"
    log:
        "logs/calling/bcftools_call_contigs.log"
    threads:
        4
    shell:
        """
        bcftools mpileup \
        -f {input.ref} \
        -R {input.targets} \
        --redo-BAQ -a "FORMAT/AD,FORMAT/DP" \
        -Ou {input.bams} 2> {log} | \
        bcftools call -mG - -Ob --threads {threads} -o {output} \
        >> {log} 2>&1
        """

# rule merge_bcfs:
#     input:
#         bcf1 = "results/calling/angsd_subset_fixref.bcf",
#         bcf2 = "results/calling/genomes_calling.bcf",
#         targets = "results/calling/targets_pruned.pos.txt",
#         ref = "results/fasta/gallo_v7.pseudohap.fasta.gz",
#     output:
#         "results/calling/merged_subset.bcf"
#     conda:
#         "../envs/calling.yaml"
#     log:
#         "logs/calling/merging.log"
#     threads:
#         4
#     shell:
#         """
#         bcftools index {input.bcf1}
#         bcftools index {input.bcf2}
#         bcftools merge -m snps -Ob -o {output} \
#         --threads {threads} {input.bcf1} {input.bcf2} > {log} 2>&1
#         """



# then merge with angsd_ref
# Then quick ADMIXTURE analysis? or PCA?