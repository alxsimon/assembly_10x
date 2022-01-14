
# Mask some adaptor errors that were detected by NCBI
# on first try of submission.
rule mask_adaptor_errors:
    input:
        "results/final/{asm}/{asm}.fa.gz",
        "resources/ncbi_errors/{asm}_NCBI_errors_adaptors.bed"
    output:
        temp("results/ncbi_sub/{asm}/{asm}_masked.fa.gz")
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        zcat {input[0]} | \
        bedtools maskfasta -fi /dev/stdin -bed {input[1]} -fo /dev/stdout | \
        gzip -c > {output}
        """
    
# Remove duplicate sequences as identified by NCBI
# on first try of submission.
rule duplicate_removal:
    input:
        "results/ncbi_sub/{asm}/{asm}_masked.fa.gz"
    output:
        temp("results/ncbi_sub/{asm}/{asm}_masked_undup.fa.gz")
    params:
        "results/ncbi_sub/{asm}/seqkit_out_duplicated_{asm}.txt"
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        seqkit rmdup -s -D {params} -i \
        -o {output} {input}
        """


# Scaffolds need to be renamed and information added
# for the first 14 scaffolds that are chromosome anchored.

rule rename_scaffolds:
    input:
        "results/ncbi_sub/{asm}/{asm}_masked_undup.fa.gz"
    output:
        temp("results/ncbi_sub/{asm}/{asm}_masked_undup_renamed_1.fa.gz")
    params:
        scaff_prefix = lambda w: f"{w.asm}",
        nr_width = 5,
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        seqkit replace -p .+ \
        -r "{params.scaff_prefix}_s{{nr}}" --nr-width {params.nr_width} \
        -o {output} {input}
        """

rule add_chr_info:
    input:
        "results/ncbi_sub/{asm}/{asm}_masked_undup_renamed_1.fa.gz"
    output:
        "results/ncbi_sub/{asm}/{asm}.fa.gz"
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        seqkit head -n 14 {input} | \
        seqkit replace -p "(.+)" -r "\$1 [location=chromosome] [chromosome=LG{{nr}}]" \
        --nr-width 2 \
        -o {output}

        seqkit fx2tab {input} | tail -n +15 | seqkit tab2fx | 
        gzip -c >> {output}
        """


rule create_name_correspondence_table:
    input:
        "results/ncbi_sub/{asm}/{asm}_masked_undup.fa.gz",
        "results/ncbi_sub/{asm}/{asm}.fa.gz",
    output:
        "results/ncbi_sub/{asm}/{asm}_scaffold_name_correspondence.txt"
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        paste \
        <(seqkit seq --name --only-id {input[0]}) \
        <(seqkit seq --name --only-id {input[1]}) \
        > {output}
        """

# Replace by new scaffold name in gff3
rule replace_seq_name_gff:
    input:
        "results/ncbi_sub/{asm}/{asm}_scaffold_name_correspondence.txt",
        "results/final/{asm}/{asm}.gff3.gz",
    output:
        "results/ncbi_sub/{asm}/{asm}.gff3.gz"
    conda:
        "../envs/ncbi_sub.yaml"
    shell:
        """
        awk -F'\\t' 'NR==FNR{{a[$1]=$2; next}}{{id=$1; if(id in a) gsub(id,a[id])}} 1' \
        {input[0]} <(zcat {input[1]}) | \
        gzip -c > {output}
        """

rule copy_rest_annotation:
    input:
        "results/final/{asm}/{asm}_pep.fa.gz",
        "results/final/{asm}/{asm}_consensus_annotation.tsv",
    output:
        "results/ncbi_sub/{asm}/{asm}_pep.fa.gz",
        "results/ncbi_sub/{asm}/{asm}_consensus_annotation.tsv",
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """