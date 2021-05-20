# Align the five soft masked Mytilus genomes using
# progressive cactus (gpu version)

RESTART=False

# keep only short name for ncbi retrived genomes
rule correct_seq_names:
    input:
        "resources/{genome}.fasta.gz"
    output:
        "resources/{genome}.fa.corrected"
    conda:
        "envs/seqkit.yaml"
    shell:
        "seqkit replace -p '\s.*$' -r '' {input} > {output}"

rule cactus_alignment:
    input:
        "resources/GCA017311375.fa.corrected",
        "results/repeats/tros_v7.fa.masked",
        "results/repeats/edu_v7.fa.masked",
        "results/repeats/gallo_v7.fa.masked",
        "resources/GCA900618805.fa.corrected",
    output:
        "results/cactus/myt_cactus.hal"
    params:
        seqfile = "results/cactus/seqfile.txt",
        restart = lambda w: '--restart' if RESTART else '',
    log:
        "logs/cactus_alignment.log"
    threads:
        workflow.cores
    container:
        "containers/cactus_v1.3.0-gpu.sif"
    shell:
        """
        echo "((((mgal_02,GCA900618805),medu_02),mtro_02),GCA017311375);" > {params.seqfile}
        echo "*GCA017311375 {input[0]}" >> {params.seqfile}
        echo "mtro_02 {input[1]}" >> {params.seqfile}
        echo "medu_02 {input[2]}" >> {params.seqfile}
        echo "mgal_02 {input[3]}" >> {params.seqfile}
        echo "GCA900618805 {input[4]}" >> {params.seqfile}

        cactus \
        {params.restart} \
	    --maxCores {threads} \
        results/cactusJobStore \
        {params.seqfile} \
        {output} \
        > {log} 2>&1
        """
