# need to move removal of all N elsewhere

rule btk_prepare_workdir:
    input:
        fa = "results/fasta/{sample}_{version}.pseudohap.fasta.gz",
        multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        fa = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta",
    conda: 
        "../envs/btk_env.yaml"
    shell:
        """
        seqkit grep -s -r -p '^N+$' -v {input.fa} > {output.fa}
        ln -s {input[1]} results/blobtoolkit/{wildcards.sample}_{wildcards.version}/GM_1.fastq.gz
        ln -s {input[2]} results/blobtoolkit/{wildcards.sample}_{wildcards.version}/GM_2.fastq.gz
        """

rule btk_prepare_conf:
    input:
        conf = "resources/blobtoolkit_{sample}_base.yaml",
        fa = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta"
    output:
        conf = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.yaml"
    shell:
        """
        cp {input.conf} {output.conf}
        SC=$( grep '>' {input.fa} | wc -l )
        SP=$( grep -v '>' {input.fa} | wc -m )
        sed -i "s/scaffold-count.*/scaffold-count: $SC/" {output.conf}
        sed -i "s/span.*/span: $SP/" {output.conf}
        sed -i "s/prefix.*/prefix: {wildcards.sample}_{wildcards.version}" {output.conf}
        """

rule btk_insdc_pipeline:
    input:
        fa = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta",
        conf = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.yaml"
    output:
        "results/blobtoolkit/{sample}_{version}/{sample}_{version}/meta.json"
    params:
        dir = lambda w: f'results/blobtoolkit/{w.sample}_{w.version}'   
    conda:
        "../envs/btk_env.yaml"
    threads:
        config['btk']['threads']
    shell:
        """
        snakemake -p \
        -s /opt/blobtoolkit/insdc-pipeline/Snakefile \
        --directory {params.dir} \
        --use-conda \
        --conda-prefix {params.dir}/.conda \
        --configfile {input.conf} \
        -j {threads} \
        --latency-wait 30 \
        --resources btk=1
        """

