rule btk_prepare_workdir:
    input:
        fa = "results/fasta/{sample}_{version}.pseudohap.fasta.gz",
        fastqs = multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        fa = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta",
        GM1 = "results/blobtoolkit/{sample}_{version}/GM_1.fastq.gz",
        GM2 = "results/blobtoolkit/{sample}_{version}/GM_2.fastq.gz"
    conda: 
        "../envs/btk_env.yaml"
    shell:
        """
        zcat {input.fa} | seqkit grep -v -s -r -p '^N+$' > {output.fa}
        ln -sr {input[1]} results/blobtoolkit/{wildcards.sample}_{wildcards.version}/GM_1.fastq.gz
        ln -sr {input[2]} results/blobtoolkit/{wildcards.sample}_{wildcards.version}/GM_2.fastq.gz
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
        sed -i "s/scaffold-count.*/scaffold-count:\ $SC/" {output.conf}
        sed -i "s/span.*/span:\ $SP/" {output.conf}
        sed -i "s/prefix.*/prefix:\ {wildcards.sample}_{wildcards.version}/" {output.conf}
        """

rule btk_insdc_pipeline:
    input:
        fa = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta",
        conf = "results/blobtoolkit/{sample}_{version}/{sample}_{version}.yaml",
        GM1 = "results/blobtoolkit/{sample}_{version}/GM_1.fastq.gz",
        GM2 = "results/blobtoolkit/{sample}_{version}/GM_2.fastq.gz"
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
        --conda-prefix .conda \
        --configfile {input.conf} \
        -j {threads} \
        --latency-wait 30 \
        --resources btk=1
        """

