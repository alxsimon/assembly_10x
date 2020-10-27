rule btk_prepare_workdir:
    input:
        fa = "results/fasta/{sample}_{version}.pseudohap.fasta.gz",
        fastqs = multiext("results/preprocessing/{sample}/{sample}_dedup_proc_fastp",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        fa = temp("results/blobtoolkit/{sample}_{version}/{sample}_{version}.fasta"),
        GM1 = "results/blobtoolkit/{sample}_{version}/GM_1.fastq.gz",
        GM2 = "results/blobtoolkit/{sample}_{version}/GM_2.fastq.gz"
    conda: 
        "../envs/btk_env.yaml"
    shell:
        """
        seqkit grep -v -s -r -p '^N+$' {input.fa} > {output.fa}
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

rule btk_add_busco:
    input:
        "results/blobtoolkit/{sample}_{version}/{sample}_{version}/meta.json",
        expand("results/busco/{{sample}}_{{version}}_{db}/run_{db}/full_table.tsv",
            db=["metazoa_odb10", "mollusca_odb10"])
    output:
        expand("results/blobtoolkit/{{sample}}_{{version}}/{{sample}}_{{version}}/{db}_busco.json",
            db=["metazoa_odb10", "mollusca_odb10"])
    params:
        blobtools_bin = config['btk']['blobtools_path'],
        blobdir = lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/btk_env.yaml"
    shell:
        """
        {params.blobtools_bin} add \
        --busco {input[2]}
        --busco {input[1]}
        {params.blobdir}
        """

rule btk_clean:
    input:
        "results/blobtoolkit/{sample}_{version}/{sample}_{version}/meta.json",
        "results/blobtoolkit/{sample}_{version}/{sample}_{version}.yaml",
        expand("results/blobtoolkit/{{sample}}_{{version}}/{{sample}}_{{version}}/{db}_busco.json",
            db=["metazoa_odb10", "mollusca_odb10"])
    output:
        "results/blobtoolkit/blobdirs/{sample}_{version}/meta.json",
        "results/blobtoolkit/{sample}_{version}.yaml",
        "results/blobtoolkit/{sample}_{version}_insdc_pipeline.tar.gz",
        touch("results/blobtoolkit/DONE_{sample}_{version}")
    params:
        indir = lambda w, input: os.path.dirname(input[0]),
        tardir = lambda w: f'results/blobtoolkit/{w.sample}_{w.version}'
    shell:
        """
        cp {input[1]} {output[1]}
        cp -r {params.indir} results/blobtoolkit/blobdirs/ && \
        tar -czf {output[2]} {params.tardir} && \
        rm -r {params.tardir}
        """
