rule kat_comp:
    input:
        expand("results/preprocessing/{{sample}}/{{sample}}_S1_L001_{R}_001.fastq.gz", 
            R=["R1", "R2"]),
        "results/fasta/{sample}_{version}.pseudohap.fasta.gz"
    output:
        "results/kat/{sample}_{version}/{sample}_{version}_comp-main.mx"
    params:
        outprefix = lambda w: f'results/kat/{w.sample}_{w.version}/{w.sample}_{w.version}_comp'
    log:
        "logs/kat_comp.{sample}_{version}.log"
    conda:
        "../envs/kat.yaml"
    threads:
        16
    shell:
        """
        kat comp -t {threads} \
        -o {params.outprefix} \
        '{input[0]} {input[1]}' \
        {input[2]} \
        > {log} 2>&1
        """

rule kat_plot_spectra:
    input:
        "results/kat/{sample}_{version}/{sample}_{version}_comp-main.mx"
    output:
        "results/kat/{sample}_{version}/{sample}_{version}_spectra.pdf"
    params:
        title = lambda w: f'{w.sample}_{w.version}'
    log:
        "logs/kat_plot.{sample}_{version}.log"
    conda:
        "../envs/kat.yaml"
    shell:
        """
        kat plot spectra-cn \
        -o {output} \
        -t {params.title} \
        {input} \
        > {log} 2>&1
        """