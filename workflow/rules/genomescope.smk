rule gs_kmer_count:
    input:
        multiext("results/preprocessing/{sample}/{sample}_S1_L001",
            "_R1_001.fastq.gz", "_R2_001.fastq.gz")
    output:
        temp(multiext("results/genomescope/{sample}/{sample}.kmcdb",
            ".kmc_pre", ".kmc_suf"))
    params:
        mem = config['genomescope']['kmc_mem'],
        k = config['genomescope']['kmer_size'],
        db = lambda w, output: output[0].replace('.kmc_pre', ''),
        files = lambda w: f'results/genomescope/{w.sample}/FILES'
    threads:
        config['genomescope']['threads']
    log:
        "logs/gs_kmc_count.{sample}.log"
    container:
        "containers/kmer_analyses.sif"
    shell:
        """
        ls {input} > {params.files}
        kmc -k{params.k} -t{threads} -m{params.mem} \
        -ci1 -cs10000 @{params.files} \
        {params.db} /tmp/ \
        > {log} 2>&1
        """

rule gs_kmer_hist:
    input:
        multiext("results/genomescope/{sample}/{sample}.kmcdb",
            ".kmc_pre", ".kmc_suf")
    output:
        "results/genomescope/{sample}/{sample}.kmc_hist"
    params:
        db = lambda w, input: input[0].replace('.kmc_pre', '')
    threads:
        config['genomescope']['threads']
    log:
        "logs/gs_kmc_hist.{sample}.log"
    container:
        "containers/kmer_analyses.sif"
    shell:
        """
        kmc_tools -t{threads} transform \
        {params.db} histogram \
        {output} -cx10000 \
        > {log} 2>&1
        """

rule gs_fit:
    input:
        "results/genomescope/{sample}/{sample}.kmc_hist"
    output:
        "results/genomescope/{sample}/genomescope_res_{sample}/{sample}_preproc_summary.txt"
    params:
        k = config['genomescope']['kmer_size'],
        name_prefix = lambda w: f'{w.sample}_preproc',
        outdir = lambda w, output: os.path.dirname(output[0])
    threads:
        config['genomescope']['threads']
    log:
        "logs/gs_fit.{sample}.log"
    container:
        "containers/kmer_analyses.sif"
    shell:
        """
        /opt/genomescope2.0/genomescope.R -i {input} \
        -o {params.outdir} -n {params.name_prefix} \
        -k {params.k} -p 2 \
        > {log} 2>&1
        """