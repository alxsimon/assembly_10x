# Assembly pipeline for Mytilus genomes

[`snakemake`](https://snakemake.readthedocs.io/en/stable/) (in a conda environnement for example) and 
[`singularity`](https://github.com/hpcng/singularity) need to be installed.

To run use:
```
conda activate snake_env
snakemake --use-conda --use-singularity -j {threads} \
--singularity-args "-B /nas_sea:/nas_sea"
```