# Assembly pipeline for Mytilus genomes

`snakemake` and `singularity` need to be installed.

To run use:
```
conda activate snake_env
snakemake --use-conda --use-singularity -j {threads}
```