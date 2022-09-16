[![DOI](https://zenodo.org/badge/258168051.svg)](https://zenodo.org/badge/latestdoi/258168051)

# Assembly pipeline for Mytilus genomes

Assembly pipeline from 10x chromium reads from the preprint
"Three new genome assemblies of blue mussel lineages: North and South European Mytilus edulis and Mediterranean Mytilus galloprovincialis" bioRxiv ([https://doi.org/10.1101/2022.09.02.506387](https://doi.org/10.1101/2022.09.02.506387 )).

[`snakemake`](https://snakemake.readthedocs.io/en/stable/) (in a conda environnement for example) and 
[`singularity`](https://github.com/hpcng/singularity) need to be installed.

## Supernova storage workarounds

Supernova use large amount of storage for temporary and final results.

The supernova results are stored on a distant NAS that needs to be mounted first on my system.
```
sshfs nas4:/share/sea/sea/projects/ref_genomes/assembly_10x/results/supernova_assemblies \
results/supernova_assemblies \
-o idmap=user,compression=no,uid=1000,gid=1000,allow_root
```

I also used a 4T disk as a temporary local storage for supernova computation
`sudo mount /dev/sd[x]1 /data/ref_genomes/assembly_10x/tmp`


## How to run

To run use:
```
conda activate snake_env

snakemake --use-conda \
--use-singularity --singularity-args "-B /nas_sea:/nas_sea" \
-j {threads} \
[either all_v6, asm_improvement, stats, repeats, annotation, finalize or ncbi_submission (see workflow/Snakefile)]
```

