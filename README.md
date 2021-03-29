# Assembly pipeline for Mytilus genomes

[`snakemake`](https://snakemake.readthedocs.io/en/stable/) (in a conda environnement for example) and 
[`singularity`](https://github.com/hpcng/singularity) need to be installed.

The supernova results are stored on a distant NAS that needs to be mounted first on my system.
```
sshfs nas4:/share/sea/sea/projects/ref_genomes/assembly_10x/results/supernova_assemblies \
results/supernova_assemblies \
-o idmap=user,compression=no,uid=1000,gid=1000,allow_root
```

I also use a 4T disk as a temporary local storage for supernova computation
`sudo mount /dev/sd[x]1 /data/ref_genomes/assembly_10x/tmp`

To run use:
```
conda activate snake_env

snakemake --use-conda --conda-frontend mamba --conda-prefix .conda \
--use-singularity --singularity-args "-B /nas_sea:/nas_sea" \
-j {threads}
```

Final versions are *_v6.pseudohap.fasta.gz and they correspond to:
- mgal_01
- medu_01
- mtro_01

Another version of mtro is done, tros_v7, also called mtro_02 which is improved by LRScaf with nanopore reads, scaffolding on the *Mytilus coruscus* reference genome and Pilon corrections.

```
conda activate snake_env

snakemake --use-conda --conda-frontend mamba --conda-prefix .conda \
--use-singularity --singularity-args "-B /nas_sea:/nas_sea" \
-j {threads} mtro_improvement
```