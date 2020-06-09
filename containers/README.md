# Containers used in the pipeline

`sudo singularity build supernova.sif supernova.def`

The `assembly_tools` container will contain the miniconda
distribution and packages not installable through conda.
`sudo singularity build assembly_tools.sif assembly_tools.def`

Other tools will be dealt with `--use-conda` directive of `snakemake`