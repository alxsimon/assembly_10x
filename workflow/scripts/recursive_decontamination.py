#! /usr/bin/env python

from snakemake.shell import shell
import os

win_step = snakemake.params.win_step
win_size = snakemake.params.win_size
pattern = snakemake.params.pattern

genome = snakemake.input[0]
clusters = snakemake.input[1:]

n_seq={}
for clust_file in clusters:
    n = 0
    with open(clust_file, 'r') as fr:
        for line in fr:
            if line.startswith('>'):
                n += 1
    n_seq[clust_file] = n

sorted_clust = sorted(n_seq, key=n_seq.__getitem__, reverse=True)

# the number of sequences is taken
# as proxy for the host cluster
# and the order of the higher to lower contamination
host = sorted_clust[0]
conta = sorted_clust[1:]

N = len(conta)
print(f'Analysing {N} clusters of potential contaminants')

tmp_genome = genome
for cl_file in conta:
    cl = os.path.basename(cl_file).replace('data_fasta_', '').replace('.fa', '')
    suffix_genome = os.path.basename(tmp_genome)

    outfile = f'{snakemake.params.wd}{suffix_genome}.mcp_hostwindows_vs_conta_data_fasta_{cl}.fa_KL.dist'
    if os.path.exists(outfile):
        print(f'Kount.py outputs already exist. Moving on to contalocate.')
    else:
        kount_cmd = f'\
            Kount.py -u {snakemake.threads} \
            -i {tmp_genome} -r {host} -c {cl_file} \
            -t {win_step} -w {win_size} -p {pattern} \
            -d {snakemake.params.dist} \
            -W {snakemake.params.wd}'
        shell(kount_cmd)

    # run contalocate on previous filtration round
    contalocate_cmd = f'\
        workflow/scripts/contalocate.R \
        -i {tmp_genome} -r {host} -c {cl_file} \
        -t {win_step} -w {win_size} \
        -d {snakemake.params.dist} \
        -W {snakemake.params.wd}'
    shell(contalocate_cmd)

    conta_gff = f'{snakemake.params.wd}/{suffix_genome}_contaminant_data_fasta_{cl}.fa.gff'

    #filter on the output
    tmp_prefix = tmp_genome.replace('.fa', '')
    new_genome = f'{tmp_prefix}-{cl}.fa'
    with open(conta_gff) as fr:
        n_lines_gff = sum(1 for x in fr)
    if (n_lines_gff > 1):
        seqkit_cmd = f'\
            seqkit grep -v -f <(tail -n +2 {conta_gff} | cut -f 1) \
            {tmp_genome} > {new_genome}'
        shell(seqkit_cmd)
    else:
        shell(f'cp {tmp_genome} {new_genome}')

    tmp_genome = new_genome

shell(f'touch {snakemake.output[0]}')

shell(f'gzip -c {tmp_genome} > {snakemake.output[1]}')
shell(f'cat {snakemake.params.wd}/*.gff \
    | grep -v "##gff-version 2" > {snakemake.output[3]}')
shell(f'seqkit grep -f <(cut -f 1 {snakemake.output[3]}) \
    | gzip -c > {snakemake.output[2]}')
