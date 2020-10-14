#! /usr/bin/env python

from snakemake.shell import shell

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
for cl in conta:
    kount_cmd = f'\
        Kount.py -u {snakemake.threads} \
        -i {tmp_genome} -r {host} -c {cl} \
        -t {win_step} -w {win_size} -p {pattern} \
        -d {snakemake.params.dist} \
        -W {snakemake.params.wd}'
    shell(kount_cmd)

    # run contalocate on previous filtration round
    contalocate_cmd = f'\
        workflow/scripts/contalocate.R \
        -i {tmp_genome} -r {host} -c {cl} \
        -t {win_step} -w {win_size} \
        -d {snakemake.params.dist} \
        -W {snakemake.params.wd}'
    shell(contalocate_cmd)

    conta_gff = f'{tmp_genome}_contaminant_{cl}.gff'

    #filter on the output
    new_genome = f'{tmp_genome}-{cl}'
    seqkit_cmd = f'\
        seqkit grep -f <(tail -n +2 {conta_gff} | cut -f 1) \
        {tmp_genome} > {new_genome}'
    shell(seqkit_cmd)

    tmp_genome = new_genome

shell(f'touch {snakemake.output[0]}')

shell(f'gzip -c {tmp_genome} > {snakemake.output[1]}')
shell(f'cat {snakemake.params.wd}*.gff > {snakemake.output[3]}')
shell(f'seqkit grep -f {snakemake.output[3]} | gzip -c > {snakemake.output[2]}')