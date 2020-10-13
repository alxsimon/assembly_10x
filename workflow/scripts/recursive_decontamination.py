#! /usr/bin/env python

genome = snakemake.input[0]
clusters = snakemake.input[1:]
N = len(clusters)
print(f'Analysing {N} clusters')

n_seq={}
for clust_file in clusters:
    n = 0
    with open(clust_file, 'r') as fr:
        for line in fr:
            if line.startswith('>'):
                n += 1
    n_seq[clust_file] = n

sorted_clust = sorted(n_seq, reverse=True)
print(sorted_clust)

# the number of sequences blasted on Mytilus contigs is taken
# as proxy for the host cluster
host = sorted_clust[0]
conta = sorted_clust[1:]

tmp_genome = genome
for cl in sorted_clust:
    kount_cmd = f'\
        Kount.py -u {snakemake.threads} \
        -i {tmp_genome} -r {host} -c {cl} \
        -d {snakemake.params.dist} \
        -W {snakemake.params.wd}'
    shell(kount_cmd)

    # run contalocate on previous filtration round
    contalocate_cmd = f'\
        ../scripts/contalocate.R \
        -i {tmp_genome} -r {host} -c {cl}'
    shell(contalocate_cmd)

    conta_gff = f'{tmp_genome}_contaminant_{cl}.gff'

    #filter on the output
    new_genome = f'{tmp_genome}-{cl}'
    seqkit_cmd = f'\
        seqkit grep -f <(tail -n +2 {conta_gff} | cut -f 1) \
        {tmp_genome} > {new_genome}'
    shell(seqkit_cmd)

    tmp_genome = new_genome

shell(f'touch {snakemake.output[1]}')