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

host = sorted_clust[0]
conta = sorted_clust[1:]

for cl in sorted_clust:
    kount_cmd = f'\
        Kount.py -u {snakemake.threads} \
        -i {genome} -r {host} -c {cl} \
        -d {snakemake.params.dist} \
        -W {snakemake.params.wd}'
    shell(kount_cmd)

    contalocate_cmd = f'\
        ../scripts/contalocate.R \
        -i {genome} -r {host} -c {cl}'
    