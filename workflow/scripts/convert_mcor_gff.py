#!python3

in_gff = snakemake.input[0]
transfer_file = snakemake.input[1]
out_gff = snakemake.output[0]

transfer = dict() # dictionary old_name: new_name
with open(transfer_file, 'r') as fr:
    for line in fr:
        L = line.split("\t")
        if not line.startswith('#'):
            transfer[L[0]] = L[4]


with open(in_gff, 'r') as fr:
    with open(out_gff, 'w') as fw:
        for line in fr:
            if not line.startswith('\n'):
                contig = line.split('\t')[0]
                new_contig = transfer[contig]
                if new_contig:
                    new_line = line.replace(contig, new_contig)
                else:
                    f"Unkown contig {contig}"
                fw.write(new_line)

