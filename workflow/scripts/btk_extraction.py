#! python

import json
import pandas as pd
from snakemake.shell import shell

files = {
    'id': 'identifiers.json',
    'gc': 'gc.json',
    'cov': 'GM_cov.json',
    'length': 'length.json',
    'Ncount': 'ncount.json',
    'superkingdom': 'bestsumorder_superkingdom.json',
    'phylum': 'bestsumorder_phylum.json'
}

with open(snakemake.input[0]) as fr:
    potential_conta = json.load(fr)

df = pd.DataFrame()
for var, file in files.items():
    with open(f'{snakemake.input[1]}/{file}') as fr:
        data = json.load(fr)
        if data['keys'] == []:
            data = data['values']
        else:
            data = [data['keys'][x] for x in data['values']]
        df[var] = data

ids_file_conta = snakemake.output[0].replace('.fa', '.ids')
ids_file_mollusca = snakemake.output[1].replace('.fa', '.ids')

df[df.phylum.isin(potential_conta)].to_csv(
    ids_file_conta,
    columns=['ids'], header=False, index=False
)

df[df.phylum=='Mollusca'].to_csv(
    ids_file_mollusca,
    columns=['ids'], header=False, index=False
)

shell(f'seqkit grep -f {ids_file_conta} {snakemake.input[2]} \
    > {snakemake.output[0]}')

shell(f'seqkit grep -f {ids_file_mollusca} {snakemake.input[2]} \
    > {snakemake.output[1]}')