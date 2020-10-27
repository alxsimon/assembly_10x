#! python3

import json
import pandas as pd
import numpy as np

files = {
    'id': 'identifiers.json',
    'gc': 'gc.json',
    'cov': 'GM_cov.json',
    'length': 'length.json',
    'Ncount': 'ncount.json',
    'superkingdom': 'bestsumorder_superkingdom.json',
    'phylum': 'bestsumorder_phylum.json'
}

blobdir = snakemake.params.blobdir

df = pd.DataFrame()
for var, file in files.items():
    with open(f'{blobdir}/{file}') as fr:
        data = json.load(fr)
        if data['keys'] == []:
            data = data['values']
        else:
            data = [data['keys'][x] for x in data['values']]
        df[var] = data

with open(f'{blobdir}/bestsumorder_positions.json') as fr:
    pos = json.load(fr)

df['besthit_length'] = [(x[0][2] - x[0][1] + 1) if x!=[] else np.nan for x in pos['values']]
df['taxid'] = [str(x[0][0]) if x!=[] else np.nan for x in pos['values']]
df['besthit_perc'] = df.besthit_length/df.length
df['N_perc'] = df.Ncount/df.length

with open(snakemake.params.mollusca_taxids) as fr:
    mollusca_taxids = set(int(x.strip()) for x in fr.readlines())
taxids = [[y[0] for y in x] if x!=[] else [] for x in pos['values']]
df['any_Mollusca_hit'] = [(len(set(x) & mollusca_taxids) > 0) for x in taxids]

df_bacteria = df[(df.superkingdom=='Bacteria')]
df_viruses = df[(df.superkingdom=='Viruses') & (df.besthit_perc > 0.1)]
df_eukaryota = df[
    (df.superkingdom=='Eukaryota') & 
    (~df.phylum.isin(['Mollusca', 'no-hit'])) &
    (df.any_Mollusca_hit == False) &
    (df.besthit_perc > 0.1)
]

conta_ids = df_bacteria.id.tolist() \
+ df_viruses.id.tolist() \
+ df_eukaryota.id.tolist()

df_kept = df[~df.id.isin(conta_ids)]

df_kept.to_csv(snakemake.output.kept, index=False)
df_bacteria.to_csv(snakemake.output.bact, index=False)
df_viruses.to_csv(snakemake.output.virus, index=False)
df_eukaryota.to_csv(snakemake.output.euka, index=False)