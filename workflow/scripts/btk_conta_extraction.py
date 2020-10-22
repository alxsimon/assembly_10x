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

df = pd.DataFrame()
for var, file in files.items():
    with open(f'{snakemake.input[0]}/{file}') as fr:
        data = json.load(fr)
        if data['keys'] == []:
            data = data['values']
        else:
            data = [data['keys'][x] for x in data['values']]
        df[var] = data

with open('bestsumorder_positions.json') as f:
    pos = json.load(f)

df['besthit_length'] = [(x[0][2] - x[0][1] + 1) if x!=[] else np.nan for x in pos['values']]
df['subject'] = [x[0][4] if x!=[] else np.nan for x in pos['values']]
df['taxid'] = [str(x[0][0]) if x!=[] else np.nan for x in pos['values']]
df['besthit_perc'] = df.besthit_length/df.length
df['N_perc'] = df.Ncount/df.length

df_bacteria = df[(df.superkingdom=='Bacteria')]
df_viruses = df[(df.superkingdom=='Viruses') & (df.besthit_perc > 0.1)]
df_eukaryota = df[
    (df.superkingdom=='Eukaryota') & 
    (~df.phylum.isin(['Mollusca', 'no-hit'])) & 
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