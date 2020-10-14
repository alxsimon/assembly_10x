#!python

import sys
import json
import pandas as pd
import math

# for now args
# - input blobdir
# - output prefix
blobdir = sys.argv[1]
out_prefix = sys.argv[2]

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
    with open(f'{blobdir}/{file}') as f:
        data = json.load(f)
        if data['keys'] == []:
            data = data['values']
        else:
            data = [data['keys'][x] for x in data['values']]
        df[var] = data

mean_gc = df.gc.mean()
sd_gc = math.sqrt(df.gc.var())
lim_gc = mean_gc + 2*sd_gc

# output lists of identifiers for each interesting category
# Mollusca
df[df.phylum=='Mollusca'].to_csv(
    f'{out_prefix}_mollusca.id', 
    columns=['id'], header=False, index=False)
# Bacteria
df[(df.superkingdom=='Bacteria') & (df.gc > lim_gc)].to_csv(
    f'{out_prefix}_bacteria.id', 
    columns=['id'], header=False, index=False)
# Viruses
df[(df.superkingdom=='Viruses') & (df.gc > lim_gc)].to_csv(
    f'{out_prefix}_viruses.id', 
    columns=['id'], header=False, index=False)