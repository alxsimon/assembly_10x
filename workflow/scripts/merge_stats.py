#! /usr/bin/env python

import pandas as pd
from pathlib import Path
import json

out_df = pd.DataFrame()
for stat_file in snakemake.input:
    with open(stat_file, 'r') as f:
        data = json.loads(f.read())
    data['C'] = data.pop('Contig Stats')
    data['S'] = data.pop('Scaffold Stats')
    tmpdf = pd.json_normalize(data)
    tmpdf.insert(0, "asm", [f"{snakemake.wildcards.sample}_{snakemake.wildcards.version}"])
    out_df = out_df.append(tmpdf)

out_df.to_csv(snakemake.output[0], sep=',', index=False)