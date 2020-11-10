#! /usr/bin/env python

import pandas as pd
import json
import re

out_df = pd.DataFrame()
for stat_file in snakemake.input:
    m = re.search('results/stats/(.+?)_(.+?).stats.json', stat_file)
    sample = m.group(1)
    version = m.group(2)
    with open(stat_file, 'r') as f:
        data = json.loads(f.read())
    data['C'] = data.pop('Contig Stats')
    data['S'] = data.pop('Scaffold Stats')
    tmpdf = pd.json_normalize(data)
    tmpdf.insert(0, "asm", [f"{sample}_{version}"])
    out_df = out_df.append(tmpdf)

out_df.to_csv(snakemake.output[0], sep=',', index=False)