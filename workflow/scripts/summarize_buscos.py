#!python3
import pathlib
import defopt
import sys
import pandas as pd

def parse_busco(filename):
    """
    Parse one short_summary BUSCO file
    """
    asm = filename.name.split('.')[3]
    db = filename.name.split('.')[2]
    data = {'asm': asm, 'db': db}
    with open(filename, 'r') as fr:
        for line in fr:
            if "Complete and single-copy BUSCOs" in line:
                data['CS'] = int(line.split("\t")[1])
            elif "Complete and duplicated BUSCOs" in line:
                data['CD'] = int(line.split("\t")[1])
            elif "Fragmented BUSCOs" in line:
                data['F'] = int(line.split("\t")[1])
            elif "Missing BUSCOs" in line:
                data['M'] = int(line.split("\t")[1])
    data['C'] = data['CS'] + data['CD']
    data['T'] = data['CS'] + data['CD'] + data['F'] + data['M']
    return data


def format_output(out, datasets):
    df = pd.DataFrame(datasets)
    df = df.melt(id_vars=['asm', 'db'], 
        value_vars=['CS', 'CD', 'F', 'M', 'C', 'T'],
        var_name='cat', value_name='N')
    if out is None:
        df.to_csv(sys.stdout, sep='\t', index=False)
    else:
        df.to_csv(out, sep='\t', index=False)


def main(files: list[pathlib.Path], out: pathlib.Path=None):
    """
    Summarize a list of short summary BUSCO files into a tidy table.

    :param files: BUSCO short_summary files (space separated list)
    :param out: Summary output as tsv, stdout by default
    """
    datasets = []
    for busco_file in files:
        datasets.append(parse_busco(busco_file))

    format_output(out, datasets)        


if __name__ == '__main__':
    defopt.run(main, strict_kwonly=False,
        short={'out': 'o'})