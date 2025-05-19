#!/usr/bin/env python3
import os
import sys
import json
import argparse
import pandas as pd

def process_classyfire(tsv_path):
    """
    Reads a ClassyFire TSV, splits sid into InChIKey_smiles_firstBlock,
    drops the original sid/smiles, fixes column names, melts & explodes ;
    entries, and tags with tool='classyfire'.
    """
    # load
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)

    # first block of sid
    df['InChIKey_smiles_firstBlock'] = df['sid'].str.partition('-')[0]

    # drop original
    df = df.drop(columns=['sid', 'smiles'])

    # fix spelling / normalize column names
    rename_map = {
        'superklass': 'superclass',
        'klass':      'class',
        'subklass':   'subclass',
        'alternative parents': 'alternative_parents'
    }
    df = df.rename(columns=rename_map)

    # melt into (InChIKey_smiles_firstBlock, variable, classification)
    id_vars = ['InChIKey_smiles_firstBlock']
    value_vars = [c for c in df.columns if c not in id_vars]
    melted = df.melt(id_vars=id_vars,
                     value_vars=value_vars,
                     var_name='variable',
                     value_name='classification')

    # drop missing or empty
    melted = melted.dropna(subset=['classification'])
    melted['classification'] = melted['classification'].astype(str).str.strip()
    melted = melted[melted['classification'] != '']

    # split on ';' and explode
    melted['classification'] = melted['classification'].str.split(';')
    melted = melted.explode('classification')
    melted['classification'] = melted['classification'].str.strip()

    # drop any empties, dedupe, tag
    melted = melted[melted['classification'] != '']
    melted = melted.drop_duplicates()
    melted['tool'] = 'classyfire'

    # Make unique by InChIKey_smiles_firstBlock, and classification
    melted = melted.drop_duplicates(subset=['InChIKey_smiles_firstBlock', 'classification'])
    melted = melted.reset_index(drop=True)

    return melted


def process_npclassifier(root_dir):
    """
    Walks root_dir for .json files named <InChIKeyBlock>.json,
    pulls out class_results, superclass_results, pathway_results, isglycoside,
    and returns a DataFrame with the same three columns + tool='npclassifier'.
    """
    records = []
    key_map = {
        'superclass_results': 'superclass',
        'class_results':      'class',
        'pathway_results':    'pathway',
        'isglycoside':        'isglycoside'
    }

    for dirpath, _, files in os.walk(root_dir):
        for fn in files:
            if not fn.lower().endswith('.json'):
                continue
            inchikey_block = os.path.splitext(fn)[0]
            fullpath = os.path.join(dirpath, fn)
            try:
                with open(fullpath, 'r', encoding='utf-8') as fh:
                    data = json.load(fh)
            except Exception as e:
                sys.stderr.write(f"Warning: failed to read {fullpath}: {e}\n")
                continue

            for json_key, varname in key_map.items():
                if json_key not in data:
                    continue
                val = data[json_key]
                # list -> multiple rows
                if isinstance(val, list):
                    for item in val:
                        records.append({
                            'InChIKey_smiles_firstBlock': inchikey_block,
                            'variable': varname,
                            'classification': str(item).strip(),
                            'tool': 'npclassifier'
                        })
                # scalar -> one row
                else:
                    records.append({
                        'InChIKey_smiles_firstBlock': inchikey_block,
                        'variable': varname,
                        'classification': str(val).strip(),
                        'tool': 'npclassifier'
                    })

    df = pd.DataFrame.from_records(records)
    # remove any blanks / duplicates
    df = df[df['classification'] != '']
    df = df.drop_duplicates()

    # adapt glycosides 
    df = df[~((df['variable'] == 'isglycoside') & (df['classification'].isin(['False', 'false'])))]
    mask = (df['variable'] == 'isglycoside') & (df['classification'].isin(['True', 'true']))
    df.loc[mask, 'classification'] = 'glycoside'

    # make unique by InChIKey_smiles_firstBlock, and classification
    melted = melted.drop_duplicates(subset=['InChIKey_smiles_firstBlock', 'classification'])
    melted = melted.reset_index(drop=True)

    return df

def process_molprops(json_root_dir):
    """
    Walk through json_root_dir (including subdirs), read each `.json` file whose name
    is the first block of an InChIKey, and assemble selected molecular properties
    into a pandas DataFrame.

    Returns:
        pd.DataFrame with columns:
          - InChIKey_smiles_firstBlock
          - alogp
          - rotatable_bond_count
          - hydrogen_bond_acceptors
          - hydrogen_bond_donors
          - aromatic_rings_count
          - qed_drug_likeliness
          - formal_charge
          - number_of_minimal_rings
          - van_der_waals_volume
          - nplikeness
    """
    # list of keys we want from each JSON
    props = [
        'alogp',
        'rotatable_bond_count',
        'hydrogen_bond_acceptors',
        'hydrogen_bond_donors',
        'aromatic_rings_count',
        'qed_drug_likeliness',
        'formal_charge',
        'number_of_minimal_rings',
        'van_der_waals_volume',
        'nplikeness'
    ]

    records = []
    for dirpath, _, files in os.walk(json_root_dir):
        for fname in files:
            if not fname.lower().endswith('.json'):
                continue
            inchikey_block = os.path.splitext(fname)[0]
            fullpath = os.path.join(dirpath, fname)
            try:
                with open(fullpath, 'r', encoding='utf-8') as f:
                    data = json.load(f)
            except Exception as e:
                # skip files that can't be parsed
                continue

            rec = {'InChIKey_smiles_firstBlock': inchikey_block}
            for key in props:
                rec[key] = data.get(key, None)
            records.append(rec)

    df = pd.DataFrame.from_records(records)
    return df

def make_classification_table(classyfire_tsv, npclassifier_dir):
    """
    Reads a ClassyFire TSV and NPClassifier JSONs, and returns a
    DataFrame with the same three columns + tool='classyfire' or
    'npclassifier'.
    """
    cf = process_classyfire(classyfire_tsv)
    npc = process_npclassifier(npclassifier_dir)

    # combine
    combined = pd.concat([cf, npc], ignore_index=True)
    return combined


def main():
    p = argparse.ArgumentParser(
        description="Merge ClassyFire TSV and NPClassifier JSONs into one table"
    )
    p.add_argument('--classyfire_tsv',
                   help="Path to your ClassyFire TSV file")
    p.add_argument('--npclassifier_dir',
                   help="Path to directory (and subfolders) of NPClassifier JSONs")
    p.add_argument('--output',
                   help="Path to output TSV file (default: stdout)")
    args = p.parse_args()

    cf = process_classyfire(args.classyfire_tsv)
    npc = process_npclassifier(args.npclassifier_dir)

    combined = pd.concat([cf, npc], ignore_index=True)
    # write out as TSV to stdout
    combined.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
