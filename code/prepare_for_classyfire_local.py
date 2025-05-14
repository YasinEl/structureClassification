#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from rdkit import Chem

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract and clean SMILES → InChIKey table from a larger CSV"
    )
    p.add_argument("input_csv", help="Path to input CSV with columns Smiles and InChIKey_smiles")
    p.add_argument("output_csv", help="Path to write cleaned CSV (inchikey,smiles)")
    p.add_argument(
        "--exclude-tsv",
        dest="exclude_tsv",
        help="Optional TSV with column 'sid' of InChIKeys to exclude",
        default=None,
    )
    return p.parse_args()

def load_excludes(path):
    if path and os.path.exists(path):
        df = pd.read_csv(path, sep="\t", usecols=["sid"])
        return set(df["sid"].dropna().astype(str))
    return set()

def main():
    args = parse_args()

    # 1) load main table
    df = pd.read_csv(args.input_csv, dtype=str, usecols=["Smiles", "InChIKey_smiles"])
    df = df.rename(columns={"Smiles": "smiles", "InChIKey_smiles": "inchikey"})

    # 2) filter: valid SMILES & inchikey format while handling missing values ie floats
    def valid_row(row):
        ik = row.inchikey
        smi = row.smiles

        # must be a non‐null string containing at least one dash
        if pd.isna(ik) or not isinstance(ik, str) or "-" not in ik:
            return False

        # must be a non‐null string that RDKit can parse
        if pd.isna(smi) or not isinstance(smi, str) or Chem.MolFromSmiles(smi) is None:
            return False

        return True
    
    df = df[df.apply(valid_row, axis=1)]

    # 3) drop duplicates by inchikey, keep first
    df = df.drop_duplicates(subset="inchikey", keep="first")

    # 4) exclude any inchikey already in the sid‑tsv
    excludes = load_excludes(args.exclude_tsv)
    if excludes:
        df = df[~df["inchikey"].isin(excludes)]

    # 5) reorder columns and write out
    df = df[["inchikey", "smiles"]]
    df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()
