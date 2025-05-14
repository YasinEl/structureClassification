#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract valid SIDs from a TSV and append them to a log file"
    )
    p.add_argument(
        "input_tsv",
        help="Path to the existing TSV file with headers including 'sid' and 'smiles'"
    )
    p.add_argument(
        "output_log",
        help="Path to the log file (one‑column, no header) to append SIDs to (may or may not exist)"
    )
    return p.parse_args()

def main():
    args = parse_args()

    # 1) verify the first TSV exists
    if not os.path.exists(args.input_tsv):
        sys.exit(f"Error: input file '{args.input_tsv}' does not exist")

    # 2) load only 'sid' and 'smiles'
    df = pd.read_csv(
        args.input_tsv,
        sep="\t",
        usecols=["sid", "smiles"],
        dtype=str
    )

    # 3) filter out invalid smiles
    df = df[~df["smiles"].isin(["Invalid", "ErrorOccur"])]

    # 4) keep only 'sid'
    new_sids = df[["sid"]].dropna()

    # 5) if output_log exists, load existing SIDs and append
    if os.path.exists(args.output_log):
        old = pd.read_csv(
            args.output_log,
            sep="\t",
            header=None,
            names=["sid"],
            dtype=str
        )
        all_sids = pd.concat([old, new_sids], ignore_index=True)
    else:
        all_sids = new_sids

    # 6) write back, no header, one column
    all_sids.to_csv(
        args.output_log,
        sep="\t",
        header=False,
        index=False
    )

if __name__ == "__main__":
    main()
