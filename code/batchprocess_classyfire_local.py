#!/usr/bin/env python3

import argparse
import csv
import os
import time
import requests
from tqdm import tqdm

CLASSYFIRE_URL = "http://classyfire.wishartlab.com"

def get_classyfire_json(inchikey):
    """
    Retrieve ClassyFire classification (JSON) for the given InChIKey.
    
    :param inchikey: The chemical's InChIKey.
    :type inchikey: str
    :return: JSON string returned by ClassyFire for that InChIKey.
    :rtype: str
    """
    # Remove any 'InChIKey=' prefix if present
    inchikey = inchikey.replace("InChIKey=", "")

    response = requests.get(f"{CLASSYFIRE_URL}/entities/{inchikey}.json")
    response.raise_for_status()
    return response.text

def main():
    parser = argparse.ArgumentParser(
        description="Fetch ClassyFire JSON for a list of InChIKeys in a CSV file."
    )
    parser.add_argument(
        "input_csv",
        help="Path to the input CSV file (must contain a column named 'InChIKey_smiles')."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to store the downloaded JSON files. Subfolders will be created."
    )
    parser.add_argument(
        "--logfile",
        required=True,
        help="Path to a log file where successfully processed InChIKey prefixes are stored."
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=10.0,
        help="Maximum requests per second (default: 10)."
    )
    args = parser.parse_args()

    # 1) Read or create the log file to track what's already downloaded
    processed_keys = set()
    if os.path.exists(args.logfile):
        with open(args.logfile, "r", encoding="utf-8") as log_fh:
            for line in log_fh:
                cleaned = line.strip()
                if cleaned:
                    processed_keys.add(cleaned)

    # Keep a counter so we know how many files have been *saved* so far
    # (this influences which subfolder to use)
    files_saved = len(processed_keys)

    # Ensure the output directory exists
    os.makedirs(args.outdir, exist_ok=True)


    # 2) Read input TSV and process each InChIKey
    with open(args.input_csv, "r", encoding="utf-8") as csv_fh:
        reader = list(csv.DictReader(csv_fh, delimiter=","))
        reversed_reader = reader[::-1]

    # setup progress bar
    progress_bar = tqdm(total=len(reader), desc="Processing rows", unit="rows")

    for row in reversed_reader:
        progress_bar.update(1)
        if "InChIKey_smiles" not in row:
            # Skip if the TSV doesn't have the 'inchikey' column
            continue

        full_inchikey = row["InChIKey_smiles"]
        if not full_inchikey:
            # Skip blank lines or missing InChIKeys
            continue

        # skip Inchikeys which are not in the expected format
        if "-" not in full_inchikey or len(full_inchikey.split("-")) != 3:
            print(f"[ERROR] Invalid InChIKey format: {full_inchikey}. Skipping.")
            continue

        # Use only the first part of the InChIKey as the “short name”
        # e.g. ABCDEFGHIJKLMA for ABCDEFGHIJKLMA-UHFFFAOYSA-N
        short_inchikey = full_inchikey.split("-", 1)[0]

        # If we've already processed this short_inchikey, skip
        if short_inchikey in processed_keys:
            continue

        # Attempt to fetch the JSON from ClassyFire
        try:
            json_text = get_classyfire_json(full_inchikey)
        except requests.RequestException as e:
            print(f"[ERROR] Could not fetch data for {full_inchikey}: {e}")
            time.sleep(3 / args.rate_limit)
            continue

        # Figure out which subfolder to store this file in
        # so as not to exceed 100 files per folder.
        folder_index = files_saved // 100
        subfolder_path = os.path.join(args.outdir, str(folder_index))
        os.makedirs(subfolder_path, exist_ok=True)

        output_filename = f"{short_inchikey}.json"
        output_path = os.path.join(subfolder_path, output_filename)

        # Write the JSON results to the file
        with open(output_path, "w", encoding="utf-8") as out_fh:
            out_fh.write(json_text)

        # Update the log file, so we won't refetch this item
        with open(args.logfile, "a", encoding="utf-8") as log_fh:
            log_fh.write(short_inchikey + "\n")

        # Add to our local set so we skip it in the same run
        processed_keys.add(short_inchikey)
        files_saved += 1

        # Respect the rate limit (e.g., 10 requests / second => 0.1s pause)
        time.sleep(1.0 / args.rate_limit)

if __name__ == "__main__":
    main()
