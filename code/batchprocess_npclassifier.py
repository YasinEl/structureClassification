#!/usr/bin/env python3

import argparse
import csv
import os
import time
import requests
from tqdm import tqdm
from urllib.parse import quote

# RDKit imports
from rdkit import Chem
from rdkit.Chem import inchi

NPCLASSIFIER_URL = "https://npclassifier.gnps2.org/classify?smiles="

def smiles_to_inchikey(smiles_str):
    """
    Converts a SMILES string to an InChIKey using RDKit.
    Returns None if SMILES is invalid or conversion fails.
    """
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    return inchi.MolToInchiKey(mol)

def get_npclassifier_json(smiles, max_retries=3):
    """
    Retrieves JSON classification from npclassifier for the given SMILES string.
    Retries up to max_retries times on failure, waiting 3s between attempts.
    
    :param smiles: A valid SMILES string
    :param max_retries: Maximum retry attempts
    :return: JSON string from npclassifier
    :raises requests.RequestException: if all retries fail
    """

    # URL-encode the SMILES
    smiles_encoded = quote(smiles, safe="")  # safe="" ensures everything is encoded

    attempt = 0
    while attempt < max_retries:
        try:
            response = requests.get(f"{NPCLASSIFIER_URL}{smiles_encoded}")
            response.raise_for_status()
            # npclassifier returns JSON, so return response as text
            return response.text
        except requests.RequestException as e:
            attempt += 1
            if attempt < max_retries:
                print(f"[WARNING] Request failed ({e}). "
                      f"Retry {attempt}/{max_retries} in 3s...")
                time.sleep(3)
            else:
                raise

def main():
    parser = argparse.ArgumentParser(
        description="Fetch npclassifier JSON for a list of SMILES in a CSV file, "
                    "logging and storing them by InChIKey prefix."
    )
    parser.add_argument(
        "input_csv",
        help="Path to the input CSV file (must contain a column named 'Smiles')."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to store downloaded JSON files. Subfolders are created as needed."
    )
    parser.add_argument(
        "--logfile",
        required=True,
        help="Path to a log file where processed short InChIKeys are stored."
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=10.0,
        help="Max requests per second (default: 10)."
    )
    args = parser.parse_args()

    delay_between_requests = 1.0 / args.rate_limit if args.rate_limit > 0 else 0

    # 1) Load previously processed InChIKeys from the log file
    processed_keys = set()
    if os.path.exists(args.logfile):
        with open(args.logfile, "r", encoding="utf-8") as log_fh:
            for line in log_fh:
                cleaned = line.strip()
                if cleaned:
                    processed_keys.add(cleaned)

    # Keep track of how many files we've saved to determine subfolder indexing
    files_saved = len(processed_keys)

    # Ensure the output directory exists
    os.makedirs(args.outdir, exist_ok=True)

    # 2) Read entire TSV so we know the length for our progress bar
    with open(args.input_csv, "r", encoding="utf-8") as csv_fh:
        reader = list(csv.DictReader(csv_fh, delimiter=","))

    # Set up progress bar
    progress_bar = tqdm(total=len(reader), desc="Processing SMILES", unit="rows")

    # 3) Process each row
    for row in reader:
        progress_bar.update(1)

        if "Smiles" not in row:
            progress_bar.write("[ERROR] No 'Smiles' column in this row. Skipping.")
            continue

        smiles_str = row["Smiles"].strip()
        if not smiles_str:
            progress_bar.write("[INFO] Blank SMILES. Skipping.")
            continue

        # check if SMILES is valid
        try:
            mol = Chem.MolFromSmiles(smiles_str)
        except Exception as e:
            progress_bar.write(f"[ERROR] Invalid SMILES '{smiles_str}': {e}. Skipping.")
            continue

        # Convert SMILES to InChIKey
        inchikey = smiles_to_inchikey(smiles_str)
        if not inchikey:
            progress_bar.write(f"[ERROR] Invalid SMILES, cannot convert to InChIKey: '{smiles_str}'. Skipping.")
            continue

        # Get the short prefix from the InChIKey
        short_inchikey = inchikey.split("-", 1)[0]

        # Skip if we already processed that short prefix
        if short_inchikey in processed_keys:
            progress_bar.write(f"[INFO] Short key '{short_inchikey}' already processed. Skipping.")
            continue

        # Attempt to retrieve JSON classification from npclassifier
        try:
            json_text = get_npclassifier_json(smiles_str, max_retries=3)
        except requests.RequestException as e:
            progress_bar.write(f"[ERROR] Could not fetch data for SMILES '{smiles_str}': {e}")
            continue

        # Determine subfolder: each folder holds up to 100 files
        folder_index = files_saved // 100
        subfolder_path = os.path.join(args.outdir, str(folder_index))
        os.makedirs(subfolder_path, exist_ok=True)

        # Save JSON using short InChIKey as filename
        output_filename = f"{short_inchikey}.json"
        output_path = os.path.join(subfolder_path, output_filename)
        with open(output_path, "w", encoding="utf-8") as out_fh:
            out_fh.write(json_text)

        # Log success
        with open(args.logfile, "a", encoding="utf-8") as log_fh:
            log_fh.write(short_inchikey + "\n")

        processed_keys.add(short_inchikey)
        files_saved += 1

        progress_bar.write(f"[SUCCESS] Fetched & saved SMILES -> {output_path}")

        # Respect rate limit
        time.sleep(delay_between_requests)

    progress_bar.close()
    print("All done!")

if __name__ == "__main__":
    main()
