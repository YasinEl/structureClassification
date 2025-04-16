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

NAT_PROD_ENDPOINT = (
    "https://api.naturalproducts.net/latest/chem/descriptors"
    "?smiles={smiles_encoded}&format=json&toolkit=rdkit"
)

def smiles_to_inchikey(smiles_str):
    """
    Converts a SMILES string to an InChIKey using RDKit.
    Returns None if SMILES is invalid or conversion fails.
    """
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    return inchi.MolToInchiKey(mol)

def get_natural_products_json(smiles, max_retries=3):
    """
    Retrieves JSON descriptor data from the Natural Products API for the given SMILES string.
    Retries up to max_retries times on failure, waiting 3s between attempts.

    :param smiles: A valid SMILES string
    :param max_retries: Maximum retry attempts
    :return: JSON string from the Natural Products endpoint
    :raises requests.RequestException: if all retries fail
    """
    # URL-encode the SMILES
    smiles_encoded = quote(smiles, safe="")  # safe="" ensures everything is encoded

    attempt = 0
    while attempt < max_retries:
        try:
            url = NAT_PROD_ENDPOINT.format(smiles_encoded=smiles_encoded)
            response = requests.get(url)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            attempt += 1
            if attempt < max_retries:
                print(f"[WARNING] Request failed for SMILES '{smiles}' ({e}). "
                      f"Retry {attempt}/{max_retries} in 3s...")
                time.sleep(3)
            else:
                raise

def main():
    parser = argparse.ArgumentParser(
        description="Fetch JSON from the Natural Products descriptors API for SMILES in CSV. "
                    "Output is saved/logged by the short InChIKey prefix."
    )
    parser.add_argument(
        "input_csv",
        help="Path to the input CSV (must contain a column named 'Smiles')."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to store JSON files. Subfolders created (max 100 files per folder)."
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

    # Convert rate-limit to a delay between requests
    delay_between_requests = 1.0 / args.rate_limit if args.rate_limit > 0 else 0

    # 1) Load already processed short InChIKeys
    processed_keys = set()
    if os.path.exists(args.logfile):
        with open(args.logfile, "r", encoding="utf-8") as log_fh:
            for line in log_fh:
                cleaned = line.strip()
                if cleaned:
                    processed_keys.add(cleaned)

    files_saved = len(processed_keys)
    os.makedirs(args.outdir, exist_ok=True)

    # 2) Read entire CSV for progress bar
    with open(args.input_csv, "r", encoding="utf-8") as csv_fh:
        reader = list(csv.DictReader(csv_fh, delimiter=","))

    progress_bar = tqdm(total=len(reader), desc="Processing rows", unit="rows")

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

        # Convert SMILES -> InChIKey
        inchikey = smiles_to_inchikey(smiles_str)
        if inchikey is None:
            progress_bar.write(f"[ERROR] Invalid SMILES '{smiles_str}'. Cannot convert to InChIKey. Skipping.")
            continue

        short_inchikey = inchikey.split("-", 1)[0]
        if short_inchikey in processed_keys:
            progress_bar.write(f"[INFO] '{short_inchikey}' already processed. Skipping.")
            continue

        # Try to fetch JSON
        try:
            json_text = get_natural_products_json(smiles_str, max_retries=3)
        except requests.RequestException as e:
            progress_bar.write(f"[ERROR] Could not fetch data for SMILES '{smiles_str}': {e}")
            continue

        # Organize into subfolders
        folder_index = files_saved // 100
        subfolder_path = os.path.join(args.outdir, str(folder_index))
        os.makedirs(subfolder_path, exist_ok=True)

        # Write JSON
        output_filename = f"{short_inchikey}.json"
        output_path = os.path.join(subfolder_path, output_filename)
        with open(output_path, "w", encoding="utf-8") as out_fh:
            out_fh.write(json_text)

        # Update log, sets
        with open(args.logfile, "a", encoding="utf-8") as log_fh:
            log_fh.write(short_inchikey + "\n")

        processed_keys.add(short_inchikey)
        files_saved += 1
        progress_bar.write(f"[SUCCESS] Saved '{short_inchikey}' -> {output_path}")

        # Respect rate limit
        time.sleep(delay_between_requests)

    progress_bar.close()
    print("All done!")

if __name__ == "__main__":
    main()
