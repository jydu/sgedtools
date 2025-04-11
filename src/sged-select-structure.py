#!/usr/bin/env python3
"""
sged-select-best-structure.py

Created on 11/04/25 by rdurak

Processes a list of protein structures (mmCIF or PDB format),
extracts structural metadata, and outputs:

- Resolution (in Å)
- Experimental method (X-ray, Cryo-EM, NMR, etc.)
- Structural completeness (number of residues modeled)
- Deposition date
- Suggestion for the protein structure

It allows filtering by resolution threshold, and optionally by experimental method.
The best structure is suggested based on resolution, completeness, and recency.

Outputs:
- A CSV summary of selected structures
- A simple .txt list of accepted CIF/PDB files (no headers)
- Printed suggestion of the best structure (if any)

Usage:
    python sged-select-best-structure.py --structures *.cif --format mmCif \
        --resolution-threshold 5.0 --csv-report structures_summary.csv \
        --preferred-method "X-RAY DIFFRACTION"
"""

import os
import argparse
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from datetime import datetime

def get_structure_info(file, format):
    if format == 'mmCif':
        parser = MMCIFParser(QUIET=True)
    elif format == 'PDB':
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError("Unsupported format: must be 'PDB' or 'mmCif'")

    try:
        structure = parser.get_structure("protein", file)

        mmcif_dict = parser._mmcif_dict if format == 'mmCif' else {}

        # Resolution
        resolution = mmcif_dict.get('_refine.ls_d_res_high')
        if isinstance(resolution, list):
            resolution = resolution[0] if resolution else None
        resolution = float(resolution) if resolution else None

        # Experimental method
        method = mmcif_dict.get('_exptl.method', 'UNKNOWN')
        if isinstance(method, list):
            method = method[0] if method else 'UNKNOWN'
        method = method.strip().upper() if method else 'UNKNOWN'

        # Deposition date
        date = mmcif_dict.get('_pdbx_database_status.recvd_initial_deposition_date')
        if isinstance(date, list):
            date = date[0] if date else None
        deposition_date = None
        try:
            deposition_date = datetime.strptime(date, "%Y-%m-%d") if date else None
        except Exception:
            pass

        # Completeness: number of residues
        residue_count = 0
        for model in structure:
            for chain in model:
                residue_count += len([res for res in chain if res.id[0] == ' '])
            break  # Use only model 0

        return {
            'file': file,
            'resolution': resolution,
            'method': method,
            'completeness': residue_count,
            'deposition_date': deposition_date.strftime("%Y-%m-%d") if deposition_date else None
        }

    except Exception as e:
        print(f"[ERROR] Could not parse {file}: {e}")
        return None

def suggest_best_structure(infos):
    valid = [s for s in infos if s['resolution'] is not None and s['completeness'] is not None]
    if not valid:
        return None

    def score(s):
        date_score = datetime.strptime(s['deposition_date'], "%Y-%m-%d") if s['deposition_date'] else datetime.min
        return (s['resolution'], -s['completeness'], -date_score.timestamp())

    return sorted(valid, key=score)[0]

def main():
    parser = argparse.ArgumentParser(description="Select and summarize structures based on resolution and completeness.")
    parser.add_argument("--structures", "-s", nargs="+", required=True, help="Input structure files (.cif or .pdb)")
    parser.add_argument("--format", "-f", default="mmCif", choices=["mmCif", "PDB"], help="File format")
    parser.add_argument("--resolution-threshold", "-r", type=float, default=8.0, help="Resolution cutoff (e.g. 5.0)")
    parser.add_argument("--preferred-method", "-m", default=None, help="Optional: restrict to a single experimental method")
    parser.add_argument("--csv-report", "-c", help="Output CSV summary report")
    parser.add_argument("--list-output", "-l", default="passed_structures.txt", help="Plain text output of accepted file names")

    args = parser.parse_args()

    print(f"Scanning {len(args.structures)} structure files...\n")
    infos = [get_structure_info(file, args.format) for file in args.structures]
    infos = [s for s in infos if s]

    # Optional filter by method
    if args.preferred_method:
        method = args.preferred_method.strip().upper()
        infos = [s for s in infos if s['method'] == method]

    # Filter by resolution
    filtered = [s for s in infos if s['resolution'] is not None and s['resolution'] <= args.resolution_threshold]

    best = suggest_best_structure(filtered)
    for s in filtered:
        s["suggested"] = "yes" if s == best else "no"

    if args.csv_report:
        df = pd.DataFrame(filtered)
        df.to_csv(args.csv_report, index=False)
        print(f"Summary CSV saved to: {args.csv_report}")

    # Save plain .txt list of files that passed
    if filtered:
        with open(args.list_output, 'w') as out:
            out.write("\n".join(s['file'] for s in filtered))
        print(f"List of accepted structures written to: {args.list_output}")

    if best:
        print(f"\nSuggested Best Structure:")
        print(f"  File: {best['file']}")
        print(f"  Resolution: {best['resolution']} Å")
        print(f"  Method: {best['method']}")
        print(f"  Completeness: {best['completeness']} residues")
        print(f"  Deposition Date: {best['deposition_date']}")
    else:
        print("No suitable structure passed the given threshold.")

if __name__ == "__main__":
    main()
