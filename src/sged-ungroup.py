#!/usr/bin/env python3

""" Created on 13/02/20 by jdutheil
    Modified on 30/07/2025 by lorenzopenone

    Convert multi-sites groups into single sites groups. 
    Allow to specify which column to replicate.
    This program is part of the SgedTools package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:o:d:g:ch"  # -m (match) removed
full_opt = [
    "sged=",
    "output=",
    "data=",
    "group=",
    "csv",
    "help"
]

def usage():
    print(
        """
sged-ungroup

    Convert multi-sites groups into single-site groups. 
    Allow to specify which column(s) to replicate.
    You can now split several columns in parallel using --group / -g with a
    comma-separated list.

Available arguments:
    --sged (-s):   Input SGED file (required).
    --output (-o): Output SGED file (required).
    --group (-g):  Column **or comma-separated list of columns** where group
                   coordinates are stored (default: Group).
    --data (-d):   Column selection (default: empty). Columns to copy to output.
                   Their values are duplicated for each emitted row.
    --csv (-c):    Input SGED uses commas instead of tabs (default: tabs).
    --help (-h):   Print this message.
"""
    )
    sys.exit()

try:
    arguments, _ = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

# Defaults
use_tsv = True  # TSV default
selected_cols = []
group_cols = ["Group"]  # may become several columns

for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print(f"SGED file: {sged_file}")
    elif arg in ("-o", "--output"):
        output_file = val
        print(f"Output ungrouped file: {output_file}")
    elif arg in ("-d", "--data"):
        selected_cols = val.split(",") if val else []
    elif arg in ("-g", "--group"):
        group_cols = [c.strip() for c in val.split(",") if c.strip()]
        print("Group coordinates are in column(s): %s" % ", ".join(group_cols))
    elif arg in ("-c", "--csv"):
        use_tsv = False
    elif arg in ("-h", "--help"):
        usage()

# Mandatory arguments check
if 'sged_file' not in globals() or 'output_file' not in globals():
    usage()

delim = "\t" if use_tsv else ","
print("SGED file is in %s format" % ("TSV" if use_tsv else "CSV"))

with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep=delim, dtype=str, comment="#")

    # Check that requested columns exist
    missing = [c for c in group_cols if c not in df.columns]
    if missing:
        print(f"ERROR: group column(s) not found: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    header_cols = group_cols + (selected_cols if selected_cols else [])

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write(delim.join(header_cols) + "\n")

        for i, g in enumerate(df[group_cols[0]]):  # primary split column drives iteration
            if not isinstance(g, str):
                continue
            tmp = g.strip()[1:-1].replace(" ", "")  # remove [ ] and spaces
            if tmp == "":
                continue
            positions = tmp.split(";")

            # Build dictionaries for all group columns
            split_values = {}
            for col in group_cols:
                raw_val = df[col].iloc[i]
                if not isinstance(raw_val, str):
                    split_values[col] = [str(raw_val)] * len(positions)
                    continue
                rv = raw_val.strip()
                if rv.startswith("[") and rv.endswith("]"):
                    inner = rv[1:-1].replace(" ", "")
                    vals = inner.split(";") if inner else [""]
                else:
                    vals = [rv]

                # Harmonise length by replication when needed
                if len(vals) == 1 and len(positions) > 1:
                    vals *= len(positions)
                if len(vals) != len(positions):
                    print(
                        f"WARNING: length mismatch in row {i}: column '{col}' has {len(vals)} items while '{group_cols[0]}' has {len(positions)}. Replicating first value.",
                        file=sys.stderr,
                    )
                    vals = [vals[0]] * len(positions)
                split_values[col] = vals

            for idx, pos in enumerate(positions):
                # Bracket each split value
                split_out = [f"[{split_values[col][idx]}]" for col in group_cols]

                if selected_cols:
                    row_vals = df[selected_cols].iloc[i].astype(str).tolist()
                    handle.write("%s%s%s\n" % (delim.join(split_out), delim, delim.join(row_vals)))
                else:
                    handle.write("%s\n" % delim.join(split_out))

print("Done.")