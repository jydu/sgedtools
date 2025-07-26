#!/usr/bin/env python3

""" Created on 13/02/20 by jdutheil
    Modified on 26/07/2025 by lorenzopenone

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

unix_opt = "s:o:d:g:cm:h"
full_opt = ["sged=", "output=", "data=", "group=", "csv", "match=", "help"]

def usage() :
    print(
"""
sged-ungroup

    Convert multi-sites groups into single sites groups. 
    Allow to specify which column to replicate.
    Optionally split a matched column in parallel.

Available arguments:
    --sged (-s):   Input SGED file (required).
    --output (-o): Output SGED file (required).
    --group (-g):  Column where group coordinates are stored (default: Group).
    --data (-d):   Column selection (default: empty). Columns to copy to output.
                   Their values are duplicated for each emitted row.
    --match (-m):  Optional column to split in parallel with --group.
                   If provided, it will be element-wise matched.
    --csv (-c):    Input SGED uses commas instead of tabs (default: tabs).
    --help (-h):   Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
selected_cols = []
group_col = "Group"
match_col = None

for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output ungrouped file: %s" % output_file)
    elif arg in ("-d", "--data"):
        selected_cols = val.split(",") if val else []
    elif arg in ("-g", "--group"):
        group_col = val
        print("Group coordinates are in column: %s" % group_col)
    elif arg in ("-m", "--match"):
        match_col = val
        print("Matched column: %s" % match_col)
    elif arg in ("-c", "--csv"):
        tabsep = False
    elif arg in ("-h", "--help"):
        usage()

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","

# Check required arguments
if not 'sged_file' in globals():
    usage()
if not 'output_file' in globals():
    usage()

# Ensure the matched column is part of the output if provided
if match_col and match_col not in selected_cols:
    selected_cols.append(match_col)

# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep=delim, dtype=str, comment='#')

    # Validate presence of required columns
    if group_col not in df.columns:
        print(f"ERROR: group column '{group_col}' not found in input.", file=sys.stderr)
        sys.exit(1)
    if match_col and match_col not in df.columns:
        print(f"ERROR: match column '{match_col}' not found in input.", file=sys.stderr)
        sys.exit(1)

    groups = df[group_col]

    with open(output_file, "w") as handle:
        # Header
        if selected_cols:
            handle.write("Group%s%s\n" % (delim, delim.join(df[selected_cols].columns)))
        else:
            handle.write("Group\n")

        # Rows
        for i, g in enumerate(groups):
            if not isinstance(g, str):
                continue
            tmp = g.strip()[1 : (len(g.strip()) - 1)]  # remove [ ]
            tmp = tmp.replace(" ", "")
            if tmp == "":
                continue
            positions = tmp.split(";")

            # Prepare match values if requested
            match_values = None
            if match_col:
                mv_raw = df[match_col].iloc[i]
                if isinstance(mv_raw, str) and mv_raw.startswith("[") and mv_raw.endswith("]"):
                    mv = mv_raw.strip()[1 : (len(mv_raw.strip()) - 1)].replace(" ", "")
                    match_values = mv.split(";") if mv != "" else []
                else:
                    # single value, replicate
                    match_values = [mv_raw] * len(positions)

                if len(match_values) != len(positions):
                    print(
                        f"WARNING: length mismatch in row {i}: "
                        f"{group_col} has {len(positions)} items, "
                        f"{match_col} has {len(match_values)} items. "
                        f"Replicating unchanged {match_col}.",
                        file=sys.stderr
                    )
                    # fall back to replicate the original field value
                    mv_raw = df[match_col].iloc[i]
                    match_values = [mv_raw] * len(positions)

            for idx, j in enumerate(positions):
                # Build the list of output values for selected columns
                if selected_cols:
                    row_vals = df[selected_cols].iloc[i].astype(str).tolist()
                    # If match_col is present, replace its value with the matched element
                    if match_col:
                        try:
                            mpos = selected_cols.index(match_col)
                            row_vals[mpos] = match_values[idx]
                        except Exception:
                            pass
                    handle.write("[%s]%s%s\n" % (j, delim, delim.join(row_vals)))
                else:
                    handle.write("[%s]\n" % j)

print("Done.")