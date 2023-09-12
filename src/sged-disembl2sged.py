#! /usr/bin/python

""" Created on 30/07/20 by jdutheil

    Convert DisEMBL scores into a SGED file.
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

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "d:o:ch"
full_opt = ["disembl=", "output=", "csv", "help"]

def usage() :
    print(
"""
sged-disembl2sged

    Convert the output of the DISEMBL program to predict intrinsically disordered regions
    into SGED format.

Available arguments:
    --disembl (-d): Input DISEMBL file (required).
    --output (-o): Output SGED file (required).
    --csv (-c): Input SGED file is with comas instead of tabs (default: tabs).
    --help (-h): Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = "\t"  # TSV by default
group_col = "Group"
for arg, val in arguments:
    if arg in ("-d", "--disembl"):
        disembl_file = val
        print("DisEMBL score file: %s" % disembl_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)
    elif arg in ("-c", "--csv"):
        tabsep = ","
    elif arg in ("-h", "--help"):
        usage()

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","

# Check required arguments

if not 'disembl_file' in globals():
    print("Error: a DISEMBL input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()

# Start parsing:

with open(disembl_file, "r") as scores:
    begin = False
    for line in scores:
        if line.startswith("# RESIDUE"):
            begin = True
            break
    with open(output_file, "w") as output:
        output.write(tabsep.join(["Group", "AA", "COILS", "REM465", "HOTLOOP"]))
        output.write("\n")
        i = 0
        for line in scores:
            i = i + 1
            line = line.strip()
            cols = line.split("\t")
            output.write(tabsep.join(["[%s]" % i, cols[0], cols[1], cols[2], cols[3]]))
            output.write("\n")

print("Done.")
