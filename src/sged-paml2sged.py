#!/usr/bin/env python3

import pandas as pd
import getopt, sys

""" Created on 03/03/23 by dyrhmd

    Converts the output of PAML positive selection test to SGED.
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

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "p:o:m:t:ch"
full_opt = ["paml=", "output=", "method=", "threshold=", "csv", "help"]

def usage() :
    print(
"""
sged-paml2sged

    Converts the output of PAML positive selection test to SGED.

Available arguments:
    --paml (-p): PAML output file (required).
    --output (-o): Output SGED file (required).
    --method (-m): Method to consider, 'bayesian' or 'naive' (default: bayesian).
    --threshold (-t): Minimum posterior probability to include a site (default: 0).
    --csv (-c): Input SGED file is with comas instead of tabs (default)
    --help (-h): Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

method = "bayesian"
tabseq = True
min_post_prob = 0.
for arg, val in arguments:
    if arg in ("-p", "--paml"):
        paml_file = val
        print("PAML result file: %s" % paml_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)
    elif arg in ("-m", "--method"):
        method = val.lower()
        print("Method: %s" % method)
    elif arg in ("-t", "--threshold"):
        min_post_prob = float(val)
        print("Minimum posterior probability: %s" % min_post_prob)
    elif arg in ("-c", "--csv"):
        tabseq = False
    elif arg in ("-h", "--help"):
        usage()
        
if tabseq:
    print("SGED file is in TSV format.")
    delim = "\t"
else:
    print("SGED file is in CSV format.")
    delim = ","

# Check required arguments

if not 'paml_file' in globals():
    print("Error: a PAML input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()

# Read the paml result file
with open(paml_file, 'r') as f:
    contents = f.readlines()
    
if method == "bayesian":
    #find the positive selected sites under the BEB
    
    start_row = [i for i, line in enumerate(contents) if 'Bayes Empirical Bayes (BEB)' in line][0] + 2
    end_row = [i for i, line in enumerate(contents) if 'The grid' in line][0] - 2

elif method == "naive":
    #find the positive selected sites under the NEB
    
    start_row = [i for i, line in enumerate(contents) if 'Naive Empirical Bayes (NEB)' in line][0] + 3
    end_row = [i for i, line in enumerate(contents) if 'Bayes Empirical Bayes (BEB)' in line][0] - 1

else:
    print("Error: Method is not valid.")
    sys.exit(2)

lines = contents[start_row:end_row]

# Extract the information
positive_sites = []
for line in lines:
    site_info = line.split()
    positive_sites.append(site_info)

# convert it to the data frame and add square brackets
df = pd.DataFrame(positive_sites, columns = ['position', 'amino_acid', 'probability'])
df['probability'] = df['probability'].astype(float)
df.insert(loc = 0, column = 'Group', value = '[' + df['position'] + ']')
df.drop(['position'], axis = 1, inplace = True)

# filter the results:
df = df[(df['probability'] >= min_post_prob)]

# converting to the csv file
df.to_csv(output_file, index = False, sep = delim)

print("Done.")
