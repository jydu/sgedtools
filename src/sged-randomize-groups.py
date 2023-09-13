#! /usr/bin/python

""" Created on 12/03/20 by jdutheil

    Generates a list of groups with characteristics similar to a given set of groups.
    The output groups have the same size and similar site properties, but sites are taken randomly.
    To remove bias due to the skewed distribution of norms, a correction is further introduced.
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


import getopt, sys, os.path, warnings
import pandas
import numpy
from progress.bar import Bar

# For convenience
def ifelse_fun(test, yes, no):
    if test:
        res = yes
    else:
        res = no
    return res


cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:r:o:m:g:k:cn:t:bh"
full_opt = [
    "sged-groups=",
    "sged-sites=",
    "output=",
    "measure=",
    "groups=",
    "sites=",
    "csv",
    "number-replicates=",
    "similarity-threshold=",
    "minimum-observations=",
    "help"
]

def usage() :
    print(
"""
sged-randomize-groups

    Generates a list of groups with characteristics similar to a given set of groups.
    The output groups have the same size and similar site properties, but sites are taken randomly.
    If no conditional variable is specified, all sites are sampled independently of their characteristics.
    To remove bias due to skewed distribution of conditional variables, a correction is further introduced.

Available arguments:
    --sged-groups (-s): SGED file with groups to be tested (required).
    --sged-sites (-r): SGED file with sites to be sampled (required).
    --group (-g): Column where group coordinates are read (default: Group).
    --site (-k): Column where site coordinates are read (default: Group).
    --output (-o): Output SGED file (required).
    --measure (-m): Column where the conditional variable should be read (default: none).
    --number-replicates (-n): How many resampling of each group should be performed (default: 100).
    --similarity-threshold (-t): Similarity threshold (%) for the considering sites with similar property.
        For instance, if N is the property considered, given my the --measure option,
        a sampled site will be considered similar to the observed site with a 10% threshold if
        (Nsample - Nobs) / Nobs < 0.1.
    --minimum-observations (-b): minimum number of sample-able sites required for an observed site to be resampled.
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

tabsep = True  # TSV by default
group_col = "Group"
site_col = "Group"
cond_var = ""
n_rep = 100
# Similarity threshold used for conditional variable:
sim_t = 0.1  # 10% difference: (Nsim - Nobs) / Nobs < 0.1
min_obs = 5  # Minimum number of matching site. If not enough sites are found for a given position, a warning will be issued. You may try to increase the threshold if too many warnings are produced.


for arg, val in arguments:
    if arg in ("-s", "--sged-groups"):
        sged_file_groups = val
        print("SGED file for groups: %s" % sged_file_groups)
    elif arg in ("-r", "--sged-sites"):
        sged_file_sites = val
        print("SGED file for sites: %s" % sged_file_sites)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output info file: %s" % output_file)
    elif arg in ("-m", "--measure"):
        cond_var = val
    elif arg in ("-g", "--group"):
        group_col = val
        print("Group coordinates are in column: %s" % group_col)
    elif arg in ("-k", "--site"):
        site_col = val
        print("Site coordinates are in column: %s" % site_col)
    elif arg in ("-n", "--number-replicates"):
        n_rep = int(val)
        if n_rep > 4294967295:
            print(
                "Number of replicates too high. Consider using uint64 instead of uint32 if needed."
            )
            sys.exit(2)
    elif arg in ("-t", "--similarity-threshold"):
        sim_t = float(val)
    elif arg in ("-b", "--minimum-observations"):
        min_obs = int(val)
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

if not 'sged_file_groups' in globals():
    print("Error: a group SGED input file should be specified.")
    usage()
if not 'sged_file_sites' in globals():
    print("Error: a site SGED input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()

print("Number of replicates per group: %i" % n_rep)
if cond_var != "":
    print("Similarity threshold for conditional variable: %f" % sim_t)

###########################################


# Read input data:
with open(sged_file_sites) as csv_file_sites:
    df_sites = pandas.read_csv(
        csv_file_sites, sep = delim, comment = "#", na_values = "NA", keep_default_na = False
    )  # NA in columns ignored
    groups = df_sites[site_col]
    groups = [g[1 : (len(g) - 1)].replace(" ", "") for g in groups]
    df_sites["Group"] = groups
    df_sites.set_index(
        "Group", drop=False, inplace=True
    )  # We use the site (as string, not number) as an index

with open(sged_file_groups) as csv_file_groups:
    df_groups = pandas.read_csv(
        csv_file_groups, sep = delim, comment = "#", na_values = "NA", keep_default_na = False
    )
    groups = df_groups[group_col]
    groups_lst = [ g[1 : (len(g) - 1)].split(";") for g in groups ]

    # Compute group sizes
    groups_sizes = [ len(g) for g in groups_lst ]


n_groups = len(df_groups.index)

# Now replicate each group:
x_rep = numpy.zeros(n_groups * n_rep, dtype = numpy.uint32)
l_grp = [[]] * (n_groups * n_rep)
x_grp = [""] * (n_groups * n_rep)
x_ave = numpy.zeros(n_groups * n_rep)  # Average of the sampled group
x_siz = numpy.zeros(n_groups * n_rep, dtype = numpy.uint32)
x_oav = numpy.zeros(n_groups * n_rep)  # Average the original group

i = 0
pbar = Bar("Simulating...", max = n_groups)
for grp in range(n_groups):
    size = groups_sizes[grp]
    print("Simulating for group %i with size %i\n" % (grp + 1, size))

    # Get all sites with adequate value for each position:
    gp = groups_lst[grp]

    gp_vals = [1] * size # If no conditional variable, we set a value of 1 for all sites 
    if cond_var != "":
        gp_vals = [df_sites.loc[gp[j], cond_var] for j in range(size)]
    x_rep[i : (i + n_rep)] = [x for x in range(1, n_rep + 1)]
    if all(~numpy.isnan(gp_vals)):

        x_siz[i : (i + n_rep)] = [size] * n_rep
        x_grp[i : (i + n_rep)] = ["["] * n_rep
        x_ave[i : (i + n_rep)] = [0] * n_rep
        x_oav[i : (i + n_rep)] = [numpy.mean(gp_vals)] * n_rep

        # Loop over each site in the group:
        for sit in range(size):
            if cond_var != "":
                # Get all sites with similar rate (or any other conditional variable):
                t = [abs(v - gp_vals[sit]) / gp_vals[sit] for v in df_sites[cond_var]]
                cond_sites = df_sites[[x <= sim_t for x in t]]
            else:
                cond_sites = df_sites #All siters are used.

            # Loop over each simulation replicate:
            for sim in range(n_rep):
                # Remove sites already present in the group
                if i + sim < len(l_grp):  # test if element already exists in list
                    tmp = cond_sites[~cond_sites[site_col].isin(l_grp[i + sim])]
                else:
                    tmp = cond_sites

                if cond_var != "":
                    # Sampling correction:
                    tmpl = tmp[tmp[cond_var] < gp_vals[sit]]
                    tmpg = tmp[tmp[cond_var] > gp_vals[sit]]
                    tmpe = tmp[tmp[cond_var] == gp_vals[sit]]
                    n = min(len(tmpl.index), len(tmpg.index))
                    n = max(
                        n, min_obs
                    )  # we add min_obs there to avoid getting no replicate site when we have extreme values for the candidate site.

                    cat_list = [
                        tmpl.sample(min(n, len(tmpl.index))),
                        tmpe,
                        tmpg.sample(min(n, len(tmpg.index))),
                    ]
                    tmp2 = pandas.concat(cat_list)
                else:
                    tmp2 = tmp

                if len(tmp2.index) == 0:
                    warnings.warn(
                        "No more similar site available for candidate site %i in group %i replicate %i"
                        % (sit, grp, sim)
                    )
                    x_grp[i + sim] = "%s%sNA" % (
                        x_grp[i + sim],
                        ifelse_fun(x_grp[i + sim] == "[", "", ";")
                    )
                    x_ave[i + sim] = numpy.nan
                else:
                    if len(tmp2.index) < min_obs:
                        warnings.warn(
                            "Minimum site frequency not matched for candidate site %i in group %i replicate %i"
                            % (sit, grp, sim)
                        )

                    # Now sample sites:
                    x = tmp2.sample(1)
                    l_grp[i + sim] = l_grp[i + sim] + x["Group"].to_numpy().tolist()

                    x_grp[i + sim] = "%s%s%s" % (
                        x_grp[i + sim],
                        ifelse_fun(x_grp[i + sim] == "[", "", ";"),
                        x["Group"].iloc[0]
                    )
                    if cond_var != "":
                        x_ave[i + sim] = x_ave[i + sim] + float(x[cond_var].iloc[0])

        for j in range(i, i + n_rep):
            x_grp[j] = x_grp[j] + "]"

    i = i + n_rep
    pbar.next()

if cond_var != "":
    x_ave = x_ave / x_siz
    results = {
        "Replicate": x_rep,
        "Group": x_grp,
        "Size": x_siz,
        "RandMean": x_ave,
        "OrigMean": x_oav,
    }
    df = pandas.DataFrame(
        results, columns=["Replicate", "Group", "Size", "RandMean", "OrigMean"]
    )
else:
    x_ave = x_ave / x_siz
    results = {
        "Replicate": x_rep,
        "Group": x_grp,
        "Size": x_siz
    }
    df = pandas.DataFrame(
        results, columns=["Replicate", "Group", "Size"]
    )

# Write results:
df.to_csv(output_file, sep = delim, na_rep = "NA", index = False)

print("\nDone.")
