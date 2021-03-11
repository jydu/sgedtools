#! /usr/bin/python

""" Created on 12/03/20 by jdutheil

    Generates a list of groups with characteristics similar to a given set of groups.
    The output groups have the same size and similar site porperties, but sites are taken randomly.
    In this version, site norms are not discretized but a similarity threshold is used.
    To remove bias due to skewed distribution of norms, a correction is further introduced.
"""


import getopt, sys, os.path
import pandas
import numpy

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:r:p:o:m:g:h:cn:t:b"
full_opt = ["sged-groups=", "sged-sites=", "pdb=", "output=", "measure=", "groups=", "sites=", "csv", "number-replicates=", "similarity-threshold=", "minimum-observations="]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

tabsep = True # TSV by default
group_col = "Group"
site_col = "Group"
cond_var = ""
n_rep = 100
# Similarity threshold used for conditional variable:
sim_t = 0.1 # 10% difference: (Nsim - Nobs) / Nobs < 0.1
min_obs = 5 # Minimum number of matching site. If not enough sites are found for a given position, a warning will be issued. You may try to increase the threshold if too many warnings are produced. 


for arg, val in arguments:
  if arg in ("-s", "--sged-groups"):
    sged_file_groups = val
    print "SGED file for groups: %s" % sged_file_groups
  elif arg in ("-r", "--sged-sites"):
    sged_file_sites = val
    print "SGED file for sites: %s" % sged_file_sites
  elif arg in ("-p", "--pdb"):
    pdb_file = val
    print "PDB file: %s" % pdb_file
  elif arg in ("-o", "--output"):
    output_file = val
    print "Output info file: %s" % output_file
  elif arg in ("-m", "--measure"):
    cond_var = val
  elif arg in ("-g", "--groups"):
    group_col = val
    print "Group coordinates are in column: %s" % group_col
  elif arg in ("-h", "--sites"):
    site_col = val
    print "Site coordinates are in column: %s" % site_col
  elif arg in ("-n", "--number-replicates"):
    n_rep = int(val)
  elif arg in ("-t", "--similarity-threshold"):
    sim_t = float(val)
  elif arg in ("-b", "--minimum-observations"):
    min_obs = int(val)
  elif arg in ("-c", "--csv"):
    tabsep = False

if tabsep:
  print "SGED file is in TSV format"
  delim = '\t'
else:
  print "SGED file is in CSV format"
  delim = ','


print "Number of replicates per group: %i" % n_rep
if cond_var != "":
  print "Similarity threshold for conditional variable: %f" % sim_t

# Input:
#sitesPath <- "Myo_sites.csv"
#groupsPath <- "Myo_stats_pvalues.csv"
# Output path:
#outputPath <- "Myo_random.csv"



###########################################


# Read input data:
with open(sged_file_sites) as csv_file_sites:
  df_sites = pandas.read_csv(csv_file_sites, sep = delim, dtype = str, keep_default_na = False) #NA in columns ignored
  groups = df_sites[site_col]
  groups = [g[1:(len(g)-1)].replace(' ', '') for g in groups]
  df_sites["Group"] = groups
  df_sites.set_index("Group", drop = True, inplace = True) # We use the site (as string, not number) as an index

with open(sged_file_groups) as csv_file_groups:
  df_groups = pandas.read_csv(csv_file_groups, sep = delim, dtype = str, keep_default_na = False)
  groups = df_groups[group_col]
  groups_lst = [ g[1:(len(g)-1)].split(";") for g in groups ]

n_groups = len(df_groups.index)

# Now replicate each group:
x_rep = [ 0 ] * (n_groups * n_rep)
l_grp = []
x_grp = [ "" ] * (n_groups * n_rep)
x_ave = [ 0 ] * (n_groups * n_rep) #Average of the sampled group
x_siz = [ 0 ] * (n_groups * n_rep)
x_oav = [ 0 ] * (n_groups * n_rep) #Average the original group

i = 0
for grp in range(n_groups):
  size = df_groups[grp, "Size"] # TODO: if the size column is not present, we should generate it
  nmin = df_groups[grp, cond_var]
  
  # Get all sites with adequate value for each position:
  gp = groups_lst[grp]
  if len(gp) != size:
    raise IOError("!!! Error in input file, group size does not match number or sites!")
  
  gp_vals = [ df_sites[gp[j], cond_var ] for j in range(size) ]

  x_rep[i:(i + nrep)] = [ x for x in range(1, n_rep + 1) ]
  x_siz[i:(i + nrep)] = [ size ] * n_rep
  x_grp[i:(i + nrep)] = [ "[" ] * n_rep
  x_ave[i:(i + nrep)] = [ 0 ] * n_rep
  x_oav[i:(i + nrep)] = [ numpy.mean(gp_vals) ] * n_rep
  
  # Loop over each site in the group:
  for (sit in 1:size) {
    # Get all sites with similar rate (or any other conditional variable):
    x <- sites[,cond.var]
    t <- abs(x - gp.vals[sit]) / gp.vals[sit]
    condSites <- subset(sites, t <= sim.t)
    # Loop over each simulation replicate:
    for (sim in 1:nrep) {    
      # Remove sites already present in the group
      if (i + sim - 1 <= length(l.grp)) { # test if element already exists in list
        tmp <- subset(condSites, ! (Group %in% l.grp[[i + sim - 1]]))
      } else {
        tmp <- condSites
      }

      # Sampling correction:
      tmpl <- tmp[tmp[,cond.var] < gp.vals[sit],]
      tmpg <- tmp[tmp[,cond.var] > gp.vals[sit],]
      tmpe <- tmp[tmp[,cond.var] == gp.vals[sit],]
      n <- min(nrow(tmpl), nrow(tmpg))
      n <- max(n, min.obs) # we add min.obs there to avoid getting no replicate site when we have extreme values for the candidate site.
      tmp2 <- rbind(tmpl[sample(1:nrow(tmpl), min(n, nrow(tmpl))),], tmpe, tmpg[sample(1:nrow(tmpg), min(n, nrow(tmpg))),])
      if (nrow(tmp2) == 0) {
        warning(paste("No more similar site available for candidate site", sit, "in group", grp, "replicate", sim))
        x.grp[i + sim - 1] <- paste(x.grp[i + sim - 1], "NA", sep = ifelse(x.grp[i + sim - 1] == "[", "", ";"))
        x.ave[i + sim - 1] <- NA
      } else {
        if (nrow(tmp2) < min.obs)
          warning(paste("Minimum site frequency not matched for candidate site", sit, "in group", grp, "replicate", sim))

        # Now sample sites:
        x <- sample(1:nrow(tmp2), 1)
        if (i + sim - 1 <= length(l.grp)) {
          l.grp[[i + sim - 1]] <- append(l.grp[[i + sim - 1]], tmp2[x, "Group"])
        } else {
          l.grp[[i + sim - 1]] <- tmp2[x, "Group"]
        }
        x.grp[i + sim - 1] <- paste(x.grp[i + sim - 1], tmp2[x, "Group"], sep = ifelse(x.grp[i + sim - 1] == "[", "", ";"))
        x.ave[i + sim - 1] <- x.ave[i + sim - 1] + tmp2[x, cond.var]
      }
    }
  }
  x.grp[i:(i + nrep - 1)] <- paste(x.grp[i:(i + nrep - 1)], "]", sep = "")
  i <- i + nrep
}

x.ave <- x.ave / x.siz
results <- data.frame(Replicate = x.rep, Group = x.grp, Size = x.siz, RandMean = x.ave, OrigMean = x.oav)

write.table(results, file = outputPath, sep = "\t", row.names = FALSE)



print "Done."

