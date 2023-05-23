# Libraries
import sys, getopt
import pandas as pd

# Variable input.
cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "r:o:"
full_opt = ["raser=", "output="]


try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

for arg, val in arguments:
    if arg in ("-r", "--raser"):
        raser_file = val
        print("RASER input file: %s" % raser_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)


# Read the raser file
with open(raser_file, "r") as f:
    contents = f.readlines()

# Select the part of the raser file we want to use
contents = contents[13:]

positive_sites = []
for line in contents:
    site_info = line.split("\t")
    positive_sites.append(site_info)

# Convert the data and add square brackets
df = pd.DataFrame(positive_sites, columns=["position", "amino_acid", "probability", "Proba > 0.95"])
df.insert(loc=0, column="Group", value="[" + df["position"] + "]")
df.drop(["position"], axis=1, inplace=True)
df.replace(to_replace=[None], value="", inplace=True)

for counter, param in enumerate(zip(df["probability"], df["Proba > 0.95"])):
    df["probability"][counter] = param[0].replace("\n", "")
    df["Proba > 0.95"][counter] = param[1].replace("\n", "")
print(df)

# Export to csv
df.to_csv(output_file, index=False)
