import pandas as pd
import getopt, sys

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "p:o:m:"
full_opt = ["paml=", "output=", "method="]


try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

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


if method == "bayesian":
    #read the paml result file
    with open(paml_file, 'r') as f:
        contents = f.readlines()
    #find the positive selected sites under the BEB
    start_row = [i for i, line in enumerate(contents) if 'Bayes Empirical Bayes (BEB)' in line][0] + 3
    end_row = [i for i, line in enumerate(contents) if 'The grid' in line][0] - 2

elif method == "naive":
    # read the paml result file
    with open(paml_file, 'r') as f:
        contents = f.readlines()
    start_row = [i for i, line in enumerate(contents) if 'Naive Empirical Bayes (NEB)' in line][0] + 3
    end_row = [i for i, line in enumerate(contents) if 'Bayes Empirical Bayes (BEB)' in line][0] - 2

else:
    print("Error: Method not specified.")
    sys.exit(2)

lines = contents[start_row:end_row]

# Extract the information
positive_sites = []
for line in lines:
    site_info = line.split()
    positive_sites.append(site_info)

# convert it to the data frame and add square brackets
df = pd.DataFrame(positive_sites, columns=['position', 'amino_acid', 'probability'])
df.insert(loc=0, column='Groups', value='[' + df['position'] + ']')
df.drop(['position'], axis=1, inplace=True)
print(df)

# converting to the csv file
df.to_csv(output_file, index=False)
