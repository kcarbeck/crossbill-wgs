#helper script to set up filenames for adapter removal
#author: katherine carbeck
#3 march 2025

# to load python on cluster
# module load python/3.12.7

import pandas as pd

#df = pd.read_csv("/workdir/kcarbeck/crossbill_metadata_masterfile_030325.csv")

df = pd.read_csv("/workdir/kcarbeck/rawdata/crossbill_metadataForAdapterRemoval.csv")


commands = [] # create a list to store commands

# for each row (i.e., each sample), gather the FASTQ files for R1 and R2
for idx, row in df.iterrows():
    # sample ID using ID column (preferred ID by crossbill group)
    sample_id = row["ID"] 
    # create a pattern for each of the r1 and r2 columns 
    r1_cols = [c for c in df.columns if "fastq_r1" in c] 
    r2_cols = [c for c in df.columns if "fastq_r2" in c]  
    # Filter out NaNs/empty cells in each row
    r1_files = [str(row[c]) for c in r1_cols if pd.notnull(row[c])]
    r2_files = [str(row[c]) for c in r2_cols if pd.notnull(row[c])]
    # now build a list of files for AdapterRemoval (requires only spaces between files)
    r1_string = " ".join(r1_files)
    r2_string = " ".join(r2_files)
    # build the AdapterRemoval command (with multiple FASTQs for some files)
    #  add whatever args needed here
    command_line = (
        f"/programs/adapterremoval_2.1.1/bin/AdapterRemoval "
        f"--file1 {r1_string} "
        f"--file2 {r2_string} "
        f"--adapter-list for_adapter_removal.txt "
        f"--basename {sample_id} "
        f"--trimns --trimqualities --minquality 35 --minlength 25 "
        f"--collapse --threads 8 --gzip"
    )
    commands.append(command_line)

# write the commands to a text file
with open("adapterRemovalCommands.txt", "w") as f:
    for cmd in commands:
        f.write(cmd + "\n")


# wc -l adapterRemovalCommands.txt
# 217 samples 
# take this txt file to run adapter removal in parallel 

## UPDATE: the adapter removal manual says you can add multiple R1s and R2s in the --file1 and --file2 arguments and it will concatenate them for you... this doesn't appear to be true or i'm doing something wrong...
## i ended up concatenating these files and running them as normal. that said this script is a bit more complicated than needed due to that 