#helper script to set up filenames for adapter removal
#author: katherine carbeck
#3 march 2025

# to load python on cluster
# module load python/3.12.7

import pandas as pd
import os

fastq_dir = "/workdir/kcarbeck/rawdata"

#df = pd.read_csv("/workdir/kcarbeck/crossbill_metadata_masterfile_030325.csv")
df = pd.read_csv("/workdir/kcarbeck/crosbill_adapterRemoval.csv")

commands = [] # create a list to store commands

for _, row in df.iterrows():
    sample_id = row["ID"] 
    r1_file = row["r1"]
    r2_file = row["r2"]
    r1_path = os.path.join(fastq_dir, r1_file)
    r2_path = os.path.join(fastq_dir, r2_file)
    # only make command if both R1 and R2 files exists in dir
    if os.path.exists(r1_path) and os.path.exists(r2_path):
        command = (
            f"/programs/adapterremoval_2.1.1/bin/AdapterRemoval "
            f"--file1 {r1_file} "
            f"--file2 {r2_file} "
            f"--adapter-list for_adapter_removal.txt "
            f"--basename {sample_id} "
            f"--trimns --trimqualities --minquality 35 --minlength 25 "
            f"--collapse --threads 8 --gzip"
        )
        commands.append(command)



# write commands to file
with open(os.path.join(fastq_dir, "adapterRemovalCommands.txt"), "w") as f:
    for cmd in commands:
        f.write(cmd + "\n")



# wc -l adapterRemovalCommands.txt
# 217 samples 
# take this txt file to run adapter removal in parallel 

## UPDATE: the adapter removal manual says you can add multiple R1s and R2s in the --file1 and --file2 arguments and it will concatenate them for you... this doesn't appear to be true or i'm doing something wrong...
## i ended up concatenating these files and running them as normal. that said this script is a bit more complicated than needed due to that 