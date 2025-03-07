# script to check all the .settings files output by adapterRemoval to see if everything processed as expected
# katherine carbeck
# 5 march 2025


import os
import pandas as pd


def parse_full_settings_file(filepath):
    """Extracts all relevant statistics from an AdapterRemoval .settings file with corrected derived calculations."""
    metrics = {
        "Filename": os.path.basename(filepath),
        "Total Read Pairs": None,
        "Unaligned Read Pairs": None,
        "Well Aligned Read Pairs": None,
        "Inadequate Alignments": None,
        "Retained Reads": None,
        "Discarded Reads": None,
        "Singleton Mate 1 Reads": None,
        "Singleton Mate 2 Reads": None,
        "Adapters Found": None,
        "Full-Length Collapsed Pairs": None,
        "Truncated Collapsed Pairs": None,
        "Retained Nucleotides": None,
        "Average Read Length": None,
    }
    discarded_1, discarded_2 = 0, 0  # Temp variables for summing discarded reads
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if "Total number of read pairs" in line:
                metrics["Total Read Pairs"] = int(line.split(":")[-1].strip())
            elif "Number of unaligned read pairs" in line:
                metrics["Unaligned Read Pairs"] = int(line.split(":")[-1].strip())
            elif "Number of well aligned read pairs" in line:
                metrics["Well Aligned Read Pairs"] = int(line.split(":")[-1].strip())
            elif "Number of inadequate alignments" in line:
                metrics["Inadequate Alignments"] = int(line.split(":")[-1].strip())
            elif "Number of retained reads" in line:
                metrics["Retained Reads"] = int(line.split(":")[-1].strip())
            elif "Number of discarded mate 1 reads" in line:
                discarded_1 = int(line.split(":")[-1].strip())
            elif "Number of discarded mate 2 reads" in line:
                discarded_2 = int(line.split(":")[-1].strip())
                metrics["Discarded Reads"] = discarded_1 + discarded_2  # sum both mate reads
            elif "Number of singleton mate 1 reads" in line:
                metrics["Singleton Mate 1 Reads"] = int(line.split(":")[-1].strip())
            elif "Number of singleton mate 2 reads" in line:
                metrics["Singleton Mate 2 Reads"] = int(line.split(":")[-1].strip())
            elif "Number of reads with adapters[0]" in line:
                metrics["Adapters Found"] = int(line.split(":")[-1].strip())
            elif "Number of full-length collapsed pairs" in line:
                metrics["Full-Length Collapsed Pairs"] = int(line.split(":")[-1].strip())
            elif "Number of truncated collapsed pairs" in line:
                metrics["Truncated Collapsed Pairs"] = int(line.split(":")[-1].strip())
            elif "Number of retained nucleotides" in line:
                metrics["Retained Nucleotides"] = int(line.split(":")[-1].strip())
            elif "Average read length of trimmed reads" in line:
                metrics["Average Read Length"] = float(line.split(":")[-1].strip())
    # derived values
    total_read_pairs = metrics["Total Read Pairs"]
    discarded_reads = metrics["Discarded Reads"] if metrics["Discarded Reads"] is not None else 0
    full_length_collapsed = metrics["Full-Length Collapsed Pairs"] if metrics["Full-Length Collapsed Pairs"] is not None else 0
    adapters_found = metrics["Adapters Found"] if metrics["Adapters Found"] is not None else 0
    if total_read_pairs:
        metrics["% Reads Retained"] = round(((total_read_pairs - discarded_reads) / total_read_pairs) * 100, 2)
        metrics["% Reads Discarded"] = round((discarded_reads / total_read_pairs) * 100, 2)
        metrics["% Collapsed Pairs"] = round((full_length_collapsed / total_read_pairs) * 100, 2)
        metrics["% Adapter Occurrence"] = round((adapters_found / total_read_pairs) * 100, 2)
    return metrics


# set dir
directory_path = "/workdir/kcarbeck/adapterRemoval" 

# make sure it's finding the dir
if not os.path.isdir(directory_path):
    print(f"Error: directory '{directory_path}' not found")
else:
    # get all settings files
    settings_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith(".settings")]
    if not settings_files:
        print("no settings files found")
    else:
        print(f"processing {len(settings_files)} settings files...")
        # process all files
        all_settings_data = [parse_full_settings_file(f) for f in settings_files]
        # convert to df
        df_full_summary = pd.DataFrame(all_settings_data)
        # save 
        output_file = os.path.join(directory_path, "adapter_removal_summary.csv")
        df_full_summary.to_csv(output_file, index=False)
        # print summary
        print(f"processed {len(settings_files)} settings files!")
        print(f"summary saved to {output_file}")
        print(df_full_summary.head())


