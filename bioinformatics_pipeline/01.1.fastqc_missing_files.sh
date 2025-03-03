#!bin/bash
# fastqc/multiqc for raw reads -- run for files that were missing from google drive
# author: katherine carbeck
# 28 feb 2025


#### all the raw files are stored in the shared crossbill dir on Cornell's BioHPC
# /lustre2/home/lc736_0001/crossbill


#### 1. copy all files in directory and subdirs without the dir structure 
find /lustre2/home/lc736_0001/crossbill -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -exec cp -t /workdir/kcarbeck {} +

# count files in current directory
ls -1 | wc -l
# 462 -- correct



#### 2. make a file containing all 462 filenames
# -maxdepth limits to no subdirs
# basename to limit to filename without path
find . -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) \
    | xargs -n 1 basename \
    | sort > all_files.txt

head all_files.txt



#### 3. use multiqc_fastqc.txt "Filename" columm to figure out what files have already been completed
cd /workdir/kcarbeck/multiqc_data


awk -F '\t' 'NR>1 {print $2}' multiqc_fastqc.txt | sort > completed_files.txt

head completed_files.txt



#### 4. figure out what files are in all_files.txt and not in completed_files.txt
cd /workdir/kcarbeck

comm -23 all_files.txt completed_files.txt > missing_files.txt



#### 5. now run fastqc and multiqc for missing files

awk -v FS="\t" '
{
    if (NR > 0){
        samp=$1
        print "fastqc"  " " samp
    }
}
' missing_files.txt > missinglist.txt


# run fastqc
parallel -j 14 < /workdir/kcarbeck/missing_fastqc/missinglist.txt


# Run above code on terminal and view output using: 
export PYTHONPATH=/programs/multiqc-1.15/lib64/python3.9/site-packages:/programs/multiqc-1.15/lib/python3.9/site-packages
export PATH=/programs/multiqc-1.15/bin:$PATH
source activate multiqc

multiqc --data-format tsv .
## was having trouble with errors with python version this time for some reason


# after done running deactivate conda 
conda deactivate


### one of the new files was corrupt