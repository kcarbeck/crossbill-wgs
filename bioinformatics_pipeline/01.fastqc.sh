#!bin/bash
# fastqc/multiqc for raw reads 
# author: katherine carbeck
# 10 feb 2025

#### all the raw files are stored in the shared crossbill dir on Cornell's BioHPC
# /lustre2/home/lc736_0001/crossbill


# 1. copy all files in directory and subdirs without the dir structure 
find /lustre2/home/lc736_0001/crossbill -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -exec cp -t /workdir/kcarbeck/fastqc {} +

# count files in current directory
ls -1 | wc -l

#find /workdir/kcarbeck/fastqc/*.gz -type f > listOfFiles.list
find *.gz -type f > listOfFiles.list


awk -v FS="\t" '
{
    if (NR > 0){
        samp=$1
        print "fastqc"  " " samp
    }
}
' listOfFiles.list > list.txt

# run fastqc
parallel -j 22 < /workdir/kcarbeck/fastqc/list.txt
#started ~10:30pm feb 11, 2025 on cbsumm23 [probably couldve done more at once but was worried about memmory because the files are so large]
# finished ~ 2:00 PM feb 12 

# Run above code on terminal and view output using: 
export PYTHONPATH=/programs/multiqc-1.15/lib64/python3.9/site-packages:/programs/multiqc-1.15/lib/python3.9/site-packages
export PATH=/programs/multiqc-1.15/bin:$PATH
source activate multiqc

# run software using command: 
multiqc --data-format tsv .

# after done running deactivate conda 
conda deactivate
    
# multiqc_data/ has output: multiqc_general_stats.txt

cp -R multiqc_data /lustre2/home/lc736_0001/crossbill/
cp multiqc_report.html /lustre2/home/lc736_0001/crossbill/multiqc_data/

