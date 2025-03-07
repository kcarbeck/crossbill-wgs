#script for trimming adapters
#author: katherine carbeck
#3 march 2025

##### ADAPTER REMOVAL

#Trimns: Trim consecutive Ns from the 5’ and 3’ termini. If quality trimming is also enabled (--trimqualities), then stretches of mixed low-quality bases and/or Ns are trimmed.
#Trimqualities: Trim consecutive stretches of low quality bases (threshold set by --minquality) from the 5’ and 3’ termini. If trimming of Ns is also enabled (--trimns), then stretches of mixed low-quality bases and Ns are trimmed.
#minquality: Set the threshold for trimming low quality bases using --trimqualities and --trimwindows. Default is 2
#collapse: In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. In single-end mode, attempt to identify templates for which the entire sequence is available. In both cases, complete “collapsed” reads are written with a ‘M_’ name prefix, and “collapsed” reads which are trimmed due to quality settings are written with a ‘MT_’ name prefix.
#adapterlist: Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to –adapter1 and –adapter2.
#minlength: Reads shorter than this length are discarded following trimming. Defaults to 15Reads shorter than this length are discarded following trimming. Defaults to 15.

# Usage: /programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 Pathto/R1.fastq.gz --file2 Pathto/R2.fastq.gz --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 8 --adapter-list for_adapter_removal_2.txt --basename LIBRARY_SAMPLEID

# output of python helper script should look like this:
/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 DP8400010194TL_L01_SP2004210287_1.fq.gz --file2 DP8400010194TL_L01_SP2004210287_2.fq.gz --adapter-list for_adapter_removal.txt --basename LC_2_170979 --trimns --trimqualities --minquality 35 --minlength 25 --collapse --threads 8 --gzip


# run adapterRemovalCommands.txt in parallel on 40 core machine
# this funnels the output to a log file and displays it on the screen
parallel -j 5 < /workdir/kcarbeck/rawdata/adapterRemovalCommands.txt 2>&1 | tee adapterRemoval.log
# started around 2PM 3 Mar 25; probably could've run more files at once?


##################     4 March 2025    ##################
# failed server ran out of storage

# sync completed files
rsync -av --progress /workdir/kcarbeck/rawdata/ /lustre2/home/lc736_0001/crossbill/rawdata/

# copy all files in directory and subdirs without the dir structure 
find /lustre2/home/lc736_0001/crossbill -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -exec cp -t /workdir/kcarbeck/rawdata/ {} +

ls -1 | wc -l

#remove files that have already been completed
xargs rm < 02.X.filesToDelete.txt


## rerun adapter removal helper script to generate new text file with the files that are left to be completed
# updated this file to reflect the files that were successfully completed: crossbill_metadataForAdapterRemoval.csv



##run adapterRemovalCommands.txt in parallel on 64 core machine
#this funnels the output to a log file and displays it on the screen
parallel -j 10 < /workdir/kcarbeck/rawdata/adapterRemovalCommands.txt 2>&1 | tee adapterRemoval.log
# started around 4PM 4 Mar 25



##########   March 5 2025    #############
# finished by the morning but had a few errors regarding gzipped files and samples that had more than one fastq associated with them
# gzip_paired_fastq::~gzip_paired_fastq: data error


## RUNNING CHECKS...
# to see if fastq files are not corrupted
# and to see if the adapter removal output is as expected

## update: ultimately found that one preferred sample name was duplicated so updated those IDs and ran again; 


