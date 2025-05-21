#script for trimming adapters
#author: katherine carbeck
#3 march 2025

cd /lustre2/home/lc736_0001/crossbill/rawdata
#checksums at source
find . -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0 \
  | xargs -0 md5sum > crossbill_fastqs.md5 &

#make a list of FASTQ files (ignore subdirs)
find . -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) > filelist.txt

#copy them using rsync to one folder in workdir 
rsync -avW --progress --partial \
  --files-from=filelist_filtered.txt \
  --no-relative \
  /lustre2/home/lc736_0001/crossbill/rawdata/ \
  /workdir/kcarbeck/rawdata/ &

#checksums
cd /workdir/kcarbeck/rawdata
md5sum -c /lustre2/home/lc736_0001/crossbill/crossbill_fastqs.md5

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


# run adapterRemovalCommands.txt in parallel on 88 core machine
# this funnels the output to a log file and displays it on the screen
parallel -j 11 < /workdir/kcarbeck/rawdata/adapterRemovalCommands.txt 2>&1 | tee adapterRemoval.log
# started around 10 AM March 24
# finished remaining samples next day



#### transfer files out:
# dry run
rsync -avh --dry-run --progress --delete /workdir/kcarbeck/adapterRemoval/ /lustre2/home/lc736_0001/crossbill/adapterRemoval/

nohup rsync -avh --progress --partial --inplace --delete /workdir/kcarbeck/adapterRemoval/ /lustre2/home/lc736_0001/crossbill/adapterRemoval/ > rsync_log.txt 2>&1 &
# -avh: archive mode, verbose, human readable file sizes
# --partial: keeps partially transferred files and can pick back up if interrupted
# --inplace: writes directly to destination (not sure if this is necessary?)
# --delete: deletes extra files at the destination

# tail -f rsync_log.txt
# Ctrl + C
