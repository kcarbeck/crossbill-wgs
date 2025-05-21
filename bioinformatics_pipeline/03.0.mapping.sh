#!/bin/bash
# align to ref genome new batch of samples (both zipped and unzipped)
#7 march 2025
#katherine carbeck 

## this was done on large 80 core machine --

### INDEX REFERENCE GENOME
bowtie2-build -f /workdir/kcarbeck/reference/red_crossbill_24Feb2018_V5eIH.fasta red_crossbill_index


### ALIGN USING BOWTIE2 > output piped into samtools to avoid sam file output

#load your sample names into an array.
#in this example I will get the sample names from the adapter removal files I generated after a previous step
mkdir /workdir/kcarbeck/align
mkdir /workdir/kcarbeck/logs





### test run
bash -c 'echo Aligning LC_2_117 with bowtie; bowtie2 --mm -p 8 --phred33 --very-sensitive-local -x /workdir/kcarbeck/reference/red_crossbill_index -I 149 -X 900 --rg-id LC_2_117 --rg SM:LC_2_117 -1 /workdir/kcarbeck/adapterRemoval/LC_2_117.pair1.truncated.gz -2 /workdir/kcarbeck/adapterRemoval/LC_2_117.pair2.truncated.gz -U /workdir/kcarbeck/adapterRemoval/LC_2_117.collapsed.gz,/workdir/kcarbeck/adapterRemoval/LC_2_117.collapsed.truncated.gz,/workdir/kcarbeck/adapterRemoval/LC_2_117.singleton.truncated.gz | samtools view -bS > /workdir/kcarbeck/align/LC_2_117.bam' 2> /workdir/kcarbeck/logs/LC_2_117.log



INDS=($(for i in /workdir/kcarbeck/adapterRemoval/*.settings; do echo $(basename -s .settings "$i"); done))
#basename - will remove the directory path and returns the file name. -s tell it which suffix to remove from the end of the file name (in this case .settings)
#note: the variable INDS will now contain an array of the sample names extracted from the files names

echo "${INDS[@]}"
#If you want to see what is stored in this variable you can type:
#(@=all the elements in the array)

### use loop to write all commands to a text file, then use text file to run in parallel:
REFERENCE=/workdir/kcarbeck/reference/red_crossbill_index

# Modified this loop to map files and pipe to a bam file
[ -f /workdir/kcarbeck/align/bowtie2Commands.txt ] && rm /workdir/kcarbeck/align/bowtie2Commands.txt

for SAMPLEID in "${INDS[@]}"; do
  #declare variables. This makes it easier and neater to write your command line and you just have to change these for future projects.
  ONESEQ=/workdir/kcarbeck/adapterRemoval/${SAMPLEID}.pair1.truncated.gz
  TWOSEQ=/workdir/kcarbeck/adapterRemoval/${SAMPLEID}.pair2.truncated.gz
  USEQ=/workdir/kcarbeck/adapterRemoval/${SAMPLEID}.collapsed.gz,/workdir/kcarbeck/adapterRemoval/${SAMPLEID}.collapsed.truncated.gz,/workdir/kcarbeck/adapterRemoval/${SAMPLEID}.singleton.truncated.gz
  OUTPUT=/workdir/kcarbeck/align/${SAMPLEID}.bam
  LOGFILE=/workdir/kcarbeck/logs/${SAMPLEID}.log

  # skip samples already aligned --> helpful if you need to interrupt and restart process
  [ -f "$OUTPUT" ] && continue
  
  #this just writes a line telling you which sample is being worked on, aligns with bowtie and the output is piped directly into samtools to avoid having the intermediate sam file
  echo "bash -c 'echo Aligning $SAMPLEID with bowtie; bowtie2 --mm -p 8 --phred33 --very-sensitive-local -x $REFERENCE -I 149 -X 900 --rg-id $SAMPLEID --rg SM:$SAMPLEID -1 $ONESEQ -2 $TWOSEQ -U $USEQ | samtools view -bS > $OUTPUT' 2> $LOGFILE" >> /workdir/kcarbeck/align/bowtie2Commands.txt
done

### NOW run using parallel & send email when complete
parallel -j 10 < /workdir/kcarbeck/align/bowtie2Commands.txt ; python /workdir/kcarbeck/job_complete_email.py


##!-------  MONITORING AS SCRIPT IS RUNNING  -------##
# monitor I/O
iostat -xm 5

#to see status as job is running
watch -n 5 'ls -lh /workdir/kcarbeck/align/*.bam' #every 5 seconds
Ctrl + C #to escape

# track failed jobs
grep -L "overall alignment rate" /workdir/kcarbeck/logs/*.log

## count files
ls | grep -E '^[^ ]+\.bam$' | grep -v '_sorted\.bam$' | wc -l

##!------------------------------------------------------##


# reminder:
# LC_T3_171252-3 = corrupted file



# LC_2_117
100617932 reads; of these:
  88372742 (87.83%) were paired; of these:
    5553989 (6.28%) aligned concordantly 0 times
    65770554 (74.42%) aligned concordantly exactly 1 time
    17048199 (19.29%) aligned concordantly >1 times
    ----
    5553989 pairs aligned concordantly 0 times; of these:
      322909 (5.81%) aligned discordantly 1 time
    ----
    5231080 pairs aligned 0 times concordantly or discordantly; of these:
      10462160 mates make up the pairs; of these:
        5472280 (52.31%) aligned 0 times
        1629776 (15.58%) aligned exactly 1 time
        3360104 (32.12%) aligned >1 times
  12245190 (12.17%) were unpaired; of these:
    361064 (2.95%) aligned 0 times
    8589003 (70.14%) aligned exactly 1 time
    3295123 (26.91%) aligned >1 times
96.91% overall alignment rate



### OLD NOTES
### had to use new version of parallel to pipe zcat in parallel: 
# installed to my home dir; run this each log in to use it:
# echo 'alias parallel="$HOME/bin/parallel"' >> ~/.bashrc 
# source ~/.bashrc
# parallel --version # should be 20250222
#### UPDATE: this isn't necessary anymore but may be useful in future