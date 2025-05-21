#!/bin/bash
# sort and index bam files using sambamba 
# 10 march 2025
#katherine carbeck 

## this was done on medium 64 core machine --

#https://www.basepairtech.com/blog/sorting-bam-files-samtools-vs-sambamba/

mkdir -p /workdir/kcarbeck/tmp

## coordinate based sorting for mpileup
[ -f /workdir/kcarbeck/align/SambambaSortCommands.txt ] && rm /workdir/kcarbeck/align/SambambaSortCommands.txt

cd /workdir/kcarbeck/align/


for BAM in *.bam; do
  BASENAME=${BAM%.bam}
  OUT=${BASENAME}_sorted.bam

  # create sort file
  echo "/programs/sambamba-0.7.1/sambamba sort -t 12 -m 100G --tmpdir /workdir/kcarbeck/tmp -o $OUT $BAM" >> /workdir/kcarbeck/align/SambambaSortCommands.txt
done


parallel -j 5 < /workdir/kcarbeck/align/SambambaSortCommands.txt ; python /workdir/kcarbeck/job_complete_email.py


# total bytes in sorted dir: 2,949,050,693




# check if bam files are sorted 
# Will show SO:coordinate if sorted
samtools view -H *sorted.bam | grep @HD



####
/programs/sambamba-0.7.1/sambamba sort -t 8 -m 150G -o LC_2_117_sorted.bam LC_2_117.bam 

samtools view LC_2_117_sorted.bam | awk '{print $5}' | sort | uniq -c

# first run?
6170753 0
 482738 1
6330229 11
 474924 12
3224503 14
 138535 16
2037127 17
1435340 18
  27700 19
7413337 2
1012925 21
2945834 22
1237197 24
 545440 25
1188303 28
  14266 31
  18900 32
  33599 33
  62215 34
  33322 35
 986915 36
 494082 37
 486830 38
 549375 39
 524923 40
 840292 41
1764733 42
146789765 44
1726572 9

# second run?
6312177 0
 505134 1
6604256 11
 507305 12
3342800 14
 152971 16
2100314 17
1478628 18
  33001 19
8076937 2
1050171 21
3001667 22
1250102 24
 577639 25
1210193 28
  14713 31
  19375 32
  34233 33
  63086 34
  34540 35
1007980 36
 520069 37
 511658 38
 582209 39
 581986 40
 860636 41
1788885 42
148410758 44
1843540 9







samtools view -c -F 4 LC_2_117_sorted.bam
183157330 # first run
186530636 # second run



