# script for prepping mapped reads for variant caller
# katherine carbeck
# 18 may 2025

# A new Bam file will be produced for each step. In the end, you only need to keep samplename_sorted_mark.bam

# Reserved a medium machine (48 cores)
# Bam files from mapping script are indexed and sorted already
# copy into workdir:
cd /lustre2/home/lc736_0001/crossbill/sortedBam/
find . -name "*.bam*" | parallel -j 40 cp -a {} /workdir/kcarbeck/{/}

# done copying?
find /workdir/kcarbeck -maxdepth 1 -type f | wc -l 


#### Prepare the genome: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R=red_crossbill_24Feb2018_V5eIH.fasta O=red_crossbill_24Feb2018_V5eIH.dict
samtools faidx red_crossbill_24Feb2018_V5eIH.fasta



###### mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

# first, check header line and RG info
samtools view -H LC_2_170979_sorted.bam | grep ^@RG
#@RG     ID:LC_2_170979  SM:LC_2_170979

samtools view LC_2_170979_sorted.bam | head -n 100 | grep 'RG:Z:'


# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000
# uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

cd /workdir/kcarbeck

#if file exists, remove it
rm -f markDuplicatesCommands.txt

for BAM in *_sorted.bam; do
  OUT="marked/${BAM/_sorted.bam/_sorted_mark.bam}"
  METRICS="marked/${BAM/_sorted.bam/.metrics.txt}"

  [[ -f "$OUT" ]] && echo "skipping $BAM (already processed)" && continue

  echo "java -Xmx15g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=$BAM OUTPUT=$OUT METRICS_FILE=$METRICS MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" >> markDuplicatesCommands.txt
done

parallel -j 15 < /workdir/kcarbeck/markDuplicatesCommands.txt


find *.bam -maxdepth 1 -type f | wc -l 
find *.txt -maxdepth 1 -type f | wc -l 

for bam in *_sorted_mark.bam; do 
    base=${bam%%_sorted_mark.bam}
    [[ ! -f "${base}.metrics.txt" ]] && echo "Missing: ${base}.metrics.txt"
done
# Missing: LC_GUI_713107.metrics.txt


# java -Xmx50g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=LC_GUI_713107_sorted.bam OUTPUT=marked/LC_GUI_713107_sorted_mark.bam METRICS_FILE=marked/LC_GUI_713107.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
# LC_GUI_713107_sorted_mark.bam
