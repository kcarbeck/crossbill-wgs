#!/bin/bash
# get mapping stats using qualimap
# 3 april 2025
#katherine carbeck 

## this was done on large 88 core machine --

## first, copy out unsorted bam files and make sure they were copied corectly
## this took forever...def make this script more efficient for future use
SRC_DIR="/workdir/kcarbeck/align"
DEST_DIR="/lustre2/home/lc736_0001/crossbill/align"
CHECKSUM_FILE="original_bam_checksums.md5"


echo "generating checksums in source dir..."
cd "$SRC_DIR" || { echo "source directory not found"; exit 1; }

# only create if not already done
if [ ! -f "$CHECKSUM_FILE" ]; then
  md5sum *.bam > "$CHECKSUM_FILE"
  echo "saved to $CHECKSUM_FILE"
else
  echo "using existing $CHECKSUM_FILE"
fi

echo "copying file to destination"
cp "$CHECKSUM_FILE" "$DEST_DIR/"

echo "verifying copied files at destination"
cd "$DEST_DIR" || { echo "destination directory not found"; exit 1; }

# Run checksum verification
md5sum -c "$CHECKSUM_FILE"

echo "verification complete :p"



############    QUALIMAP   ############
# may have to increase mem limit:
#/programs/qualimap_v2.2.1/qualimap bamqc -bam SAMPLE.bam --java-mem-size=30G -outfile SAMPLE.sorted.pdf

# set up custom tmpdir to avoid filling... was getting weird silent failing and this fixed it
export TMPDIR=/workdir/kcarbeck/tmp
mkdir -p $TMPDIR



[ -f /workdir/kcarbeck/sorted/qualimapCommands.txt ] && rm /workdir/kcarbeck/sorted/qualimapCommands.txt

cd /workdir/kcarbeck/sorted/

for BAM in *.bam; do

  #string replacement command
  BASENAME=${BAM%_sorted.bam}
  OUTDIR=/workdir/kcarbeck/bamqc/${BASENAME}_bamqc

  # qualimap
  echo "/programs/qualimap_v2.2.1/qualimap bamqc -bam $BAM -outdir $OUTDIR -outformat PDF:HTML --java-mem-size=40G" >> /workdir/kcarbeck/sorted/qualimapCommands.txt

done

parallel -j 11 < /workdir/kcarbeck/sorted/qualimapCommands.txt ; python /workdir/kcarbeck/job_complete_email.py



# to run and get log output:
# bash qualimap.sh 2>&1 | tee qualimap_may1_$(date +%Y%m%d-%Hh%Mm%Ss).log





