#### index .bam files for 

cd /workdir/kcarbeck/marked/
rm -f IndexCommands.txt

for MARKBAM in *mark.bam; do

  #create index file
  [[ -f "${MARKBAM}.bai" ]] && echo "skipping $MARKBAM (index exists)" && continue
  echo "java -Xmx2g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=$MARKBAM" >> IndexCommands.txt

done

parallel -j 24 < IndexCommands.txt