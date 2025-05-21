#### Validate files
# since running HaplotypeCaller, don't need to realign or fix mates unless there is an error in ValidateSamFile

[ -f /workdir/kcarbeck/marked/validateCommands.txt ] && rm /workdir/kcarbeck/marked/validateCommands.txt

cd /workdir/kcarbeck/marked/

for MARKBAM in *mark.bam; do

  OUTFILE=${MARKBAM/%.bam/_validate}

  echo "java -Xmx30g -jar /programs/picard-tools-2.8.2/picard.jar ValidateSamFile I=$MARKBAM OUTPUT=$OUTFILE MODE=SUMMARY" >> /workdir/kcarbeck/marked/ValidateCommands.txt

done

parallel -j 9 < /workdir/kcarbeck/marked/ValidateCommands.txt

# concatenate all validate output into one text file to see if there are any errors
cat *validate > summary_validate.txt
