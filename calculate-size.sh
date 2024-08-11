#!/bin/bash 
echo "Downloading and measuring genomes"
#download genomes where hits were found
parallel -a $1 -C '\t' -j8 --delay 1.5 "efetch -db nuccore -id {1} -format fasta >{1}.fasta"  


#get sizes of each chromosome with faidx
VAR2=$(echo *.fasta | xargs ls)
for y in ${VAR2}
   do
      faidx ${y} -i chromsizes > ${y}.size
   done
#remove whole chromosomes
rm *.fasta
rm *.fai

cat *.size>>chromosome.sizes

