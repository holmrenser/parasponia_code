#!/bin/bash
#Usage: bash convert_gff-to-embl.ah [input gff] [input fasta] [output directory name] [species name] 
#dependencies: genometools (http://genometools.org/), gff3_to_embl.lua 

GFF=`readlink -f "${1}"`
FAS=`readlink -f "${2}"`
NAME=${GFF##*\/}

mkdir $3
cd $3

echo "tidying gff3 input file"
gt gff3 -tidy -retainids -sort  ${GFF} > ${NAME/.gff/.clean.gff} 

echo "converting gff3 file to embl format"
mkdir embl
cd embl
gff3_to_embl.lua ../${NAME/.gff/.clean.gff} '${4}' ${FAS}
