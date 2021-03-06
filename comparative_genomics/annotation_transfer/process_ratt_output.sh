#!/bin/bash
#Usage: bash process_ratt_output.sh outbasename
#Dependencies: genometools, find, awk

NAME=$1

#convert all to gff3
echo "converting embl to gff"
mkdir gff3
find embl -type f -name '*.final.embl' -exec bash -c 'seqret -sequence ${1} -feature -outseq gff3::${1//embl/gff3}' - {} \;

echo "reformat and concatenate gffs"
find gff3 -maxdepth 1 -type f -name "*.final.gff3" -exec cat '{}' ';' |dd obs=1M |
awk -v NAME="^${NAME}." '$0 ~ NAME {gsub(NAME,"",$0) gsub(".final","",$0); print $0}' |
awk  -v OFS='\t' '!/locus_tag/ && $3=="CDS" {print} /locus_tag/ && $3=="CDS" {$3="mRNA";print;$3="CDS";gsub(";.+","",$9);gsub("ID=","Parent=",$9);print} /locus_tag/ && $3=="biological_region" {$3="mRNA";print}' |
awk '$4!=$5' |
gt gff3 -tidy -retainids -sort \
>${NAME}.final.gff3

echo "create filter"
cd embl
rm ${NAME}.filter
find . -maxdepth 1 -type f  -name '*final.embl' -exec awk 'BEGIN{RS="FT   CDS"}($3~/\.1\.2/){gsub(/\.1\..+\"/, ".1\"",$3);print $3}' {} >>${NAME}.filter \;

echo "filter output" 
AWK='(FNR==NR){a[$0]; next}{RS="\nFT   CDS"} ($1~"ID"||$1~"FH"){print $0;next} {gsub(/\.1\..+\"/, ".1\"")} !($3 in a) {print "FT   CDS" $0} END{if ($3 in a) {n=split($0,E,"SQ   "); print "SQ   " E[n]}}'
find . -maxdepth 1 -type f -name '*final.embl' -exec bash -c 'awk "$1" "$2" "$3" >"${3/final./final.filtered.}"' - "$AWK" "$NAME.filter" {} \;
cd ..

echo "converting filtered embl to gff"
find embl -maxdepth 1 -type f -name '*.final.filtered.embl' -exec bash -c 'seqret -sequence ${1} -feature -outseq gff3::${1//embl/gff3}' - {} \;

echo "reformat and concatenate filtered gffs"
find gff3 -maxdepth 1 -type f -name "*.final.filtered.gff3" -exec cat '{}' ';' |dd obs=1M|
awk -v NAME="^${NAME}." '$0 ~ NAME {gsub(NAME,"",$0) gsub(".final","",$0); print $0}' |
awk  -v OFS='\t' '!/locus_tag/ && $3=="CDS" {print} /locus_tag/ && $3=="CDS" {$3="mRNA";print;$3="CDS";gsub(";.+","",$9);gsub("ID=","Parent=",$9);print} /locus_tag/ && $3=="biological_region" {$3="mRNA";print}' |
awk '$4!=$5' |
gt gff3 -tidy -retainids -sort \
>${NAME}.final.filtered.gff3

exit

echo "clean up files"
find Sequences -maxdepth 1 -type f -name '*' | xargs -s 20000 zip sequences.zip -@ && rm -r Sequences

find reports -maxdepth 1 -type f -name '*' | xargs -s 20000 zip reports.zip -@ && rm -r reports

cd embl
find . -maxdepth 1 -type f -name '*final.embl' | xargs -s 20000 zip final.embl.zip -@ && \
find . -maxdepth 1 -type f -name '*final.filtered.embl' | xargs -s 20000 zip final.filtered.embl.zip -@ && \
find . -maxdepth 1 -type f -name '*NOTTransfered.embl' | xargs -s 20000 zip nottransfered.embl.zip -@ && \
find . -maxdepth 1 -type f -name '*final*embl' | xargs -s 20000 rm && \
find . -maxdepth 1 -type f -name '*NOTTransfered.embl*' | xargs -s 20000 rm && \
find . -maxdepth 1 -type f -name '*.embl' | xargs -s 20000 zip embl.zip -@ && \
find . -maxdepth 1 -type f -name '*.embl' | xargs -s 20000 rm
cd ../

cd gff3
find . -maxdepth 1 -type f -name '*final.gff3' | xargs -s 20000 zip final.gff3.zip -@ && \
find . -maxdepth 1 -type f -name '*final.filtered.gff3' | xargs -s 20000 zip final.filtered.gff3.zip -@ && \
find . -maxdepth 1 -type f -name '*final*gff3' | xargs -s 20000 rm
cd ../

rm -r Query/
rm -r Reference/

exit
