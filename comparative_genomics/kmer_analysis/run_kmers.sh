#!/bin/bash

usage="$(basename "$0") [-h] input-reads -- script to run kmer analysis using jellyfish on .fastq sequence reads

where:
    -h  show this help text"

while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

FILES=($*)

for F in "${FILES[@]}"
do
echo $F
#test if input is a file
if ! [ -f "$1" ]
then
    echo "$1 is not a file"
    exit
fi

done


/mnt/nexenta/TremaParasponia/software/jellyfish-2.2.0/bin/jellyfish count -m 21 -s 100M -t 16 -C $* &&\
/mnt/nexenta/TremaParasponia/software/jellyfish-2.2.0/bin/jellyfish histo mer_counts.jf > mer_histo.txt &&\

R --no-save </mnt/nexenta/TremaParasponia/scripts/parse_jellyfish_output.R

exit
