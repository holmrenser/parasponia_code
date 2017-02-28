#!/bin/bash
#Usage: bash run_ratt.sh embl_dir query.fasta outbasename
#Dependencies: RATT (http://ratt.sourceforge.net/)

RATT_HOME=/home/velze002/progs_nobackup/ratt-code/; export RATT_HOME


EMBL=`readlink -f "${1}"`
FAS=`readlink -f "${2}"`
mkdir $3
cd $3


#$RATT_HOME/start.ratt.sh <Directory with embl-files> <Query-fasta sequence> <Resultname> <Transfer type>
$RATT_HOME/start.ratt.sh $EMBL ${FAS} $3 \
Custom

#cleanup
printf '%s\0' *tmp2* | xargs -0 rm

mkdir reports
printf '%s\0' *.Report.* | xargs -0 mv -t reports/

mkdir embl
printf '%s\0' *.embl | xargs -0 mv -t embl/


