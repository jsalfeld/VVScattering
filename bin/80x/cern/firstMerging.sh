#!/bin/sh

export theRND=file_merging_$RANDOM;

export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/user/ceballos/Nero/skim_80x;
export MERGINGDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/user/ceballos/Nero/merging_80x;

mkdir -p $MERGINGDIR;

find $SKIMDIR/* -name '000?'| awk '{split($1,a,ENVIRON["SKIMDIR"]);split(a[2],b,"/");printf("hadd -f %s/%s.root %s/*.root\n",ENVIRON["MERGINGDIR"],b[2],$1);}' > $theRND.sh;
chmod a+x $theRND.sh;./$theRND.sh;rm -f $theRND.sh;
