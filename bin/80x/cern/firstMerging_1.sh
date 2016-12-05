#!/bin/sh

export theRND=file_merging1_$RANDOM;

export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/skim_80x;
export MERGINGDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/merging_80x;

mkdir -p $MERGINGDIR;

find $SKIMDIR/* -name '0000'| awk '{split($1,a,ENVIRON["SKIMDIR"]);split(a[2],b,"/");printf("hadd -f %s/%s.root %s/*.root\n",ENVIRON["MERGINGDIR"],b[3],$1);}' > $theRND.sh;
chmod a+x $theRND.sh;./$theRND.sh;rm -f $theRND.sh;

find $SKIMDIR/* -name '0001'| awk '{split($1,a,ENVIRON["SKIMDIR"]);split(a[2],b,"/");printf("hadd -f %s/%s_0001.root %s/*.root\n",ENVIRON["MERGINGDIR"],b[3],$1);}' > $theRND.sh;
chmod a+x $theRND.sh;./$theRND.sh;rm -f $theRND.sh;

find $SKIMDIR/* -name '0002'| awk '{split($1,a,ENVIRON["SKIMDIR"]);split(a[2],b,"/");printf("hadd -f %s/%s_0002.root %s/*.root\n",ENVIRON["MERGINGDIR"],b[3],$1);}' > $theRND.sh;
chmod a+x $theRND.sh;./$theRND.sh;rm -f $theRND.sh;
