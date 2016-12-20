#!/bin/sh

export theRND=file_merging_$RANDOM;

export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/skim_80x;
export MERGINGDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/merging_80x;
if [ $# == 1 ] && [ $1 == 1 ]; then
  export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x;
  export MERGINGDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
elif [ $# == 1 ] && [ $1 == 1 ]; then
  export SKIMDIR=/data/t3home000/ceballos/Nero/skim_80x;
  export MERGINGDIR=/data/t3home000/ceballos/Nero/merging_80x;
fi

mkdir -p $MERGINGDIR;

rm -f $theRND.sh;
for folder in 0000 0001 0002 0003 0004 0005;
do
    export theFolder=${folder}
    echo "checking "${theFolder}
    find $SKIMDIR/* -name ${theFolder}| awk '{split($1,a,ENVIRON["SKIMDIR"]);split(a[2],b,"/");printf("hadd -f %s/%s_%s.root %s/*.root\n",ENVIRON["MERGINGDIR"],b[3],ENVIRON["theFolder"],$1);}' >> $theRND.sh;
done
sort -u $theRND.sh > $theRND_sorted.sh; mv $theRND_sorted.sh $theRND.sh;
#chmod a+x $theRND.sh;./$theRND.sh;rm -f $theRND.sh;
