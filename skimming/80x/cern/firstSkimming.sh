#!/bin/sh

# MC samples used user/ceballos

export theRND=file_skim_$RANDOM;

export INPUTDIR=/eos/cms/store/user/ceballos/Nero/v2.1;
export SKIMDIR=/eos/cms/store/caf/user/ceballos/Nero/skim_80x;
export TYPE="dm"

if [ $# == 2 ] && [ $2 == 1 ]; then
  echo "AREA 1";
  export INPUTDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/v2.1;
  export SKIMDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x;
  export TYPE="dm"

elif [ $# == 2 ] && [ $2 == 2 ]; then
  echo "AREA 2";
  export INPUTDIR=/eos/cms/store/group/phys_higgs/ceballos/setup80x_ichep/Data/Nero/v2.0/$1;
  export  SKIMDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/$1;
  export TYPE="data"

elif [ $# == 2 ] && [ $2 == 3 ]; then
  echo "AREA 3";
  export INPUTDIR=/mnt/hadoop/cms/store/user/dhsu/Nero/v1.3/$1;
  export SKIMDIR=/data/t3home000/ceballos/Nero/skim_80x/$1;
  export TYPE="data"

else
  echo "AREA 0";
fi

find $INPUTDIR/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > $theRND.txt

for datasetName in `cat $theRND.txt`; do
  echo "Merging dataset: "$INPUTDIR/$datasetName
  mkdir -p $SKIMDIR/$datasetName;
  ls $INPUTDIR/$datasetName|grep root > $SKIMDIR/$datasetName.txt;
  for fileName in `cat $SKIMDIR/$datasetName.txt`; do
    root -l -q -b VVScattering/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/$datasetName/$fileName\",\"$SKIMDIR/$datasetName/$fileName\",\"${TYPE}\",-1,0\);
  done
  rm -f $SKIMDIR/$datasetName.txt;
done

rm -f $theRND.txt;
