#!/bin/sh

export theRND=file_skim_$RANDOM;

export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/setup80x_ichep/Data/Nero/v2.0/$1;
export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/$1;

find $INPUTDIR/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > $theRND.txt

for datasetName in `cat $theRND.txt`; do
  echo "Merging dataset: "$INPUTDIR/$datasetName
  mkdir -p $SKIMDIR/$datasetName;
  ls $INPUTDIR/$datasetName|grep root > $SKIMDIR/$datasetName.txt;
  for fileName in `cat $SKIMDIR/$datasetName.txt`; do
    root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/$datasetName/$fileName\",\"$SKIMDIR/$datasetName/$fileName\",\"data\",0,0,0\);
  done
  rm -f $SKIMDIR/$datasetName.txt;
done

rm -f $theRND.txt;
