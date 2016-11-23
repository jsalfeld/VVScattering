#!/bin/sh

export theRND=file_skim_$RANDOM;

export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/user/ceballos/Nero/v1.4;
export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/user/ceballos/Nero/skim_80x;

find $INPUTDIR/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > $theRND.txt

for datasetName in `cat $theRND.txt`; do
  echo "Merging dataset: "$INPUTDIR/$datasetName
  ls $INPUTDIR/$datasetName|grep root > $theRND_$INPUTDIR/$datasetName.txt;
  mkdir -p $SKIMDIR/$datasetName;
  for fileName in `cat $theRND_$INPUTDIR/$datasetName.txt`; do
    root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/$datasetName/$fileName\",\"$SKIMDIR/$datasetName/$fileName\",\"data\",0,0,0\);
  done
  rm -f $theRND_$INPUTDIR/$datasetName.txt;
done
