#!/bin/sh

export theRND=file_skim_$RANDOM;

export INPUTDIR=/mnt/hadoop/cms/store/user/dhsu/Nero/v1.3/$1;
export SKIMDIR=/data/t3home000/ceballos/Nero/skim_80x;

find $INPUTDIR/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > $theRND.txt

for datasetName in `cat $theRND.txt`; do
  echo "Merging dataset: "$INPUTDIR/$datasetName
  ls $INPUTDIR/$datasetName|grep root > $INPUTDIR/$datasetName.txt;
  mkdir -p $SKIMDIR/$datasetName;
  for fileName in `cat $INPUTDIR/$datasetName.txt`; do
    root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/$datasetName/$fileName\",\"$SKIMDIR/$datasetName/$fileName\",\"data\",0,0,0\);
  done
  rm -f $INPUTDIR/$datasetName.txt;
done

rm -f $theRND.txt;
