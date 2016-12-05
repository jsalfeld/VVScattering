#!/bin/sh
export INPUTDIR=/data/t3home000/ceballos/unmerged_ntuples;
export SKIMTEMPDIR=/data/t3home000/ceballos/ntuples_skim_80x/temp;
export SKIMDIR=/data/t3home000/ceballos/ntuples_skim_80x;
export GOODRUNDIR=/data/t3home000/ceballos/ntuples_goodrun_80x;
export JSONFILE=MitAnalysisRunII/json/80x/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt;

export theRND=file_Run2016H_$RANDOM;

cat > $theRND.txt <<EOF
DoubleEG+Run2016H-PromptReco-v1+AOD
DoubleEG+Run2016H-PromptReco-v2+AOD
DoubleEG+Run2016H-PromptReco-v3+AOD
DoubleMuon+Run2016H-PromptReco-v1+AOD
DoubleMuon+Run2016H-PromptReco-v2+AOD
DoubleMuon+Run2016H-PromptReco-v3+AOD
MuonEG+Run2016H-PromptReco-v2+AOD
MuonEG+Run2016H-PromptReco-v3+AOD
SingleElectron+Run2016H-PromptReco-v2+AOD
SingleMuon+Run2016H-PromptReco-v2+AOD
SinglePhoton+Run2016H-PromptReco-v2+AOD
EOF

for datasetName in `cat $theRND.txt`; do
  echo "Merging dataset: "$datasetName
  ls $INPUTDIR/$datasetName|grep root|grep nero > $theRND_$datasetName.txt;
  mkdir -p $SKIMTEMPDIR/$datasetName;
  for fileName in `cat $theRND_$datasetName.txt`; do
    root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/$datasetName/$fileName\",\"$SKIMTEMPDIR/$datasetName/$fileName\",\"data\"\);
  done
  rm -f $theRND_$datasetName.txt
done

export ALEASTONEPD=0
export msg="hadd -f $SKIMTEMPDIR/data_Run2016H.root "
for datasetName in `cat $theRND.txt`; do
  if [ "$(ls -A $SKIMTEMPDIR/$datasetName)" ] && [[ $datasetName != *"SinglePhoton"* ]]; then
    export msg=$msg" $SKIMTEMPDIR/$datasetName/*.root"
    export ALEASTONEPD=1
  elif [ "$(ls -A $SKIMTEMPDIR/$datasetName)" ] && [[ $datasetName == *"SinglePhoton"* ]]; then
    echo "$SKIMTEMPDIR/$datasetName is photon"
    hadd -f $SKIMTEMPDIR/photon_Run2016H.root $SKIMTEMPDIR/$datasetName/*.root;
    root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$SKIMTEMPDIR/photon_Run2016H.root\",\"$GOODRUNDIR/photon_Run2016H.root\",\"$JSONFILE\"\);
    mv $GOODRUNDIR/photon_Run2016H.root $SKIMDIR/photon_Run2016H.root;
  else
    echo "$SKIMTEMPDIR/$datasetName is empty"
  fi
done

if [ $ALEASTONEPD -eq 1 ]; then
  echo $msg
  $msg

  root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$SKIMTEMPDIR/data_Run2016H.root\",\"$GOODRUNDIR/data_Run2016H.root\",\"$JSONFILE\"\);  
  mv $GOODRUNDIR/data_Run2016H.root $SKIMDIR/data_Run2016H.root;

  root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$SKIMDIR/data_Run2016H.root\",\"$SKIMDIR/met_data_Run2016H.root\",\"data\",2\);
  else
    echo "NO PD TO PROCESS"
fi

rm -f $theRND.txt;
exit;
