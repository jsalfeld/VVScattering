#!/bin/sh

export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/merging_80x;
export OUTPUTDIR=/afs/cern.ch/work/c/ceballos/test;

#ls $INPUTDIR|grep MINIAODSIM|awk '{printf(" root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\\(\\\"$INPUTDIR/%s\\\",\\\"$OUTPUTDIR/${PREFIX}%s\\\",\\\"data\\\",$1\\)\n",$1,$1)}'

export PREFIX="";
if [ $1 == 1 ]
then
  export PREFIX="qcd_";
elif [ $1 == 2 ]
then
  export PREFIX="met_";
elif [ $1 == 3 ]
then
  export PREFIX="zmet_";
elif [ $1 == 4 ]
then
  export PREFIX="pho_";
fi

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_mllmin01_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_mllmin01_13TeV-powheg-pythia8.root\",\"wz3ln1_powheg\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root\",\"tt\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"zll50_emu\",$1\)
