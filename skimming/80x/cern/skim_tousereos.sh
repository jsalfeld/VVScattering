#!/bin/sh

export INPUTDIR=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/mc;
export OUTPUTDIR=/eos/user/c/ceballos/Nero/output_80x;

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

if [ $2 == 0 ]
then
root -l -q -b VVScattering/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTTo2L2Nu_13TeV-powheg.root\",\"tt2l\",$1\)
elif [ $2 == 1 ]
then
root -l -q -b VVScattering/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTToSemiLeptonic_13TeV-powheg.root\",\"ttqql\",$1\)
fi
