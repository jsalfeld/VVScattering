#!/bin/sh

echo "Arguments: $*"

export INPUTDIR=$1;
export SKIMDIR=$2;
export datasetName=$3
export fileName=$4
export TYPE=$5

/cvmfs/cms.cern.ch/cmsset_default.sh;
alias eosmount='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount'
if [[ ! -e ~/eos/cms/store ]]; then
  eosmount ~/eos;
fi
cd ~/releases/CMSSW_8_0_20/src/;
alias cmsenv='eval `scramv1 runtime -sh`'

if [[ ! -e $CMSSW_BASE/src ]]; then
   cmsenv
fi

if [[ -e ~/eos/cms/store ]] &&  [[ -e $CMSSW_BASE/src ]]; then
   root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"${INPUTDIR}/${datasetName}/${fileName}\",\"${SKIMDIR}/${datasetName}/${fileName}\",\"${TYPE}\",0,0,0\);
else
   echo "INITIALIZATION FAILED";
fi
