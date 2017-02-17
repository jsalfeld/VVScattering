#!/bin/sh

echo "Arguments: $*"

export INPUTDIR=$1;
export SKIMDIR=$2;
export datasetName=$3
export fileName=$4
export TYPE=$5

export thePWD=$PWD;

#alias eosmount='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount'
#if [[ ! -e ${thePWD}/eos/cms/store ]]; then
#  eosmount ${thePWD}/eos;
#fi
#if [[ ! -e $CMSSW_BASE/src ]]; then
   echo 'begin setting CMSSW'
   cd ~/releases/CMSSW_8_0_26_patch1/src/;
   /cvmfs/cms.cern.ch/cmsset_default.sh;
   eval `scramv1 runtime -sh`
   echo 'end setting CMSSW'
#else
#   echo 'CMSSW found: '$CMSSW_BASE/src
#   exit;
#fi

cd ~/releases/CMSSW_8_0_26_patch1/src/;

echo 'Working folder: '${thePWD}

#if [[ -e /eos/cms/store ]] &&  [[ -e $CMSSW_BASE/src ]]; then
if [[ -e $CMSSW_BASE/src ]]; then
   ###root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"root://eoscms.cern.ch//eos/${INPUTDIR}/${datasetName}/${fileName}\",\"${thePWD}/${SKIMDIR}/${datasetName}/${fileName}\",\"${TYPE}\",0,0,0\);
   mkdir -p ${thePWD}/eos/${SKIMDIR}/${datasetName};
   root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"root://eoscms.cern.ch//eos/${INPUTDIR}/${datasetName}/${fileName}\",\"${thePWD}/eos/${SKIMDIR}/${datasetName}/${fileName}\",\"${TYPE}\",0,0,0\);
   eos rm /eos/${SKIMDIR}/${datasetName}/${fileName};
   xrdcp ${thePWD}/eos/${SKIMDIR}/${datasetName}/${fileName} root://eoscms.cern.ch//eos/${SKIMDIR}/${datasetName}/${fileName};
   rm -f ${thePWD}/eos/${SKIMDIR}/${datasetName}/${fileName};
   rm -f core.*;
else
   echo "INITIALIZATION FAILED";
fi

#pkill -u ceballos eosfsd;
