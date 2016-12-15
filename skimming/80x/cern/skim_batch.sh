#!/bin/sh

echo "Arguments: $*"

export INPUTDIR=$1;
export SKIMDIR=$2;
export datasetName=$3
export fileName=$4
export TYPE=$5

export thePWD=$PWD;

/cvmfs/cms.cern.ch/cmsset_default.sh;
alias eosmount='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount'
#if [[ ! -e ${thePWD}/eoslink2/cms/store ]]; then
#  eosmount ${thePWD}/eoslink2;
#fi
cd ~/releases/CMSSW_8_0_20/src/;
alias cmsenv='eval `scramv1 runtime -sh`'
if [[ ! -e $CMSSW_BASE/src ]]; then
   cmsenv
fi

#if [[ -e ~/eos/cms/store ]] &&  [[ -e $CMSSW_BASE/src ]]; then
if [[ -e $CMSSW_BASE/src ]]; then
   #root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"${thePWD}/eos/${INPUTDIR}/${datasetName}/${fileName}\",\"${thePWD}/eos/${SKIMDIR}/${datasetName}/${fileName}\",\"${TYPE}\",0,0,0\);
   root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"${INPUTDIR}/${datasetName}/${fileName}\",\"${SKIMDIR}/${datasetName}/${fileName}\",\"${TYPE}\",0,0,0\);
   #mkdir -p ${SKIMDIR}/${datasetName};
   #touch ${SKIMDIR}/${datasetName}/${fileName}
   #eos rm /eos/${SKIMDIR}/${datasetName}/${fileName};
   #xrdcp ${SKIMDIR}/${datasetName}/${fileName} root://eoscms.cern.ch//eos/${SKIMDIR}/${datasetName}/${fileName};
   #rm -f ${SKIMDIR}/${datasetName}/${fileName};
else
   echo "INITIALIZATION FAILED";
fi

#pkill -u ceballos eosfsd;
