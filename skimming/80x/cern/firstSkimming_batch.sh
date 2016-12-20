#!/bin/sh

export theRND=file_skim_batch_$RANDOM;

export INPUTDIR=$1;
export SKIMDIR=$2;
export TYPE=$3

export workDir=/afs/cern.ch/work/c/ceballos/test

alias eosmount='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount'
if [[ ! -e ~/eos/cms/store ]]; then
  eosmount ~/eos;
fi

if [[ ! -e $CMSSW_BASE/src ]]; then
   echo 'begin setting CMSSW'
   cd ~/releases/CMSSW_8_0_24_patch1/src/;
   /cvmfs/cms.cern.ch/cmsset_default.sh;
   eval `scramv1 runtime -sh`
   echo 'end setting CMSSW'
else
   echo 'CMSSW found: '$CMSSW_BASE/src
fi

cd ~/releases/CMSSW_8_0_24_patch1/src/;

find ~/eos/${INPUTDIR}/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > ${theRND}.txt

for datasetName in `cat ${theRND}.txt`; do
  echo "Merging dataset: "${INPUTDIR}/${datasetName}
  mkdir -p ~/eos/${SKIMDIR}/${datasetName};
  mkdir -p ${workDir}/${SKIMDIR}/${datasetName};
  ls ~/eos/${INPUTDIR}/${datasetName}|grep root > ~/eos/${SKIMDIR}/${datasetName}.txt;
  for fileName in `cat ~/eos/${SKIMDIR}/${datasetName}.txt`; do
    bsub -q 1nh -o ${workDir}/${SKIMDIR}/${datasetName}/${fileName}.out -J ${SKIMDIR}_${datasetName}_${fileName} $CMSSW_BASE/src/MitAnalysisRunII/skimming/80x/cern/skim_batch.sh ${INPUTDIR} ${SKIMDIR} ${datasetName} $fileName $TYPE
    #$CMSSW_BASE/src/MitAnalysisRunII/skimming/80x/cern/skim_batch.sh ${INPUTDIR} ${SKIMDIR} ${datasetName} $fileName $TYPE
  done
  rm -f ${SKIMDIR}/${datasetName}.txt;
done

rm -f ${theRND}.txt;
pkill -u ceballos eosfsd;
