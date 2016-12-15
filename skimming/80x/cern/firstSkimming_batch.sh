#!/bin/sh

export theRND=file_skim_batch_$RANDOM;

export INPUTDIR=$1;
export SKIMDIR=$2;
export TYPE=$3

export workDir=/afs/cern.ch/work/c/ceballos/test

alias eosmount='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount'
if [[ ! -e ~/eoslink2/cms/store ]]; then
  eosmount ~/eoslink2;
fi

find ~/eoslink2/${INPUTDIR}/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > ${theRND}.txt

for datasetName in `cat ${theRND}.txt`; do
  echo "Merging dataset: "${INPUTDIR}/${datasetName}
  mkdir -p ~/eoslink2/${SKIMDIR}/${datasetName};
  mkdir -p ${workDir}/${SKIMDIR}/${datasetName};
  ls ~/eoslink2/${INPUTDIR}/${datasetName}|grep root > ~/eoslink2/${SKIMDIR}/${datasetName}.txt;
  for fileName in `cat ~/eoslink2/${SKIMDIR}/${datasetName}.txt`; do
    bsub -q 1nh -o ${workDir}/${SKIMDIR}/${datasetName}/${fileName}.out -J ${SKIMDIR}_${datasetName}_${fileName} $HOME/releases/CMSSW_8_0_20/src/MitAnalysisRunII/skimming/80x/cern/skim_batch.sh ${INPUTDIR} ${SKIMDIR} ${datasetName} $fileName $TYPE
  done
  rm -f ${SKIMDIR}/${datasetName}.txt;
done

rm -f ${theRND}.txt;
