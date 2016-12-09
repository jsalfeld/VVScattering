#!/bin/sh

export theRND=file_skim_batch_$RANDOM;

export INPUTDIR=$1;
export SKIMDIR=$2;
export TYPE=$3

export workDir=$HOME/Condor_proc

find ${INPUTDIR}/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > ${theRND}.txt

for datasetName in `cat ${theRND}.txt`; do
  echo "Merging dataset: "${INPUTDIR}/${datasetName}
  mkdir -p ${SKIMDIR}/${datasetName};
  mkdir -p ${workDir}/test/${SKIMDIR}/${datasetName};
  ls ${INPUTDIR}/${datasetName}|grep root > ${SKIMDIR}/${datasetName}.txt;
  for fileName in `cat ${SKIMDIR}/${datasetName}.txt`; do
    bsub -q 1nh -o ${workDir}/test/${SKIMDIR}/${datasetName}/${fileName}.out -J ${SKIMDIR}_${datasetName}_${fileName} $HOME/releases/CMSSW_8_0_20/src/MitAnalysisRunII/bin/80x/cern/skim_batch.sh ${INPUTDIR} ${SKIMDIR} ${datasetName} $fileName $TYPE
  done
  rm -f ${SKIMDIR}/${datasetName}.txt;
done

rm -f ${theRND}.txt;
