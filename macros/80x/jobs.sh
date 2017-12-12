#!/bin/bash

suffix=$1
workDirLocal=$2
outputDir=$3

workDirRemote=`pwd`

cd ${workDirLocal}
eval `scramv1 runtime -sh`
#cd $workDirRemote

#echo $suffix
#echo (($suffix))

#cp ${workDirLocal}/* .

#eos cp $workDirRemote/*zmass*.root $outputDir
declare -i y=$1
echo $y
root -l -b -q 'MitAnalysisRunII/macros/80x/wzScatteringAnalysisBatch.C(76,106,true,"default","default","./testBatch/",false,'$y')'