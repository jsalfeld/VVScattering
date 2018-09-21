#!/bin/bash

samplesAndDatasets[0]="/WZJJ_EWK_13TeV-madgraph-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"
samplesAndDatasets[1]="/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM"

# ZG sample missing (general)
#WWZ sample missing
#WZZ sample missing
#ZZZ sample missing
#WJe sample missing
#tZq sample missing


#IFS='/' read -r -a array <<< "$samplesAndDatasets[1]"
#echo ${array[1]}


for j in {0..1}; do
	echo ${samplesAndDatasets[j]}
	IFS='/' read -r -a array <<< "${samplesAndDatasets[j]}"
	#echo ${array[1]}
	crab submit -c crab_cfgAuto.py General.requestName=${array[1]} Data.inputDataset=${samplesAndDatasets[j]} Data.inputDBS="global"
	array=0
	
done