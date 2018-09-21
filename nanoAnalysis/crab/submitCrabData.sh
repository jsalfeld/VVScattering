#!/bin/bash

samplesAndDatasets[0]="/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017B-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[1]="/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017C-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[2]="/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017D-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[3]="/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017E-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[4]="/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017F-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[5]="/DoubleEG/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017F-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[6]="/DoubleEG/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017B-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[7]="/DoubleEG/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017C-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[8]="/DoubleEG/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017D-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[9]="/DoubleEG/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017E-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER"
samplesAndDatasets[10]="/SingleMuon/arizzi-RunII2017ReReco17Nov17-94X-Nano01-e70630e8aef2c186cd650f6150c31168/USER"
samplesAndDatasets[11]="/SingleElectron/arizzi-RunII2017ReReco17Nov17-94X-Nano01-e70630e8aef2c186cd650f6150c31168/USER"
samplesAndDatasets[12]="/DoubleMuon/arizzi-RunII2017ReReco17Nov17-94X-Nano01-e70630e8aef2c186cd650f6150c31168/USER"
samplesAndDatasets[13]="/DoubleEG/arizzi-RunII2017ReReco17Nov17-94X-Nano01-e70630e8aef2c186cd650f6150c31168/USER"


# ZG sample missing (general)
#WWZ sample missing
#WZZ sample missing
#ZZZ sample missing
#WJe sample missing
#tZq sample missing


#IFS='/' read -r -a array <<< "$samplesAndDatasets[1]"
#echo ${array[1]}


for j in {0..9}; do
	echo ${samplesAndDatasets[j]}
	IFS='/' read -r -a array <<< "${samplesAndDatasets[j]}"
	#echo ${array[1]}${array[2]:42:8}
	crab submit -c crab_cfgAutoData.py General.requestName=${array[1]}${array[2]:42:8} Data.inputDataset=${samplesAndDatasets[j]} Data.outputDatasetTag=${array[1]}${array[2]:42:8}
	array=0
	
done