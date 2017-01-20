#!/bin.sh

python sendOnBatch.py -i ../test/testNero.py --mc -d mysub/TTToSemilepton -e /store/mc/RunIISummer16MiniAODv2/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1 --put-in=/store/group/phys_higgs/ceballos/test/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8 -q 8nh -n 1000
python sendOnBatch.py -i ../test/testNero.py --mc -d mysub/TTTo2L2Nu -e /store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1 --put-in=/store/group/phys_higgs/ceballos/test/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8 -q 8nh -n 1000

awk '{split($1,a,"/");print a[2]}' ../test/miniaod_samples.txt | awk '{print"python sendOnBatch.py -i ../test/testNero.py --mc -d mysub/"$1" -e /store/mc/RunIISummer16MiniAODv2/"$1"/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1 --put-in=/store/group/phys_higgs/ceballos/test/"$1" -q 1nd -n 1000"}'

python sendOnBatch.py -s -d mysub/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8

python sendOnBatch.py --only-submit -d mysub/HWplusJ_HToWW_M125_13TeV_powheg_pythia8  -q 1nd -j "fail"
