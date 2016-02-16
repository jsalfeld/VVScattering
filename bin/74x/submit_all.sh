#!/bin/sh

./submit.sh bambuToNero.py 01 DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 02 ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 03 TTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 04 WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM
./submit.sh bambuToNero.py 05 WWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 06 WW_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM
./submit.sh bambuToNero.py 07 WZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 08 ZZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM
./submit.sh bambuToNero.py 11 DoubleEG+Run2015B-PromptReco-v1+AOD        "--data"
./submit.sh bambuToNero.py 12 DoubleMuon+Run2015B-PromptReco-v1+AOD      "--data"
./submit.sh bambuToNero.py 13 SingleElectron+Run2015B-PromptReco-v1+AOD  "--data"
./submit.sh bambuToNero.py 14 SingleMuon+Run2015B-PromptReco-v1+AOD      "--data"
