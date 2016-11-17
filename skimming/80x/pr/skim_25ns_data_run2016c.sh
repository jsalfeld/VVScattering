#!/bin/sh

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/MuonEG+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/temp/MuonEG+Run2016C-PromptReco-v2+AOD.root\",\"data\"\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/DoubleMuon+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/temp/DoubleMuon+Run2016C-PromptReco-v2+AOD.root\",\"data\"\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/SingleMuon+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/temp/SingleMuon+Run2016C-PromptReco-v2+AOD.root\",\"data\"\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/DoubleEG+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/temp/DoubleEG+Run2016C-PromptReco-v2+AOD.root\",\"data\"\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/SingleElectron+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/temp/SingleElectron+Run2016C-PromptReco-v2+AOD.root\",\"data\"\)

hadd -f /scratch/ceballos/ntuples_goodrun_80x/temp/data_Run2016C.root /scratch/ceballos/ntuples_goodrun_80x/temp/*Run2016C*AOD*.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/scratch/ceballos/ntuples_goodrun_80x/temp/data_Run2016C.root\",\"/scratch/ceballos/ntuples_goodrun_80x/data_Run2016C.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/SinglePhoton+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/SinglePhoton+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_80x/data_Run2016C.root\",\"/scratch/ceballos/ntuples_weightsDA_80x/data_Run2016C.root\",\"data\",0\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_80x/data_Run2016C.root\",\"/scratch/ceballos/ntuples_weightsDA_80x/met_data_Run2016C.root\",\"data\",2\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_80x/SinglePhoton+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_weightsDA_80x/SinglePhoton+Run2016C-PromptReco-v2+AOD.root\",\"data\",0\)

#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/MuonEG+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/MuonEG+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/DoubleMuon+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/DoubleMuon+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
###root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/SingleMuon+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/SingleMuon+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/DoubleEG+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/DoubleEG+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
###root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"/data/t3home000/ceballos/merging_ntuples/SingleElectron+Run2016C-PromptReco-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_80x/SingleElectron+Run2016C-PromptReco-v2+AOD.root\",\"MitAnalysisRunII/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\"\)
