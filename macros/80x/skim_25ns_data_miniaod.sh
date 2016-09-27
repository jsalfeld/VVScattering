#!/bin/sh

root -l -q -b MitAnalysisRunII/macros/80x/makeGoodRunSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_Run2016B.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016B.root\",\"MitAnalysisRunII/json/80x/Cert_271036-279116_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt\"\)
root -l -q -b MitAnalysisRunII/macros/80x/makeGoodRunSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_Run2016C.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016C.root\",\"MitAnalysisRunII/json/80x/Cert_271036-279116_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt\"\)
root -l -q -b MitAnalysisRunII/macros/80x/makeGoodRunSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_Run2016D.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016D.root\",\"MitAnalysisRunII/json/80x/Cert_271036-279116_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt\"\)
root -l -q -b MitAnalysisRunII/macros/80x/makeGoodRunSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_Run2016E.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016E.root\",\"MitAnalysisRunII/json/80x/Cert_271036-279116_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt\"\)
root -l -q -b MitAnalysisRunII/macros/80x/makeGoodRunSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_Run2016F.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016F.root\",\"MitAnalysisRunII/json/80x/Cert_271036-279116_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt\"\)

root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016B.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/met_data_Run2016B_skim.root\",\"data\",2\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016C.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/met_data_Run2016C_skim.root\",\"data\",2\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016D.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/met_data_Run2016D_skim.root\",\"data\",2\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016E.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/met_data_Run2016E_skim.root\",\"data\",2\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016F.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/met_data_Run2016F_skim.root\",\"data\",2\)

root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016B.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_Run2016B_skim.root\",\"data\",0\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016C.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_Run2016C_skim.root\",\"data\",0\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016D.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_Run2016D_skim.root\",\"data\",0\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016E.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_Run2016E_skim.root\",\"data\",0\)
root -l -q -b MitAnalysisRunII/macros/80x/makeOneSkimSample.C+\(\"/scratch5/dhsu/ntuples_goodrun_80x/data_miniAOD_goodrun_Run2016F.root\",\"/scratch5/dhsu/ntuples_goodrun_80x/data_Run2016F_skim.root\",\"data\",0\)
