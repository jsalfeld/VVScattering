#!/bin/sh

export JSON=Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON.txt;
export MAIN=$CMSSW_BASE/src;
cd $MAIN;

echo "mueg-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_mueg.txt /scratch5/ceballos/ntuples_noweights_80x/MuonEG+Run2016B-PromptReco-v2+AOD*.root;

echo "dmu-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_dmu.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleMuon+Run2016B-PromptReco-v2+AOD*.root;

echo "smu-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_smu.txt /scratch5/ceballos/ntuples_noweights_80x/SingleMuon+Run2016B-PromptReco-v2+AOD*.root;

echo "del-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_del.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleEG+Run2016B-PromptReco-v2+AOD*.root;

echo "sel-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_sel.txt /scratch5/ceballos/ntuples_noweights_80x/SingleElectron+Run2016B-PromptReco-v2+AOD*.root;

echo "pho-Run2016B"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_pho.txt /scratch5/ceballos/ntuples_noweights_80x/SinglePhoton+Run2016B-PromptReco-v2+AOD*.root;

echo "mueg-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_mueg.txt /scratch5/ceballos/ntuples_noweights_80x/MuonEG+Run2016C-PromptReco-v2+AOD*.root;

echo "dmu-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_dmu.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleMuon+Run2016C-PromptReco-v2+AOD*.root;

echo "smu-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_smu.txt /scratch5/ceballos/ntuples_noweights_80x/SingleMuon+Run2016C-PromptReco-v2+AOD*.root;

echo "del-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_del.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleEG+Run2016C-PromptReco-v2+AOD*.root;

echo "sel-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_sel.txt /scratch5/ceballos/ntuples_noweights_80x/SingleElectron+Run2016C-PromptReco-v2+AOD*.root;

echo "pho-Run2016C"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016c_pho.txt /scratch5/ceballos/ntuples_noweights_80x/SinglePhoton+Run2016C-PromptReco-v2+AOD*.root;

cd -;
