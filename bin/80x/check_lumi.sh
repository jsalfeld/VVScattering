#!/bin/sh

export JSON=Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt;
export MAIN=$CMSSW_BASE/src;
cd $MAIN;

echo "mueg"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_mueg.txt /scratch5/ceballos/ntuples_noweights_80x/MuonEG+Run2016B-PromptReco-v2+AOD.root;

echo "dmu"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_dmu.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleMuon+Run2016B-PromptReco-v2+AOD.root;

echo "smu"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_smu.txt /scratch5/ceballos/ntuples_noweights_80x/SingleMuon+Run2016B-PromptReco-v2+AOD.root;

echo "del"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_del.txt /scratch5/ceballos/ntuples_noweights_80x/DoubleEG+Run2016B-PromptReco-v2+AOD.root;

echo "sel"
$MAIN/MitAnalysisRunII/bin/makeJson.py -m $MAIN/MitAnalysisRunII/json/80x/$JSON -o lumis_r2016b_sel.txt /scratch5/ceballos/ntuples_noweights_80x/SingleElectron+Run2016B-PromptReco-v2+AOD.root;

cd -;
