#!/bin/sh

export JSON=Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt;
export MAIN=~ceballos/cms/cmssw/042/CMSSW_7_4_6/src;
cd $MAIN;

echo "mueg"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_mueg.txt /scratch5/ceballos/ntuples_noweights/MuonEG+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_mueg.txt /scratch5/ceballos/ntuples_noweights/MuonEG+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_mueg.txt /scratch5/ceballos/ntuples_noweights/MuonEG+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_mueg.txt /scratch5/ceballos/ntuples_noweights/MuonEG+Run2015D-PromptReco-v4+AOD.root

echo "dmu"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_dmu.txt /scratch5/ceballos/ntuples_noweights/DoubleMuon+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_dmu.txt /scratch5/ceballos/ntuples_noweights/DoubleMuon+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_dmu.txt /scratch5/ceballos/ntuples_noweights/DoubleMuon+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_dmu.txt /scratch5/ceballos/ntuples_noweights/DoubleMuon+Run2015D-PromptReco-v4+AOD.root

echo "smu"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_smu.txt /scratch5/ceballos/ntuples_noweights/SingleMuon+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_smu.txt /scratch5/ceballos/ntuples_noweights/SingleMuon+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_smu.txt /scratch5/ceballos/ntuples_noweights/SingleMuon+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_smu.txt /scratch5/ceballos/ntuples_noweights/SingleMuon+Run2015D-PromptReco-v4+AOD.root

echo "del"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_del.txt /scratch5/ceballos/ntuples_noweights/DoubleEG+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_del.txt /scratch5/ceballos/ntuples_noweights/DoubleEG+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_del.txt /scratch5/ceballos/ntuples_noweights/DoubleEG+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_del.txt /scratch5/ceballos/ntuples_noweights/DoubleEG+Run2015D-PromptReco-v4+AOD.root

echo "sel"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_sel.txt /scratch5/ceballos/ntuples_noweights/SingleElectron+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_sel.txt /scratch5/ceballos/ntuples_noweights/SingleElectron+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_sel.txt /scratch5/ceballos/ntuples_noweights/SingleElectron+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_sel.txt /scratch5/ceballos/ntuples_noweights/SingleElectron+Run2015D-PromptReco-v4+AOD.root

echo "spho"
###$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015b1_spho.txt /scratch5/ceballos/ntuples_noweights/SinglePhoton+Run2015B-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015c1_spho.txt /scratch5/ceballos/ntuples_noweights/SinglePhoton+Run2015C-PromptReco-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d3_spho.txt /scratch5/ceballos/ntuples_noweights/SinglePhoton+Run2015D-PromptReco-v3+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/$JSON -o lumis_r2015d4_spho.txt /scratch5/ceballos/ntuples_noweights/SinglePhoton+Run2015D-PromptReco-v4+AOD.root

cd -;
