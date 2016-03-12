#!/bin/sh

export JSON=Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt;
export MAIN=$CMSSW_BASE/src;
cd $MAIN;

echo "mueg"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_mueg.txt /scratch5/ceballos/ntuples_noweights_76x/MuonEG+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_mueg.txt /scratch5/ceballos/ntuples_noweights_76x/MuonEG+Run2015D-16Dec2015-v1+AOD.root

echo "dmu"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_dmu.txt /scratch5/ceballos/ntuples_noweights_76x/DoubleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_dmu.txt /scratch5/ceballos/ntuples_noweights_76x/DoubleMuon+Run2015D-16Dec2015-v1+AOD.root

echo "smu"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_smu.txt /scratch5/ceballos/ntuples_noweights_76x/SingleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_smu.txt /scratch5/ceballos/ntuples_noweights_76x/SingleMuon+Run2015D-16Dec2015-v1+AOD.root

echo "del"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_del.txt /scratch5/ceballos/ntuples_noweights_76x/DoubleEG+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_del.txt /scratch5/ceballos/ntuples_noweights_76x/DoubleEG+Run2015D-16Dec2015-v2+AOD.root

echo "sel"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_sel.txt /scratch5/ceballos/ntuples_noweights_76x/SingleElectron+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_sel.txt /scratch5/ceballos/ntuples_noweights_76x/SingleElectron+Run2015D-16Dec2015-v1+AOD.root

echo "spho"
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015c_spho.txt /scratch5/ceballos/ntuples_noweights_76x/SinglePhoton+Run2015C_25ns-16Dec2015-v1+AOD.root
$MAIN/MonoX/common/makeJson.py -m $MAIN/MitAnalysisRunII/json/76x/$JSON -o lumis_r2015d_spho.txt /scratch5/ceballos/ntuples_noweights_76x/SinglePhoton+Run2015D-16Dec2015-v1+AOD.root

cd -;
