#!/bin/sh

export JSONFILE=MitAnalysisRunII/json/80x/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt;
export MERGINGDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
#export OUTPUTDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x;
export MERGINGDIREOS=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
export OUTPUTDIR=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x;

export WORKDIR=/afs/cern.ch/work/c/ceballos/test

export ERA=$1;

hadd -O -k -f $WORKDIR/data1_${ERA}.root $MERGINGDIR/SingleElectron*${ERA}*_????.root $MERGINGDIR/DoubleEG*${ERA}*_????.root $MERGINGDIR/SingleElectron*${ERA}*failed.root $MERGINGDIR/DoubleEG*${ERA}*failed.root 
hadd -O -k -f $WORKDIR/data2_${ERA}.root $MERGINGDIR/DoubleMuon*${ERA}*_????.root $MERGINGDIR/SingleMuon*${ERA}*_????.root $MERGINGDIR/DoubleMuon*${ERA}*failed.root $MERGINGDIR/SingleMuon*${ERA}*failed.root;
hadd -O -k -f $WORKDIR/data3_${ERA}.root $MERGINGDIR/MuonEG*${ERA}*_????.root $MERGINGDIR/MuonEG*${ERA}*failed.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$WORKDIR/data1_${ERA}.root\",\"$WORKDIR/data1_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$WORKDIR/data2_${ERA}.root\",\"$WORKDIR/data2_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$WORKDIR/data3_${ERA}.root\",\"$WORKDIR/data3_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
rm -f $WORKDIR/data?_${ERA}.root;

hadd -O -k -f $WORKDIR/data_${ERA}.root $WORKDIR/data?_GOODRUN_${ERA}.root
xrdcp $WORKDIR/data_${ERA}.root $MERGINGDIREOS/data_${ERA}.root;
rm -f $WORKDIR/dat*_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/data_${ERA}.root\",\"$WORKDIR/data_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
xrdcp $WORKDIR/data_GOODRUN_${ERA}.root $OUTPUTDIR/data_GOODRUN_${ERA}.root;
rm -f $WORKDIR/data_GOODRUN_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_GOODRUN_${ERA}.root\",\"$WORKDIR/data_${ERA}.root\",\"data\",0\);
xrdcp $WORKDIR/data_${ERA}.root $OUTPUTDIR/data_${ERA}.root;
rm -f $WORKDIR/data_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_GOODRUN_${ERA}.root\",\"$WORKDIR/met_data_${ERA}.root\",\"data\",2\);
xrdcp $WORKDIR/met_data_${ERA}.root $OUTPUTDIR/met_data_${ERA}.root;
rm -f $WORKDIR/met_data_${ERA}.root;
