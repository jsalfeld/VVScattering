#!/bin/sh

export JSONFILE=MitAnalysisRunII/json/80x/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt;
export MERGINGDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
#export OUTPUTDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x;
export MERGINGDIREOS=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
export OUTPUTDIR=root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x;

export WORKDIR=/afs/cern.ch/work/c/ceballos/test

export ERA=$1;

hadd -O -k -f $MERGINGDIR/MET_${ERA}.root $MERGINGDIR/MET*${ERA}*_????.root $MERGINGDIR/MET*${ERA}*failed.root;
hadd -O -k -f $MERGINGDIR/photon_${ERA}.root $MERGINGDIR/SinglePhoton*${ERA}*_????.root $MERGINGDIR/SinglePhoton*${ERA}*failed.root;
hadd -O -k -f $MERGINGDIR/data_${ERA}.root $MERGINGDIR/SingleElectron*${ERA}*_????.root $MERGINGDIR/DoubleEG*${ERA}*_????.root $MERGINGDIR/DoubleMuon*${ERA}*_????.root $MERGINGDIR/MuonEG*${ERA}*_????.root $MERGINGDIR/SingleMuon*${ERA}*_????.root $MERGINGDIR/SingleElectron*${ERA}*failed.root $MERGINGDIR/DoubleEG*${ERA}*failed.root $MERGINGDIR/DoubleMuon*${ERA}*failed.root $MERGINGDIR/MuonEG*${ERA}*failed.root $MERGINGDIR/SingleMuon*${ERA}*failed.root;

#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/MET_${ERA}.root\",\"$OUTPUTDIR/MET_${ERA}.root\",\"$JSONFILE\"\);
#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/photon_${ERA}.root\",\"$OUTPUTDIR/photon_${ERA}.root\",\"$JSONFILE\"\);
#root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/data_${ERA}.root\",\"$OUTPUTDIR/data_${ERA}.root\",\"$JSONFILE\"\);
#root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_${ERA}.root\",\"$OUTPUTDIR/met_data_${ERA}.root\",\"data\",2\);

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/MET_${ERA}.root\",\"$WORKDIR/MET_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
xrdcp $WORKDIR/MET_GOODRUN_${ERA}.root $OUTPUTDIR/MET_GOODRUN_${ERA}.root;
rm -f $WORKDIR/MET_GOODRUN_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/photon_${ERA}.root\",\"$WORKDIR/photon_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
xrdcp $WORKDIR/photon_GOODRUN_${ERA}.root $OUTPUTDIR/photon_GOODRUN_${ERA}.root;
rm -f $WORKDIR/photon_GOODRUN_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIREOS/data_${ERA}.root\",\"$WORKDIR/data_GOODRUN_${ERA}.root\",\"$JSONFILE\"\);
xrdcp $WORKDIR/data_GOODRUN_${ERA}.root $OUTPUTDIR/data_GOODRUN_${ERA}.root;
rm -f $WORKDIR/data_GOODRUN_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/MET_GOODRUN_${ERA}.root\",\"$WORKDIR/MET_${ERA}.root\",\"data\",0\);
xrdcp $WORKDIR/MET_${ERA}.root $OUTPUTDIR/MET_${ERA}.root;
rm -f $WORKDIR/MET_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/photon_GOODRUN_${ERA}.root\",\"$WORKDIR/photon_${ERA}.root\",\"data\",0\);
xrdcp $WORKDIR/photon_${ERA}.root $OUTPUTDIR/photon_${ERA}.root;
rm -f $WORKDIR/photon_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_GOODRUN_${ERA}.root\",\"$WORKDIR/data_${ERA}.root\",\"data\",0\);
xrdcp $WORKDIR/data_${ERA}.root $OUTPUTDIR/data_${ERA}.root;
rm -f $WORKDIR/data_${ERA}.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_GOODRUN_${ERA}.root\",\"$WORKDIR/met_data_${ERA}.root\",\"data\",2\);
xrdcp $WORKDIR/met_data_${ERA}.root $OUTPUTDIR/met_data_${ERA}.root;
rm -f $WORKDIR/met_data_${ERA}.root;
