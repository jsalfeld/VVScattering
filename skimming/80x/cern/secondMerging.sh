#!/bin/sh

export JSONFILE=MitAnalysisRunII/json/80x/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt;
export MERGINGDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x;
export OUTPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x;

export ERA=$1;

hadd -f $MERGINGDIR/data_${ERA}.root $MERGINGDIR/DoubleEG-${ERA}-*_????.root $MERGINGDIR/DoubleMuon-${ERA}-*_????.root $MERGINGDIR/MuonEG-${ERA}-*_????.root $MERGINGDIR/SingleMuon-${ERA}-*.root $MERGINGDIR/SingleElectron-${ERA}-*_????.root;
hadd -f $MERGINGDIR/photon_${ERA}.root $MERGINGDIR/SinglePhoton-${ERA}-*_????.root;

root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIR/photon_${ERA}.root\",\"$OUTPUTDIR/photon_${ERA}.root\",\"$JSONFILE\"\);
root -l -q -b MitAnalysisRunII/skimming/80x/makeGoodRunSample.C+\(\"$MERGINGDIR/data_${ERA}.root\",\"$OUTPUTDIR/data_${ERA}.root\",\"$JSONFILE\"\);
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$OUTPUTDIR/data_${ERA}.root\",\"$OUTPUTDIR/met_data_${ERA}.root\",\"data\",2\);
