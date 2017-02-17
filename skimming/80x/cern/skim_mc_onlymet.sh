#!/bin/sh

export INPUTDIR=/eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet;
export OUTPUTDIR=/eos/cms/store/caf/user/ceballos/Nero/output_80x;

#ls $INPUTDIR|grep MINIAODSIM|awk '{printf(" root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\\(\\\"$INPUTDIR/%s\\\",\\\"$OUTPUTDIR/${PREFIX}%s\\\",\\\"data\\\",$OPTION\\)\n",$1,$1)}'

export OPTION=2;

export PREFIX="";
if [ $OPTION == 1 ]
then
  export PREFIX="qcd_";
elif [ $OPTION == 2 ]
then
  export PREFIX="met_";
elif [ $OPTION == 3 ]
then
  export PREFIX="zmet_";
elif [ $OPTION == 4 ]
then
  export PREFIX="pho_";
fi

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_mllmin01_13TeV-powheg-pythia8.root\",\"wz3ln1_powheg\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"zll50_emu\",$OPTION\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt050To100\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt100To250\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt250To400\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt400To650\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt650ToInf\",$OPTION\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy1050_1j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy1050_2j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy1050_3j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy1050_4j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy50_1j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy50_2j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy50_3j\",$OPTION\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"dy50_4j\",$OPTION\)
