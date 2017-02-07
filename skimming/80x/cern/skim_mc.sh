#!/bin/sh

export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/merging_80x;
export OUTPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/output_80x;

#ls $INPUTDIR|grep MINIAODSIM|awk '{printf(" root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\\(\\\"$INPUTDIR/%s\\\",\\\"$OUTPUTDIR/${PREFIX}%s\\\",\\\"data\\\",$1\\)\n",$1,$1)}'

export PREFIX="";
if [ $1 == 1 ]
then
  export PREFIX="qcd_";
elif [ $1 == 2 ]
then
  export PREFIX="met_";
elif [ $1 == 3 ]
then
  export PREFIX="zmet_";
elif [ $1 == 4 ]
then
  export PREFIX="pho_";
fi

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"zll1050\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"zll50\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt050To100\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt100To250\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt250To400\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt400To650\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"dyll_pt650ToInf\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root\",\"ggzz2mu2nu\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root\",\"ggzz2e2nu\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root\",\"ggzz2e2mu\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root\",\"ggzz2e2tau\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root\",\"ggzz2mu2tau\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root\",\"ggzz4e\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root\",\"ggzz4mu\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root\",\"ggzz4tau\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluHToTauTau_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluHToTauTau_M125_13TeV_powheg_pythia8.root\",\"gghtautau125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_herwigpp_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_herwigpp.root\",\"gghwwlnln125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root\",\"gghwwlnln125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root\",\"gghzz4l125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/VBFHToTauTau_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}VBFHToTauTau_M125_13TeV_powheg_pythia8.root\",\"vbfhtautau125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8.root\",\"vbfhwwlnln125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root\",\"vbfhzz4l125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/HWminusJ_HToWW_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}HWminusJ_HToWW_M125_13TeV_powheg_pythia8.root\",\"wmhww125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/HWplusJ_HToWW_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}HWplusJ_HToWW_M125_13TeV_powheg_pythia8.root\",\"wphww125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/HZJ_HToWW_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}HZJ_HToWW_M125_13TeV_powheg_pythia8.root\",\"zhww125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root\",\"wmhzz4l125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root\",\"wphzz4l125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WminusHToTauTau_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WminusHToTauTau_M125_13TeV_powheg_pythia8.root\",\"wmhtautau125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WplusHToTauTau_M125_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WplusHToTauTau_M125_13TeV_powheg_pythia8.root\",\"wphtautau125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluWWTo2L2Nu_HInt_MCFM_13TeV_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluWWTo2L2Nu_HInt_MCFM_13TeV.root\",\"ggwwlnln_int\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/GluGluWWTo2L2Nu_MCFM_13TeV_0000.root\",\"$OUTPUTDIR/${PREFIX}GluGluWWTo2L2Nu_MCFM_13TeV.root\",\"ggwwlnln\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_0000.root\",\"$OUTPUTDIR/${PREFIX}ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\",\"tbarw\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_0000.root\",\"$OUTPUTDIR/${PREFIX}ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root\",\"tw\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root\",\"vhnob125\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root\",\"tthnob125\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root\",\"ttg\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root\",\"ttwln\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root\",\"ttwqq\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"ttzlnln\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"ttzqq\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/tZq_ll_4f_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}tZq_ll_4f_13TeV-amcatnlo-pythia8.root\",\"tzqll\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph.root\",\"wlng130\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph.root\",\"zllg130\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph.root\",\"znng40130\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph.root\",\"znng130\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WGstarToLNuEE_012Jets_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}WGstarToLNuEE_012Jets_13TeV-madgraph.root\",\"wgslee\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WGstarToLNuMuMu_012Jets_13TeV-madgraph_0000.root\",\"$OUTPUTDIR/${PREFIX}WGstarToLNuMuMu_012Jets_13TeV-madgraph.root\",\"wgslmm\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1_0000.root\",\"$OUTPUTDIR/${PREFIX}WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1.root\",\"wpwp_ewk_pwg\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1_0000.root\",\"$OUTPUTDIR/${PREFIX}WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1.root\",\"wmwm_ewk_pwg\",$1\)
###root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WLLJJToLNu_M-60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wlljj60_ewk_qcd\",$1\)
###root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WLLJJToLNu_M-4to60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WLLJJToLNu_M-4to60_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wlljj4_60_ewk_qcd\",$1\)
###root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8.root\",\"wlljj60_ewk\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wpwp_ewk_qcd\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wpwp_ewk\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wpwp_qcd\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WGJJToLNu_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WGJJToLNu_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root\",\"wlgjj\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWJJToLNuLNu_EWK_QCD_noTop_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWJJToLNuLNu_EWK_QCD_noTop_13TeV-madgraph-pythia8.root\",\"wpwmjj_ewk_qcd_notop\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWJJToLNuLNu_EWK_noTop_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWJJToLNuLNu_EWK_noTop_13TeV-madgraph-pythia8.root\",\"wpwmjj_ewk_notop\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWJJToLNuLNu_QCD_noTop_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWJJToLNuLNu_QCD_noTop_13TeV-madgraph-pythia8.root\",\"wpwmjj_qcd_notop\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZJJTo4L_EWK_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZJJTo4L_EWK_13TeV-madgraph-pythia8.root\",\"zzjj_ewk\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZJJTo4L_QCD_13TeV-madgraph-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZJJTo4L_QCD_13TeV-madgraph-pythia8.root\",\"zzjj_qcd\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8.root\",\"wlljj_ewk\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_0Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_0Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"wz3l_0j\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_1Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_1Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"wz3l_1j\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_2Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_2Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"wz3l_2j\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_3Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_3Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root\",\"wz3l_3j\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWTo2L2Nu_DoubleScattering_13TeV-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWTo2L2Nu_DoubleScattering_13TeV-pythia8.root\",\"wwlnln_dps\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"wlg\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"zllg\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root\",\"wln\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWTo2L2Nu_13TeV-powheg_0000.root\",\"$OUTPUTDIR/${PREFIX}WWTo2L2Nu_13TeV-powheg.root\",\"wwlnln\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWTo2L2Nu_13TeV-powheg-herwigpp_0000.root\",\"$OUTPUTDIR/${PREFIX}WWTo2L2Nu_13TeV-powheg-herwigpp.root\",\"wwlnln\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"www\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"wwz\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"wzz\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root\",\"zzz\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root\",\"wz1l1n2q4\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8.root\",\"wz1l3n4\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root\",\"wz2q2l4\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root\",\"wz3ln4_powheg\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZTo2L2Nu_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZTo2L2Nu_13TeV_powheg_pythia8.root\",\"zz2l2n4\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root\",\"zz2l2q4\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/ZZTo4L_13TeV_powheg_pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}ZZTo4L_13TeV_powheg_pythia8.root\",\"zz4l4\",$1\)

root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTTo2L2Nu_13TeV-powheg.root\",\"tt2l\",$1\)
root -l -q -b MitAnalysisRunII/skimming/80x/makeOneSkimSample.C+\(\"$INPUTDIR/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root\",\"$OUTPUTDIR/${PREFIX}TTToSemiLeptonic_13TeV-powheg.root\",\"ttqql\",$1\)
