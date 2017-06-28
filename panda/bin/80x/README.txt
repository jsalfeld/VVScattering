--> setup
source /home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/MitAnalysisRunII/panda/bin/80x/setup_leptonic.sh

--> run macro 
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/inputs/skim_lep_tmpl.py

--> Datasets to be included
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/bin/catalogT2Prod.py --outfile ~/public_html/$USER/catalog/test.cfg --catalog ~cmsprod/catalog/t2mit/pandaf/004 --include SingleElectron SingleMuon DoubleEG DoubleMuon MuonEG ZJets_nlo ZJets_lo --exclude EWK

/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/bin/catalogT2Prod.py --outfile ~/public_html/$USER/catalog/test.cfg --catalog ~cmsprod/catalog/t2mit/pandaf/004 \
--include SingleElectron SingleMuon DoubleEG DoubleMuon MuonEG ZZ WZ WW DYJetsToLL_M-50_Tune DYJetsToLL_M-10to50_Tune tZq GluGluH VBFH VBF_H ttHToNonbb VHToNonbb \
TTG TTZ TTW ST_tW TTTo2L2Nu WGstarTo WGToLNuG ZGTo2LG JetsToLL DYJetsToTauTau NNPDF30_13TeV-powheg --exclude ZpWW_med JetsToLL_M-50_HT

--> Processes must be included in 
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaCore/Tools/python/processes/*.py

--> Define how many files to be ran per job
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/bin/buildMergedInputs.sh -t -n 30

--> Submit
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/bin/submit.py

--> Check status
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/bin/checkMissingFiles.py


--> Merge files once all jobs are done
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py ZZ		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py WZ		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py qqWW		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py ggWW		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py WWdps  	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py VVV		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py ttV		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py TT2L		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py TW		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py WGstar 	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py VG		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py H125		    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYNJetsToLL	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYJetsToLL_Pt      
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYJetsToTauTau     
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYJetsToLL_M-10to50
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYJetsToLL_M-50_NLO
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DYJetsToLL_M-50_LO 
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py SingleMuon	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py SingleElectron     
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py MuonEG 	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DoubleEG	    
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py DoubleMuon
/home/ceballos/cms/cmssw/047/CMSSW_8_0_26_patch1/src/PandaAnalysis/T3/merging/merge.py data_overlaps

