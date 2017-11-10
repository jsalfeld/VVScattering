#!/bin/sh

export PREFIX="";
if [ $1 == 1 ]
then
  export PREFIX="qcd_";
elif [ $1 == 2 ]
then
  export PREFIX="met_";
elif [ $1 == 3 ]
then
  export PREFIX="wln_";
elif [ $1 == 4 ]
then
  export PREFIX="pho_";
fi

if [ $1 == 4 ]
then
  echo "PHO"
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SinglePhoton+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_weightsDA_76x/SinglePhoton+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SinglePhoton+Run2015D-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_weightsDA_76x/SinglePhoton+Run2015D-16Dec2015-v1+AOD.root\",\"data\",$1\)

  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/DYJetsToNuNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}DYJetsToNuNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"dynunu\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"gjets40_100\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"gjets100_200\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"gjets200_400\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"gjets400_600\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch5/ceballos/ntuples_noweights_76x/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"/scratch5/ceballos/ntuples_weightsMC_76x/${PREFIX}GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root\",\"gjets600\",$1\)
else
  echo "LEPC"
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/MuonEG+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}MuonEG+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/DoubleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SingleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleMuon+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/DoubleEG+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleEG+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SingleElectron+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleElectron+Run2015C_25ns-16Dec2015-v1+AOD.root\",\"data\",$1\)

  hadd -f /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}data_AOD_Run2015C_25ns.root \
          /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}MuonEG*Run2015C*+AOD.root \
	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleM*Run2015C*+AOD.root  \
	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleM*Run2015C*+AOD.root \
	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleE*Run2015C*+AOD.root \
	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleE*Run2015C*+AOD.root;

  root -l -q -b VVScattering/macros/76x/makeGoodRunSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}data_AOD_Run2015C_25ns.root\",\"/scratch/ceballos/ntuples_weightsDA_76x/${PREFIX}data_AOD_Run2015C_25ns.root\",\"VVScattering/json/76x/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt\"\)

  echo "LEPD"
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/MuonEG+Run2015D-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}MuonEG+Run2015D-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/DoubleMuon+Run2015D-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleMuon+Run2015D-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SingleMuon+Run2015D-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleMuon+Run2015D-16Dec2015-v1+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/DoubleEG+Run2015D-16Dec2015-v2+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleEG+Run2015D-16Dec2015-v2+AOD.root\",\"data\",$1\)
  root -l -q -b VVScattering/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/SingleElectron+Run2015D-16Dec2015-v1+AOD.root\",\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleElectron+Run2015D-16Dec2015-v1+AOD.root\",\"data\",$1\)

  hadd -f /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}data_AOD_Run2015D_25ns.root \
  	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}MuonEG*Run2015D*+AOD.root \
  	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleM*Run2015D*+AOD.root  \
  	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleM*Run2015D*+AOD.root \
  	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}DoubleE*Run2015D*+AOD.root \
  	  /scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}SingleE*Run2015D*+AOD.root;

  root -l -q -b VVScattering/macros/76x/makeGoodRunSample.C+\(\"/scratch/ceballos/ntuples_goodrun_76x/temp/${PREFIX}data_AOD_Run2015D_25ns.root\",\"/scratch/ceballos/ntuples_weightsDA_76x/${PREFIX}data_AOD_Run2015D_25ns.root\",\"VVScattering/json/76x/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt\"\)
fi
