#!/bin/sh

ls ~/cms/hist/??_*/*/*/*/*|grep AOD|grep -v root|awk '{split($1,a,":");print"hadd -f "a[1]".root "a[1]"/*root"}' > kkk0;
chmod a+x kkk0;./kkk0;rm -f kkk0;

ls ~/cms/hist/??_*/*/*/*/*root|awk '{print"mv "$1" /scratch5/ceballos/ntuples_noweights/"}' > kkk1;
chmod a+x kkk1;./kkk1;rm -f kkk1;

export INPUTDIR="/scratch5/ceballos/ntuples_noweights";
export GOODRUNDIR="/scratch5/ceballos/ntuples_goodrun";
export OUTPUTDIR="/scratch5/ceballos/ntuples_weights";
export ANADIR="/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src";

ls $INPUTDIR|grep -v AODSIM|grep -v 25ns|grep AOD.root|awk '{printf("root -l -q -b %s/VVScattering/macros/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/VVScattering/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' >> kkk2;
chmod a+x kkk2;./kkk2;rm -f kkk2;

#hadd -f $INPUTDIR/data_AOD_Run2015B1_25ns.root $GOODRUNDIR/MuonEG*Run2015B-PromptReco-v1+AOD.root $GOODRUNDIR/DoubleM*Run2015B-PromptReco-v1+AOD.root  $GOODRUNDIR/SingleM*Run2015B-PromptReco-v1+AOD.root $GOODRUNDIR/DoubleE*Run2015B-PromptReco-v1+AOD.root $GOODRUNDIR/SingleE*Run2015B-PromptReco-v1+AOD.root;
hadd -f $INPUTDIR/data_AOD_Run2015C1_25ns.root $GOODRUNDIR/MuonEG*Run2015C-PromptReco-v1+AOD.root $GOODRUNDIR/DoubleM*Run2015C-PromptReco-v1+AOD.root  $GOODRUNDIR/SingleM*Run2015C-PromptReco-v1+AOD.root $GOODRUNDIR/DoubleE*Run2015C-PromptReco-v1+AOD.root $GOODRUNDIR/SingleE*Run2015C-PromptReco-v1+AOD.root;
hadd -f $INPUTDIR/data_AOD_Run2015D3_25ns.root $GOODRUNDIR/MuonEG*Run2015D-PromptReco-v3+AOD.root $GOODRUNDIR/DoubleM*Run2015D-PromptReco-v3+AOD.root  $GOODRUNDIR/SingleM*Run2015D-PromptReco-v3+AOD.root $GOODRUNDIR/DoubleE*Run2015D-PromptReco-v3+AOD.root $GOODRUNDIR/SingleE*Run2015D-PromptReco-v3+AOD.root;
hadd -f $INPUTDIR/data_AOD_Run2015D4_25ns.root $GOODRUNDIR/MuonEG*Run2015D-PromptReco-v4+AOD.root $GOODRUNDIR/DoubleM*Run2015D-PromptReco-v4+AOD.root  $GOODRUNDIR/SingleM*Run2015D-PromptReco-v4+AOD.root $GOODRUNDIR/DoubleE*Run2015D-PromptReco-v4+AOD.root $GOODRUNDIR/SingleE*Run2015D-PromptReco-v4+AOD.root;
hadd -f $INPUTDIR/data_AOD_spho_25ns.root      $GOODRUNDIR/SingleP*.root;

ls $INPUTDIR|grep -v AODSIM|grep AOD|grep 25ns|awk '{printf("root -l -q -b %s/VVScattering/macros/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/VVScattering/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' >> kkk3;
chmod a+x kkk3;./kkk3;rm -f kkk3;

#ls $GOODRUNDIR|grep AOD|grep 25ns|awk '{printf("root -l -q -b %s/VVScattering/macros/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["GOODRUNDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk4;
#chmod a+x kkk4;./kkk4;rm -f kkk4;

#ls $INPUTDIR|grep AODSIM|awk '{printf("root -l -q -b %s/VVScattering/macros/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk5;
#chmod a+x kkk5;
#echo "kkk5 must be edited"
