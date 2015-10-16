#!/bin/sh

export INPUTDIR="/scratch5/ceballos/ntuples_noweights";
export GOODRUNDIR="/scratch5/ceballos/ntuples_goodrun";
export OUTPUTDIR="/scratch5/ceballos/ntuples_weights";
export ANADIR="/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src";
export FILTER=0;

ls $INPUTDIR|grep -v AODSIM|grep AOD|grep 50ns|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/MitAnalysisRunII/json/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_v2.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' > kkk0;
ls $INPUTDIR|grep -v AODSIM|grep AOD|grep 25ns|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/MitAnalysisRunII/json/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' >> kkk0;
chmod a+x kkk0;./kkk0;rm -f kkk0;

#ls $GOODRUNDIR|grep AOD|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["GOODRUNDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk1;
#chmod a+x kkk1;./kkk1;rm -f kkk1;

#ls $INPUTDIR|grep AODSIM|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk2;
#chmod a+x kkk2;
#echo "kkk2 must be edited"
