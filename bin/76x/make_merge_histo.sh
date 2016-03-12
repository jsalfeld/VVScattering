#!/bin/sh

ls ~/cms/hist/??_*/*/*/*/*|grep AOD|grep -v root|awk '{split($1,a,":");print"hadd -f "a[1]".root "a[1]"/*root"}' > kkk0;
chmod a+x kkk0;./kkk0;rm -f kkk0;

ls ~/cms/hist/??_*/*/*/*/*root|awk '{print"mv "$1" /scratch5/ceballos/ntuples_noweights_76x/"}' > kkk1;
chmod a+x kkk1;./kkk1;rm -f kkk1;

export INPUTDIR="/scratch5/ceballos/ntuples_noweights_76x";
export GOODRUNDIR="/scratch/ceballos/ntuples_goodrun_76x";
export OUTPUTDADIR="/scratch/ceballos/ntuples_weightsDA_76x";
export OUTPUTMCDIR="/scratch5/ceballos/ntuples_weightsMC_76x";
export ANADIR="/home/ceballos/cms/cmssw/043/CMSSW_7_6_3/src";

ls $INPUTDIR|grep -v AODSIM|grep AOD.root|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/76x/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/MitAnalysisRunII/json/76x/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' > kkk2;
chmod a+x kkk2;./kkk2;rm -f kkk2;

#ls $GOODRUNDIR|grep AOD|grep 25ns|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/76x/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["GOODRUNDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk4;
#chmod a+x kkk4;./kkk4;rm -f kkk4;

#ls $INPUTDIR|grep AODSIM|awk '{printf("root -l -q -b %s/MitAnalysisRunII/macros/76x/makeSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"data\\\",%d\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["OUTPUTDIR"],$1,ENVIRON["ANADIR"],ENVIRON["FILTER"])}' > kkk5;
#chmod a+x kkk5;
#echo "kkk5 must be edited"
