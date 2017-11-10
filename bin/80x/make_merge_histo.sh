#!/bin/sh

#ls ~/cms/hist/??_*/*/*/*/*|grep AOD|grep -v root|awk '{split($1,a,":");print"hadd -f "a[1]".root "a[1]"/*root"}' > kkk0;
ls ~/cms/hist/??_*/*/*/*/*|grep AOD|grep -v root|awk '{split($1,a,":");split(a[1],b,"/");print"hadd -f /tmp/"b[10]".root "a[1]"/*root;mv /tmp/"b[10]".root "a[1]".root;rm -rf "a[1]";"}' > kkk0;
chmod a+x kkk0;./kkk0;rm -f kkk0;

ls ~/cms/hist/??_*/*/*/*/*root|awk '{print"mv "$1" /scratch5/ceballos/ntuples_noweights_80x/"}' > kkk1;
chmod a+x kkk1;./kkk1;rm -f kkk1;

export INPUTDIR="/scratch5/ceballos/ntuples_noweights_80x";
export GOODRUNDIR="/scratch/ceballos/ntuples_goodrun_80x";
export OUTPUTDADIR="/scratch/ceballos/ntuples_weightsDA_80x";
export OUTPUTMCDIR="/scratch5/ceballos/ntuples_weightsMC_80x";
export ANADIR=$CMSSW_BASE/src;

ls $INPUTDIR|grep -v AODSIM|grep AOD.root|awk '{printf("root -l -q -b %s/VVScattering/macros/80x/makeGoodRunSample.C+\\(\\\"%s/%s\\\",\\\"%s/%s\\\",\\\"%s/VVScattering/json/80x/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["GOODRUNDIR"],$1,ENVIRON["ANADIR"])}' > kkk2;
chmod a+x kkk2;./kkk2;rm -f kkk2;

#ls $INPUTDIR|grep -v AODSIM|grep AOD|awk '{printf("root -l -q -b %s/VVScattering/macros/80x/makeOneSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/${PREFIX}%s\\\",\\\"data\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["GOODRUNDIR"],$1,ENVIRON["OUTPUTDADIR"],$1,ENVIRON["ANADIR"])}' > kkk4;
#chmod a+x kkk4;./kkk4;rm -f kkk4;

#ls $INPUTDIR|grep AODSIM|awk '{printf("root -l -q -b %s/VVScattering/macros/80x/makeOneSkimSample.C+\\(\\\"%s/%s\\\",\\\"%s/${PREFIX}%s\\\",\\\"data\\\"\\)\n",ENVIRON["ANADIR"],ENVIRON["INPUTDIR"],$1,ENVIRON["OUTPUTMCDIR"],$1,ENVIRON["ANADIR"])}' > kkk5;
#chmod a+x kkk5;
#echo "kkk5 must be edited"
