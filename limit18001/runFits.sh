#!/bin/sh

export PARAM="";

if [ $# != 2 ]; then
   echo "TOO FEW PARAMETERS"
   exit
fi

if [ $1 = "exp" ]; then
  PARAM="--expectSignal=1 -t -1";
fi

if [ $2 = "impacts" ]; then
  text2workspace.py datacard.text -m 400;
  combineTool.py -M Impacts -d datacard.text.root -m 400 --robustFit 1 --doInitialFit --allPars $PARAM
  combineTool.py -M Impacts -d datacard.text.root -m 400 --robustFit 1 --doFits --allPars $PARAM
  combineTool.py -M Impacts -d datacard.text.root -m 400 -o impacts_datacard.json --allPars $PARAM
  plotImpacts.py -i impacts_datacard.json -o impacts_datacard;

elif [ $2 = "limits" ]; then
  combine -M AsymptoticLimits datacard.text -m 1500 -n hzhg $PARAM

elif [ $2 = "sig" ]; then
  combine -M Significance datacard.text $PARAM

elif [ $2 = "mu" ]; then
   #combine -M MultiDimFit datacard.text -n hzhg --algo=singles --robustFit=1 $PARAM -m 1500
    combine -M MultiDimFit datacard.text --algo=singles --robustFit=1 $PARAM -m 150
else
   echo "WRONG OPTION"

fi

rm -f combine_logger.out;
rm -f higgsCombine_paramFit_Test_*;
