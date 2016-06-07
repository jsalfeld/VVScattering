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
else
  echo "LEPC"
  root -l -q -b MitAnalysisRunII/macros/76x/makeSkimSample.C+\(\"/scratch/ceballos/ntuples_goodrun_80x/data.root\",\"/scratch/ceballos/ntuples_weightsDA_80x/${PREFIX}data.root\",\"data\",$1\)
fi
