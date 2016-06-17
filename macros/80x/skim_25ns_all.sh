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

./MitAnalysisRunII/macros/80x/skim_25ns_data.sh $PREFIX;
./MitAnalysisRunII/macros/80x/skim_25ns_mc.sh $PREFIX;
./MitAnalysisRunII/macros/80x/skim_25ns_dm.sh $PREFIX;
