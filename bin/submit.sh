#!/bin/sh

export MIT_CATALOG=~cmsprod/catalog;
export MIT_PROD_LOGS=$HOME/cms/logs;
export MIT_PROD_HIST=$HOME/cms/hist;
export MIT_DATA=$CMSSW_BASE/src/MitPhysics/data;
cd $CMSSW_BASE/src/MitAna/bin;
./runOnDatasets.py --name=test_task_$2 --book=t2mit/filefi/041 \
--dataset=$3 $4 \
--analysis=$CMSSW_BASE/src/MitAnalysisRunII/python/$1;
# --overwrite
cd -;
