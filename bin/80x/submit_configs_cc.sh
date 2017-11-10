#!/bin/sh

#export MIT_CATALOG=~cmsprod/catalog;
export MIT_CATALOG=~ceballos/catalog;
export MIT_PROD_LOGS=$HOME/cms/logs;
export MIT_PROD_HIST=$HOME/cms/hist;
#export MIT_DATA=$CMSSW_BASE/src/MitPhysics/data;
export MIT_DATA=/cvmfs/cvmfs.cmsaf.mit.edu/hidsk0001/cmsprod/cms/MitPhysics/data
cd $CMSSW_BASE/src/MitAna/bin;

./runOnDatasets.py --cfg=$CMSSW_BASE/src/VVScattering/configs/80x/cc_all.txt \
--condor-template=$CMSSW_BASE/src/NeroProducer/Bambu/config/condor_template.jdl \
--analysis=$CMSSW_BASE/src/VVScattering/python/80x/bambuToNero.py;

#--update --name=cc_all \
# --overwrite
cd -;
