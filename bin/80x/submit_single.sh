#!/bin/sh

export MIT_CATALOG=~cmsprod/catalog;
export MIT_PROD_LOGS=$HOME/cms/logs;
export MIT_PROD_HIST=$HOME/cms/hist;
#export MIT_DATA=$CMSSW_BASE/src/MitPhysics/data;
export MIT_DATA=/cvmfs/cvmfs.cmsaf.mit.edu/hidsk0001/cmsprod/cms/MitPhysics/data

$CMSSW_BASE/src/MitAna/bin/analysis.py \
--nentries=$1 --book=t2mit/filefi/046 --dataset=$2 $3 --custom bx=25ns \
--fileset=0013 --output=test.root $CMSSW_BASE/src/VVScattering/python/80x/bambuToNero.py;

#$CMSSW_BASE/src/MitAna/bin/analysis.py \
#--nentries=$1 --book=t2mit/fullsm/044 --dataset=$2 $3 --custom bx=25ns \
#--fileset=0000 --output=test.root $CMSSW_BASE/src/VVScattering/python/80x/bambuToNero.py;
