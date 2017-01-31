#!/bin/sh

voms-proxy-init --valid 168:00 -voms cms;

for datasetName in `cat list_of_datasets.txt`; do
  export theFiles=list_of_files_$RANDOM.txt;
  das_client --query "file dataset=${datasetName}" | grep "store/mc" > ${theFiles}.txt;
  for fileName in `cat ${theFiles}.txt`; do
    cmsRun ana.py inputFiles=${fileName} maxEvents=-1;
  done
  rm -f ${theFiles}.txt;
done
