#!/bin/bash

#outputDir=/store/user/jsalfeld/vdmAnalysis/beamImaging/statisticsEffect/
#outputDir=/afs/cern.ch/work/j/jsalfeld/MIT/vdmAnalysis/CMSSW_7_4_2/src/vdm_2582015_Analysis/outFiles2//
#outputDir=/store/user/jsalfeld/vdmScan_2482025/BI_histosNEW7_ext/
#outputDir=/store/user/jsalfeld/vdmScan_2482025/BI_histosNEW7_ext_withLSC_z/
#outputDir=/store/user/jsalfeld/vdmScan_2482025/ReRecoVdMPixVtxhistos5/
outputDir=/store/user/jsalfeld/
workDirLocal=`pwd`
echo $workDirLocal

script=jobs.sh

for j in {0..26}
	do
	       
		bsub -o out.%J -q 8nh  $script $j $workDirLocal $outputDir
done