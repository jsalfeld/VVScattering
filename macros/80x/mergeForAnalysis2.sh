#!/bin/bash
outputDir=/store/user/jsalfeld/
workDirLocal=`pwd`
echo $workDirLocal

for j in {0..59}
do
for k in {0..5}
do
targetfile=histowz_nice_${k}_${j}_test.root
sourcefiles=`ls *_histowz_nice_${k}_${j}.root`
echo `ls *_histowz_nice_${k}_${j}.root`
hadd -f $targetfile *_histowz_nice_${k}_${j}.root
done
done