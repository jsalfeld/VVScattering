#!/bin/sh

ls ~/cms/hist/ww_*/*/*/*/*|grep AOD|grep -v root|awk '{split($1,a,":");print"hadd -f "a[1]".root "a[1]"/*root"}' > kkk0;
chmod a+x kkk0;./kkk0;rm -f kkk0;

ls ~/cms/hist/ww_*/*/*/*/*root|awk '{print"mv "$1" /scratch5/ceballos/ntuples_noweights/"}' > kkk1;
chmod a+x kkk1;./kkk1;rm -f kkk1;

cd /scratch5/ceballos/ntuples_noweights/;
mv *AOD.root backup/;
hadd -f data_AOD_50ns.root backup/*.root;
hadd -f data_AOD_25ns.root backup/*.root;

cd -;
