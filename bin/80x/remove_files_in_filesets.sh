#!/bin/sh

grep Aborted  ~/cms/logs/ww_all/t2mit/filefi/044/*/*err|wc;
grep Aborted  ~/cms/logs/ww_all/t2mit/filefi/044/*/*err|awk '{split($1,a,":");print"grep root "a[1]"|grep -e Opening|tail -1"}' > yyy0;chmod a+x yyy0;wc yyy0;
./yyy0|awk '{split($6,a,"043/");split(a[2],b,"/");print b[1],b[2]}' > yyy1;
awk '{print "grep -vwE "$2" ~ceballos/cms/condor/ww_all/t2mit/filefi/044/"$1"/run.py > yyy3;mv yyy3 ~ceballos/cms/condor/ww_all/t2mit/filefi/044/"$1"/run.py"}' yyy1 > yyy2;chmod a+x yyy2;
wc yyy2;
./yyy2;
#rm -f yyy?;
