#!/bin/sh

grep -e Successfully -e "After filter" log|awk '{split($7,a,"root://eoscms.cern.ch//eos/cms/store/mc/RunIISummer16MiniAODv2/");split(a[2],b,"/"); if($4=="Successfully")printf("%s ",b[1]);if($1=="After")printf("%20.15f %20.15f\n",$7,$9);}' > log_xs;

for datasetName in `cat list_of_datasets.txt`; do
   grep ${datasetName} list_of_datasets.txt|awk '{split($1,a,"/");printf("grep %s log_xs\n",a[2]);}' > kkkxs_1;
   chmod a+x kkkxs_1; ./kkkxs_1| awk '{print$2,$3,$1}'|awk -f compute_average.awk;
done
rm -f kkkxs_1;
