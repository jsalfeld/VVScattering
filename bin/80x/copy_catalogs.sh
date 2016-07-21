#!/bin/sh

### 044 ###
export CATALOGA=~cmsprod/catalog/t2mit/filefi/044;
export CATALOGB=~ceballos/catalog/t2mit/filefi/044;

rm -f list_copy_datasets;

cat > list_copy_datasets <<EOF
MuonEG+Run2016B-PromptReco-v2+AOD 
DoubleMuon+Run2016B-PromptReco-v2+AOD 
DoubleEG+Run2016B-PromptReco-v2+AOD 
SingleMuon+Run2016B-PromptReco-v2+AOD 
SingleElectron+Run2016B-PromptReco-v2+AOD 
SinglePhoton+Run2016B-PromptReco-v2+AOD 
MuonEG+Run2016C-PromptReco-v2+AOD 
DoubleMuon+Run2016C-PromptReco-v2+AOD 
DoubleEG+Run2016C-PromptReco-v2+AOD 
SingleMuon+Run2016C-PromptReco-v2+AOD 
SingleElectron+Run2016C-PromptReco-v2+AOD 
SinglePhoton+Run2016C-PromptReco-v2+AOD 
EOF

for datasetName in `cat list_copy_datasets`; do
echo ${CATALOGB}/${datasetName};
rm -rf ${CATALOGB}/${datasetName};
cp -r ${CATALOGA}/${datasetName} ${CATALOGB}/.;
done

rm -f list_copy_datasets;

### 045 ###
export CATALOGA=~cmsprod/catalog/t2mit/filefi/045;
export CATALOGB=~ceballos/catalog/t2mit/filefi/045;

rm -f list_copy_datasets;

cat > list_copy_datasets <<EOF
MuonEG+Run2016D-PromptReco-v2+AOD 
DoubleMuon+Run2016D-PromptReco-v2+AOD 
DoubleEG+Run2016D-PromptReco-v2+AOD 
SingleMuon+Run2016D-PromptReco-v2+AOD 
SingleElectron+Run2016D-PromptReco-v2+AOD 
SinglePhoton+Run2016D-PromptReco-v2+AOD 
EOF

for datasetName in `cat list_copy_datasets`; do
echo ${CATALOGB}/${datasetName};
rm -rf ${CATALOGB}/${datasetName};
cp -r ${CATALOGA}/${datasetName} ${CATALOGB}/.;
done

rm -f list_copy_datasets;
