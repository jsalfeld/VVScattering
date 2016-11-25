#!/bin/sh

### 046 ###
export CATALOGA=~cmsprod/catalog/t2mit/filefi/046;
export CATALOGB=~ceballos/catalog/t2mit/filefi/046;

rm -f list_copy_datasets;

cat > list_copy_datasets <<EOF
DoubleEG+Run2016B-23Sep2016-v2+AOD
DoubleEG+Run2016B-23Sep2016-v3+AOD
DoubleEG+Run2016C-23Sep2016-v1+AOD
DoubleEG+Run2016D-23Sep2016-v1+AOD
DoubleEG+Run2016E-23Sep2016-v1+AOD
DoubleEG+Run2016F-23Sep2016-v1+AOD
DoubleEG+Run2016G-23Sep2016-v1+AOD
DoubleEG+Run2016H-PromptReco-v1+AOD
DoubleEG+Run2016H-PromptReco-v2+AOD
DoubleEG+Run2016H-PromptReco-v3+AOD
DoubleMuon+Run2016B-23Sep2016-v1+AOD
DoubleMuon+Run2016B-23Sep2016-v3+AOD
DoubleMuon+Run2016C-23Sep2016-v1+AOD
DoubleMuon+Run2016D-23Sep2016-v1+AOD
DoubleMuon+Run2016E-23Sep2016-v1+AOD
DoubleMuon+Run2016F-23Sep2016-v1+AOD
DoubleMuon+Run2016G-23Sep2016-v1+AOD
DoubleMuon+Run2016H-PromptReco-v1+AOD
DoubleMuon+Run2016H-PromptReco-v2+AOD
DoubleMuon+Run2016H-PromptReco-v3+AOD
MuonEG+Run2016B-23Sep2016-v2+AOD
MuonEG+Run2016B-23Sep2016-v3+AOD
MuonEG+Run2016C-23Sep2016-v1+AOD
MuonEG+Run2016D-23Sep2016-v1+AOD
MuonEG+Run2016E-23Sep2016-v1+AOD
MuonEG+Run2016F-23Sep2016-v1+AOD
MuonEG+Run2016G-23Sep2016-v1+AOD
MuonEG+Run2016H-PromptReco-v1+AOD
MuonEG+Run2016H-PromptReco-v2+AOD
MuonEG+Run2016H-PromptReco-v3+AOD
SingleElectron+Run2016B-23Sep2016-v2+AOD
SingleElectron+Run2016B-23Sep2016-v3+AOD
SingleElectron+Run2016C-23Sep2016-v1+AOD
SingleElectron+Run2016D-23Sep2016-v1+AOD
SingleElectron+Run2016E-23Sep2016-v1+AOD
SingleElectron+Run2016F-23Sep2016-v1+AOD
SingleElectron+Run2016G-23Sep2016-v1+AOD
SingleElectron+Run2016H-PromptReco-v1+AOD
SingleElectron+Run2016H-PromptReco-v2+AOD
SingleElectron+Run2016H-PromptReco-v3+AOD
SingleMuon+Run2016B-23Sep2016-v1+AOD
SingleMuon+Run2016B-23Sep2016-v3+AOD
SingleMuon+Run2016C-23Sep2016-v1+AOD
SingleMuon+Run2016D-23Sep2016-v1+AOD
SingleMuon+Run2016E-23Sep2016-v1+AOD
SingleMuon+Run2016F-23Sep2016-v1+AOD
SingleMuon+Run2016G-23Sep2016-v1+AOD
SingleMuon+Run2016H-PromptReco-v1+AOD
SingleMuon+Run2016H-PromptReco-v2+AOD
SingleMuon+Run2016H-PromptReco-v3+AOD
SinglePhoton+Run2016B-23Sep2016-v1+AOD
SinglePhoton+Run2016B-23Sep2016-v3+AOD
SinglePhoton+Run2016C-23Sep2016-v1+AOD
SinglePhoton+Run2016D-23Sep2016-v1+AOD
SinglePhoton+Run2016E-23Sep2016-v1+AOD
SinglePhoton+Run2016F-23Sep2016-v1+AOD
SinglePhoton+Run2016G-23Sep2016-v1+AOD
SinglePhoton+Run2016H-PromptReco-v1+AOD
SinglePhoton+Run2016H-PromptReco-v2+AOD
SinglePhoton+Run2016H-PromptReco-v3+AOD
EOF

for datasetName in `cat list_copy_datasets`; do
echo ${CATALOGB}/${datasetName};
rm -rf ${CATALOGB}/${datasetName};
cp -r ${CATALOGA}/${datasetName} ${CATALOGB}/.;
done

rm -f list_copy_datasets;
