#!/bin/sh

export CATALOGA=~cmsprod/catalog/t2mit/filefi/042;
export CATALOGB=~ceballos/catalog/t2mit/filefi/042;

rm -f list_copy_datasets;

cat > list_copy_datasets <<EOF
MuonEG+Run2015D-PromptReco-v4+AOD								     
DoubleMuon+Run2015D-PromptReco-v4+AOD								     
SingleMuon+Run2015D-PromptReco-v4+AOD								     
DoubleEG+Run2015D-PromptReco-v4+AOD								     
SingleElectron+Run2015D-PromptReco-v4+AOD
SinglePhoton+Run2015D-PromptReco-v4+AOD
EOF

for datasetName in `cat list_copy_datasets`; do
echo ${CATALOGB}/${datasetName};
rm -rf ${CATALOGB}/${datasetName};
cp -r ${CATALOGA}/${datasetName} ${CATALOGB}/.;
done

rm -f list_copy_datasets;
