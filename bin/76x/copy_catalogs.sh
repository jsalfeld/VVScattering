#!/bin/sh

export CATALOGA=~cmsprod/catalog/t2mit/filefi/043;
export CATALOGB=~ceballos/catalog/t2mit/filefi/043;

rm -f list_copy_datasets;

cat > list_copy_datasets <<EOF
MuonEG+Run2015C_25ns-16Dec2015-v1+AOD
MuonEG+Run2015D-16Dec2015-v1+AOD
DoubleMuon+Run2015C_25ns-16Dec2015-v1+AOD
DoubleMuon+Run2015D-16Dec2015-v1+AOD
DoubleEG+Run2015C_25ns-16Dec2015-v1+AOD
DoubleEG+Run2015D-16Dec2015-v2+AOD
SingleMuon+Run2015C_25ns-16Dec2015-v1+AOD
SingleMuon+Run2015D-16Dec2015-v1+AOD
SingleElectron+Run2015C_25ns-16Dec2015-v1+AOD
SingleElectron+Run2015D-16Dec2015-v1+AOD
SinglePhoton+Run2015C_25ns-16Dec2015-v1+AOD
SinglePhoton+Run2015D-16Dec2015-v1+AOD
EOF

for datasetName in `cat list_copy_datasets`; do
echo ${CATALOGB}/${datasetName};
rm -rf ${CATALOGB}/${datasetName};
cp -r ${CATALOGA}/${datasetName} ${CATALOGB}/.;
done

rm -f list_copy_datasets;
