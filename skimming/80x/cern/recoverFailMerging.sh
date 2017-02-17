#!/bin/sh

if [ $# == 1 ] && [ $1 == 1 ]; then

grep -i err ppp|wc
grep -i err ppp|grep "Failed to read data"|wc
grep -i err ppp|grep "Input/output error"|wc
grep -i err ppp|grep "Error in <TFile::ReadBuffer>:"|grep -v "Input/output error"|wc
grep -i err ppp|grep -v "Failed to read data"|grep -v "Error in <TFile::ReadBuffer>:"|wc

grep -i err ppp|grep "Input/output error"|awk '{print$8}' > ppp1
grep -i err ppp|grep "Error in <TFile::ReadBuffer>:"|grep -v "Input/output error"|awk '{split($11,a,",");print a[1];}' >> ppp1
grep -i err ppp|grep -v "Failed to read data"|grep -v "Error in <TFile::ReadBuffer>:"|awk '{print$5}' >> ppp1

sort -u ppp1 > ppp2;
wc ppp?;

for PD in DoubleEG DoubleMuon MuonEG SingleMuon SingleElectron SinglePhoton MET;
#for PD in  MuonEG SinglePhoton MET;
do
  for ERA in Run2016B Run2016C Run2016D Run2016E Run2016F Run2016G Run2016H;
  do
    echo ${PD} ${ERA}
    export thePD=${PD}
    export theERA=${ERA}
    grep ${PD}-${ERA} ppp2|awk '{if(NR==1) printf("hadd -O -f /eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/%s-%s_failed.root \\\n",ENVIRON["thePD"],ENVIRON["theERA"]);printf("%s ",$1)}END{printf("\n")}'>  failed_${PD}_${ERA}.sh;
    grep ${PD}_${ERA} ppp2|awk '{if(NR==1) printf("hadd -O -f /eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/%s_%s_failed.root \\\n",ENVIRON["thePD"],ENVIRON["theERA"]);printf("%s ",$1)}END{printf("\n")}'>> failed_${PD}_${ERA}.sh;
    chmod a+x failed_${PD}_${ERA}.sh;
  done
done
rm -f ppp ppp1 ppp2;

elif [ $# == 1 ] && [ $1 == 2 ]; then

#for PD in DoubleEG DoubleMuon MuonEG SingleMuon SingleElectron SinglePhoton MET;
for PD in SingleMuon;
do
  echo "***********${PD}***********"
  find /eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/${PD} -name '000?'| awk '{split($1,a,"/");printf("du -k /eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/%s/%s/*/%s;\n",a[10],a[11],a[13]);}'|sort -u > ppp1;
  find /eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/${PD} -name '000?'| awk '{split($1,a,"/");printf("du -k /eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/%s*%s*.root\n",a[11],a[13]);}'|sort -u > ppp2;
  chmod a+x ppp1;./ppp1;
  chmod a+x ppp2;./ppp2;
  du -k /eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/${PD}*Run2016*failed.root;
  chmod a+x ppp1;./ppp1|awk -f ~/public/bin/sum.awk;rm -f ppp1;
  chmod a+x ppp2;./ppp2|awk -f ~/public/bin/sum.awk;rm -f ppp2;
  du -k /eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/${PD}*Run2016*failed.root|awk -f ~/public/bin/sum.awk;
done

elif [ $# == 2 ] && [ $1 == 1 ]; then

for ERA in Run2016B Run2016C Run2016D Run2016E Run2016F Run2016G Run2016H;
do
  echo "***********${ERA}***********"
  ./failed_$2_${ERA}.sh
done

else
  echo "Wrong option "$#;
fi
