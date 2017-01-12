#!/bin/sh

if [ $# == 1 ] && [ $1 == 1 ]; then

export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/setup80x_ichep/Data/Nero/v2.0/;
export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x;

for PD in DoubleEG DoubleMuon MuonEG SingleMuon SingleElectron SinglePhoton;
#for PD in SingleMuon;
do
  for ERA in Run2016B Run2016C Run2016D Run2016E Run2016F Run2016G Run2016H;
  #for ERA in Run2016H;
  do
    echo ${PD} ${ERA}
    export thePD=${PD}
    export theERA=${ERA}
    ls -l $INPUTDIR/${PD}/*${ERA}*/*/*/*.root|awk '{sp=ENVIRON["INPUTDIR"]"/"ENVIRON["thePD"];split($9,a,sp);print a[2];}' |sort -u > comp_inp_${PD}_${ERA}.txt;
    ls -l  $SKIMDIR/${PD}/*${ERA}*/*/*/*.root|awk '{sp=ENVIRON["SKIMDIR"]"/"ENVIRON["thePD"]; split($9,a,sp);print a[2];}' |sort -u > comp_out_${PD}_${ERA}.txt;
    wc comp_inp_${PD}_${ERA}.txt; wc comp_out_${PD}_${ERA}.txt;
    diff comp_inp_${PD}_${ERA}.txt comp_out_${PD}_${ERA}.txt|grep "<"|awk '{split($2,a,"NeroNtuples_");printf("./MitAnalysisRunII/skimming/80x/cern/skim_batch.sh %s/%s %s/%s %s NeroNtuples_%s data\n",ENVIRON["INPUTDIR"],ENVIRON["thePD"],ENVIRON["SKIMDIR"],ENVIRON["thePD"],a[1],a[2])}' > diff_${PD}_${ERA}.sh;
    rm -f comp_inp_${PD}_${ERA}.txt comp_out_${PD}_${ERA}.txt;
  done
done

elif [ $# == 2 ] && [ $1 == 2 ]; then

if [ $2 == 1 ]; then
export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/user/ceballos/Nero/v2.1;
export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/caf/user/ceballos/Nero/skim_80x;
else
export INPUTDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/v2.1;
export SKIMDIR=/afs/cern.ch/user/c/ceballos/eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x;
fi

export theRND=file_testskim_$RANDOM;

find $INPUTDIR/* -name '000?'| awk '{split($1,a,ENVIRON["INPUTDIR"]);print a[2]}' > ${theRND}.txt
awk '{printf("ls -l %s%s/*.root\n",ENVIRON["INPUTDIR"],$1)}' ${theRND}.txt > ${theRND}_sh1.sh; chmod a+x ${theRND}_sh1.sh;
awk '{printf("ls -l %s%s/*.root\n",ENVIRON["SKIMDIR"],$1)}' ${theRND}.txt > ${theRND}_sh2.sh;  chmod a+x ${theRND}_sh2.sh;

./${theRND}_sh1.sh|awk '{sp=ENVIRON["INPUTDIR"];split($9,a,sp);print a[2];}' |sort -u > ${theRND}_comp_inp.txt;
./${theRND}_sh2.sh|awk '{sp=ENVIRON["SKIMDIR"]; split($9,a,sp);print a[2];}' |sort -u > ${theRND}_comp_out.txt;
wc ${theRND}_comp_inp.txt ${theRND}_comp_out.txt;
diff ${theRND}_comp_inp.txt ${theRND}_comp_out.txt|grep "<"|awk '{split($2,a,"NeroNtuples_");printf("./MitAnalysisRunII/skimming/80x/cern/skim_batch.sh %s %s %s NeroNtuples_%s dm\n",ENVIRON["INPUTDIR"],ENVIRON["SKIMDIR"],a[1],a[2])}' > diff_${theRND}.sh;

sed -i 's|/afs/cern.ch/user/c/ceballos/eos/cms|cms|' diff_${theRND}.sh;
sed -i 's|/afs/cern.ch/user/c/ceballos/eos/cms|cms|' diff_${theRND}.sh;
chmod a+x diff_${theRND}.sh;
rm -f ${theRND}*;

else
  echo "Wrong option";
fi
