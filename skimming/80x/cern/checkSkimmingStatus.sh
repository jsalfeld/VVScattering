
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
