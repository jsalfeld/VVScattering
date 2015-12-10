#!/bin/sh

condor_q|grep H|awk '{printf("condor_q -long %s|grep logs|grep Err\n",$1)}' > lll0;
chmod a+x lll0;./lll0|awk '{print$3}' > lll1;
sed -i 's|logs|hist|' lll1
sed -i 's|.err||' lll1;
sed -i 's|AOD/0|AOD/*0|' lll1;
sed -i 's|AODSIM/0|AODSIM/*0|' lll1;
sed -i 's|"||' lll1;
sed -i 's|"||' lll1;

awk '{printf("rm %s.root\n",$1)}' lll1 > lll2;
cat lll2;
chmod a+x lll2; ./lll2;

condor_q|grep H|awk '{printf("condor_rm %s\n",$1)}' > lll3;
chmod a+x lll3;./lll3;
rm -f lll0 lll1 lll2 lll3;

cp /home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/bin/submit_configs_ww.sh submit_configs_temp.sh;
sed -i 's|--analysis|--update --name=ww_all --analysis|' submit_configs_temp.sh;
./submit_configs_temp.sh;
rm -f submit_configs_temp.sh;
