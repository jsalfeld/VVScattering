#!/bin/sh

ls /eos/cms/store/caf/user/ceballos/Nero/output_80x/|grep met_ > metfiles.txt;
ls /eos/cms/store/caf/user/ceballos/Nero/output_80x/|grep -e DarkMatter -e ADDMonoZ -e Unpart -e Inv > dmfiles.txt;
ls /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/|grep met_ > datafiles.txt;
ls /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/|grep photon_Run20 > dataphotonfiles.txt;
ls /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/|grep met_ > topfiles.txt;

awk '{printf("gfal-copy --timeout 200000 --force /eos/cms/store/caf/user/ceballos/Nero/output_80x/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/%s\n",$1,$1);}' metfiles.txt > copy_metfiles.sh;
awk '{printf("gfal-copy --timeout 200000 --force /eos/cms/store/caf/user/ceballos/Nero/output_80x/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/dm/%s\n",$1,$1);}' dmfiles.txt > copy_dmfiles.sh;
awk '{printf("gfal-copy --timeout 200000 --force /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/data/%s\n",$1,$1);}' datafiles.txt > copy_datafiles.sh;
awk '{printf("gfal-copy --timeout 200000 --force /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/data/%s\n",$1,$1);}' dataphotonfiles.txt > copy_dataphotonfiles.sh;
awk '{printf("gfal-copy --timeout 200000 --force /eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/%s\n",$1,$1);}' topfiles.txt > copy_topfiles.sh;

chmod a+x *.sh;
