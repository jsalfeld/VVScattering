#/!/bin/sh

cat > list_bad_files <<EOF
5408510F-A43E-E611-ADD4-02163E011E69.root MuonEG+Run2016C-PromptReco-v2+AOD filefi/044
1078906E-473F-E611-9C05-02163E0144A4.root MuonEG+Run2016C-PromptReco-v2+AOD filefi/044
3CEAF931-463E-E611-BBFC-02163E011D88.root SingleMuon+Run2016C-PromptReco-v2+AOD filefi/044
A6733DC7-1241-E611-9FEC-02163E01453D.root SingleMuon+Run2016C-PromptReco-v2+AOD filefi/044
2EA177B1-CD3D-E611-BA0F-02163E014137.root SingleMuon+Run2016C-PromptReco-v2+AOD filefi/044
A6733DC7-1241-E611-9FEC-02163E01453D.root SingleMuon+Run2016C-PromptReco-v2+AOD filefi/044
529C2AA6-DF3C-E611-8B09-02163E011A4F.root DoubleEG+Run2016C-PromptReco-v2+AOD filefi/044
6837746D-213E-E611-AA1D-02163E013873.root DoubleEG+Run2016C-PromptReco-v2+AOD filefi/044
FACA8BB7-9F3D-E611-8DAC-02163E011D08.root DoubleEG+Run2016C-PromptReco-v2+AOD filefi/044
5E656321-5A3E-E611-BB92-02163E0133A9.root SingleElectron+Run2016C-PromptReco-v2+AOD filefi/044
541E37F2-C93E-E611-8409-02163E01253C.root SingleElectron+Run2016C-PromptReco-v2+AOD filefi/044
6A4F943F-853C-E611-9A87-02163E0129EA.root SingleElectron+Run2016C-PromptReco-v2+AOD filefi/044

D0D57485-9834-E611-A573-02163E0139A3.root SingleElectron+Run2016B-PromptReco-v2+AOD filefi/044
881BEEB6-0434-E611-8D83-02163E0135DE.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
7ADB852B-B736-E611-BBD1-02163E011A27.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
D2782C96-BD38-E611-AA28-02163E014681.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
B67A8234-CB2D-E611-B3BB-02163E01270F.root DoubleMuon+Run2016B-PromptReco-v2+AOD filefi/044
98C4E8B0-5233-E611-8212-02163E014544.root SingleMuon+Run2016B-PromptReco-v2+AOD filefi/044
F4C8F439-2F2A-E611-98CF-02163E014621.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
8A7131E1-2B29-E611-8326-02163E0146FD.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
CA2DA46C-2B2A-E611-B098-02163E0137BA.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
0A6D6415-1B34-E611-A379-02163E014418.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
4A5BC2F1-E62C-E611-A7E1-02163E012114.root SingleMuon+Run2016B-PromptReco-v2+AOD filefi/044
76EA4789-CB1F-E611-8FF9-02163E013526.root SingleMuon+Run2016B-PromptReco-v2+AOD filefi/044
F8BD7C56-4337-E611-98A1-02163E0137C4.root DoubleEG+Run2016B-PromptReco-v2+AOD filefi/044
EOF

grep Run2016 list_bad_files > list_bad_files_aux; mv list_bad_files_aux list_bad_files;

awk '{print"grep -vwE "$1" ~/catalog/t2mit/"$3"/"$2"/Files       > ooo;wc ooo ~/catalog/t2mit/"$3"/"$2"/Files;      mv ooo ~/catalog/t2mit/"$3"/"$2"/Files;      "}' list_bad_files  > zzz
awk '{print"grep -vwE "$1" ~/catalog/t2mit/"$3"/"$2"/RawFiles.00 > ooo;wc ooo ~/catalog/t2mit/"$3"/"$2"/RawFiles.00;mv ooo ~/catalog/t2mit/"$3"/"$2"/RawFiles.00;"}' list_bad_files >> zzz
chmod a+x zzz;./zzz;rm -f zzz;

grep Run2016B list_bad_files > list_bad_files_b;
grep Run2016C list_bad_files > list_bad_files_c;
grep Run2016D list_bad_files > list_bad_files_d;

awk '{print"grep -vwE "$1" ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/Files > ooo;wc ooo ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/Files;mv ooo ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/Files;"}' list_bad_files_b  > zzz
awk '{print"grep -vwE "$1" ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/Files > ooo;wc ooo ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/Files;mv ooo ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/Files;"}' list_bad_files_c >> zzz
awk '{print"grep -vwE "$1" ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/Files > ooo;wc ooo ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/Files;mv ooo ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/Files;"}' list_bad_files_d >> zzz
chmod a+x zzz;./zzz;rm -f zzz;

awk '{print"grep -vwE "$1" ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/run.py > ooo;wc ooo ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/run.py;mv ooo ~ceballos/cms/condor/aa_all/t2mit/"$3"/"$2"/run.py;"}' list_bad_files_b  > zzz
awk '{print"grep -vwE "$1" ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/run.py > ooo;wc ooo ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/run.py;mv ooo ~ceballos/cms/condor/bb_all/t2mit/"$3"/"$2"/run.py;"}' list_bad_files_c >> zzz
awk '{print"grep -vwE "$1" ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/run.py > ooo;wc ooo ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/run.py;mv ooo ~ceballos/cms/condor/cc_all/t2mit/"$3"/"$2"/run.py;"}' list_bad_files_d >> zzz
chmod a+x zzz;./zzz;rm -f zzz;
rm -f list_bad_files*;
