#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

enum ERANumber {Run2016B=0,Run2016C,Run2016D,Run2016E,Run2016F,Run2016G,Run2016H};

void checkEvents(int nsel = 0){

  vector<TString> filesPathInput;
  if(nsel == 0){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016B-03Feb2017_ver2/170210_071458/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016B-03Feb2017_ver2/170214_110026/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016B-03Feb2017_ver2/170210_071458/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016B-03Feb2017_ver2/170210_071458/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016C-03Feb2017/170209_121904/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016C-03Feb2017/170214_110149/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016D-03Feb2017/170209_101642/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016D-03Feb2017/170214_110318/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016D-03Feb2017/170209_101642/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016E-03Feb2017/170209_104815/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016E-03Feb2017/170215_072818/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016E-03Feb2017/170209_104815/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016F-03Feb2017/170209_101907/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016F-03Feb2017/170214_110545/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016G-03Feb2017/170209_102159/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016G-03Feb2017/170214_110823/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016G-03Feb2017/170209_102159/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016G-03Feb2017/170209_102159/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016H-03Feb2017_ver2/170213_135256/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016H-03Feb2017_ver2/170215_072936/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016H-03Feb2017_ver2/170213_135256/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016H-03Feb2017_ver2/170213_135256/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleEG/DoubleEG-Run2016H-03Feb2017_ver3/170213_135452/0000/");
  } else if(nsel == 1){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016B-03Feb2017_ver2/170210_072814/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016B-03Feb2017_ver2/170215_074020/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016B-03Feb2017_ver2/170210_072814/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016B-03Feb2017_ver2/170210_072814/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016C-03Feb2017/170210_072930/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016D-03Feb2017/170214_171332/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016D-03Feb2017/170214_171332/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016D-03Feb2017/170214_171332/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016E-03Feb2017/170210_071326/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016E-03Feb2017/170214_111042/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016E-03Feb2017/170210_071326/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016F-03Feb2017/170209_102828/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016G-03Feb2017/170209_103142/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016G-03Feb2017/170214_111248/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016G-03Feb2017/170209_103142/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016G-03Feb2017/170209_103142/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016H-03Feb2017_ver2/170209_105213/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016H-03Feb2017_ver2/170214_111541/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016H-03Feb2017_ver2/170209_105213/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016H-03Feb2017_ver2/170209_105213/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/DoubleMuon/DoubleMuon-Run2016H-03Feb2017_ver3/170209_103504/0000/");
  } else if(nsel == 2){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016B-03Feb2017_ver1-v2/170210_111857/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016B-03Feb2017_ver1-v2/170210_111857/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016B-03Feb2017_ver1-v2/170210_111857/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016C-03Feb2017-v1/170209_221032/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016D-03Feb2017-v1/170209_221244/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016D-03Feb2017-v1/170209_221244/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016E-03Feb2017-v1/170209_221358/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016E-03Feb2017-v1/170209_221358/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016F-03Feb2017-v1/170209_221522/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016G-03Feb2017-v1/170210_112039/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016G-03Feb2017-v1/170210_112039/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016G-03Feb2017-v1/170210_112039/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016H-03Feb2017_ver2-v1/170209_222722/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016H-03Feb2017_ver2-v1/170209_222722/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016H-03Feb2017_ver2-v1/170209_222722/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MET/MET-Run2016H-03Feb2017_ver3-v1/170209_222513/0000/");
  } else if(nsel == 3){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016B_ver2/170209_204840/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016B_ver2/170213_110951/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016B_ver2/170209_204840/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016B_ver2/170209_204840/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016C/170209_204948/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016C/170213_124540/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016D/170209_205054/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016D/170213_113744/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016D/170209_205054/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016E/170210_131800/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016E/170213_114010/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016E/170210_131800/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016F/170209_205158/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016G/170210_131903/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016G/170213_114119/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016G/170210_131903/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016G/170210_131903/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170210_161755/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170213_114423/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170210_161755/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170213_114423/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170210_161755/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170210_161755/0003/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver2/170210_161755/0004/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver3/170209_205338/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/MuonEG/MuonEG_Run2016H_ver3/170213_173654/0000/");
  } else if(nsel == 4){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016B_ver2/170209_230126/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016B_ver2/170209_230126/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016B_ver2/170209_230126/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016C/170209_230248/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016D/170209_230417/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016D/170209_230417/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016E/170209_230535/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016E/170209_230535/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016F/170212_143458/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016G/170209_230712/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016G/170209_230712/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016G/170209_230712/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016H_ver2/170211_100128/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016H_ver2/170211_100128/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016H_ver2/170211_100128/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleElectron/SingleElectron-Run2016H_ver3/170209_231102/0000/");
  } else if(nsel == 5){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016C/170209_225556/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016D/170211_094307/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016D/170211_094307/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016E/170209_225717/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016E/170209_225717/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016F/170211_094423/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016G/170212_143244/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016G/170212_143244/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016G/170212_143244/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016H_ver2/170209_225838/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016H_ver2/170209_225838/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016H_ver2/170209_225838/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SingleMuon/SingleMuon-Run2016H_ver3/170209_225958/0000/");
  } else if(nsel == 6){
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016B-03Feb2017_ver2-v2/170209_220042/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016B-03Feb2017_ver2-v2/170209_220042/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016B-03Feb2017_ver2-v2/170209_220042/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016C-03Feb2017-v1/170212_212932/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016D-03Feb2017-v1/170209_220206/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016D-03Feb2017-v1/170209_220206/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016E-03Feb2017-v1/170209_220319/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016E-03Feb2017-v1/170209_220319/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016F-03Feb2017-v1/170209_220431/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016G-03Feb2017-v1/170212_213049/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016G-03Feb2017-v1/170212_213049/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016G-03Feb2017-v1/170212_213049/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016H-03Feb2017_ver2-v1/170209_220633/0000/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016H-03Feb2017_ver2-v1/170209_220633/0001/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016H-03Feb2017_ver2-v1/170209_220633/0002/");
  filesPathInput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/skim_80x/SinglePhoton/SinglePhoton-Run2016H-03Feb2017_ver3-v1/170209_220749/0000/");
  } else {
    printf("Wrong option\n"); return;
  }

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamevInput;
  for(unsigned int nd=0; nd<filesPathInput.size(); nd++){
    int edges[2] = {1, 1000};
    if     (filesPathInput[nd].Contains("/0001/")) {edges[0] = 1001; edges[1] = 2000;}
    else if(filesPathInput[nd].Contains("/0002/")) {edges[0] = 2001; edges[1] = 3000;}
    else if(filesPathInput[nd].Contains("/0003/")) {edges[0] = 3001; edges[1] = 4000;}
    else if(filesPathInput[nd].Contains("/0004/")) {edges[0] = 4001; edges[1] = 5000;}
    for(int i=edges[0]; i<=edges[1]; i++){
      infilenamevInput.push_back(Form("%sNeroNtuples_%d.root",filesPathInput[nd].Data(),i));
    }
  }

  int sumInput[7][3];
  for(int i=0; i<7; i++){
    for(int j=0; j<3; j++){
      sumInput[i][j] = 0;
    }
  }

  for(UInt_t ifile=0; ifile<infilenamevInput.size(); ifile++) {
    printf("sampleNamesInput(%d): %s\n",ifile,infilenamevInput[ifile].Data());
    TFile *the_input_file = TFile::Open(infilenamevInput[ifile].Data());
    if(the_input_file){
      TTree *the_input_tree  = (TTree*)the_input_file->FindObjectAny("events");
      TTree *the_input_all   = (TTree*)the_input_file->FindObjectAny("all");
      TTree *the_SelBit_tree = (TTree*)the_input_file->FindObjectAny("SelBit_tree");
      if(the_SelBit_tree->GetEntries() != the_input_tree->GetEntries()) printf("BIG PROBLEM!: %lld != %lld\n",the_SelBit_tree->GetEntries(),the_input_tree->GetEntries());
      int theCat = -1;
      if     (infilenamevInput[ifile].Contains("Run2016B")) theCat = Run2016B;
      else if(infilenamevInput[ifile].Contains("Run2016C")) theCat = Run2016C;
      else if(infilenamevInput[ifile].Contains("Run2016D")) theCat = Run2016D;
      else if(infilenamevInput[ifile].Contains("Run2016E")) theCat = Run2016E;
      else if(infilenamevInput[ifile].Contains("Run2016F")) theCat = Run2016F;
      else if(infilenamevInput[ifile].Contains("Run2016G")) theCat = Run2016G;
      else if(infilenamevInput[ifile].Contains("Run2016H")) theCat = Run2016H;
      else {printf("Wrong ERA\n"); return;}
      sumInput[theCat][0] = sumInput[theCat][0] + the_input_all->GetEntries();
      sumInput[theCat][1] = sumInput[theCat][1] + the_input_tree->GetEntries();
      sumInput[theCat][2] = sumInput[theCat][2] + the_SelBit_tree->GetEntries();
      the_input_file->Close();
    }
  }
  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamevOutput;
  if       (nsel == 0){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016B-03Feb2017_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016B-03Feb2017_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016B-03Feb2017_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016B-03Feb2017_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016C-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016D-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016D-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016E-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016E-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016F-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016G-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016G-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016G-03Feb2017_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016H-03Feb2017_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016H-03Feb2017_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016H-03Feb2017_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016H-03Feb2017_ver3_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleEG-Run2016H_failed.root");
  } else if(nsel == 1){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016B-03Feb2017_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016B-03Feb2017_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016B-03Feb2017_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016C-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016D-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016D-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016D-03Feb2017_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016E-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016E-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016F-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016G-03Feb2017_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016G-03Feb2017_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016H-03Feb2017_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016H-03Feb2017_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016H-03Feb2017_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016H-03Feb2017_ver3_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/DoubleMuon-Run2016H_failed.root");
  } else if(nsel == 2){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016B-03Feb2017_ver1-v2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016B-03Feb2017_ver1-v2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016B-03Feb2017_ver1-v2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016C-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016D-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016D-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016E-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016E-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016F-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016G-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016G-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016G-03Feb2017-v1_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016H-03Feb2017_ver2-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016H-03Feb2017_ver2-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016H-03Feb2017_ver2-v1_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016H-03Feb2017_ver3-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MET-Run2016H_failed.root");
  } else if(nsel == 3){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016B_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016B_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016B_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016C_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016D_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016D_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016E_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016E_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016F_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016G_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016G_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016G_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver2_0003.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver2_0004.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/MuonEG_Run2016H_ver3_0000.root");
  } else if(nsel == 4){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016B_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016B_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016B_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016C_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016D_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016D_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016E_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016E_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016F_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016G_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016G_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016G_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016H_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016H_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016H_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016H_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleElectron-Run2016H_ver3_0000.root");
  } else if(nsel == 5){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016C_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016D_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016D_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016E_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016E_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016F_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016G_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016G_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016G_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016H_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016H_ver2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016H_ver2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016H_ver2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SingleMuon-Run2016H_ver3_0000.root");
  } else if(nsel == 6){
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016B-03Feb2017_ver2-v2_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016B-03Feb2017_ver2-v2_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016B-03Feb2017_ver2-v2_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016B_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016C-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016C_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016D-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016D-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016D_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016E-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016E-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016E_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016F-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016F_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016G-03Feb2017-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016G-03Feb2017-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016G-03Feb2017-v1_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016G_failed.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016H-03Feb2017_ver2-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016H-03Feb2017_ver2-v1_0001.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016H-03Feb2017_ver2-v1_0002.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016H-03Feb2017_ver3-v1_0000.root");
  infilenamevOutput.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/SinglePhoton-Run2016H_failed.root");
  }

  int sumOutput[7][3];
  for(int i=0; i<7; i++){
    for(int j=0; j<3; j++){
      sumOutput[i][j] = 0;
    }
  }

  for(UInt_t ifile=0; ifile<infilenamevOutput.size(); ifile++) {
    printf("sampleNamesOutput(%d): %s\n",ifile,infilenamevOutput[ifile].Data());
    TFile *the_input_file = TFile::Open(infilenamevOutput[ifile].Data());
    if(the_input_file){
      TTree *the_input_tree  = (TTree*)the_input_file->FindObjectAny("events");
      TTree *the_input_all   = (TTree*)the_input_file->FindObjectAny("all");
      TTree *the_SelBit_tree = (TTree*)the_input_file->FindObjectAny("SelBit_tree");
      if(the_SelBit_tree->GetEntries() != the_input_tree->GetEntries()) printf("BIG PROBLEM!: %lld != %lld\n",the_SelBit_tree->GetEntries(),the_input_tree->GetEntries());
      int theCat = -1;
      if     (infilenamevOutput[ifile].Contains("Run2016B")) theCat = Run2016B;
      else if(infilenamevOutput[ifile].Contains("Run2016C")) theCat = Run2016C;
      else if(infilenamevOutput[ifile].Contains("Run2016D")) theCat = Run2016D;
      else if(infilenamevOutput[ifile].Contains("Run2016E")) theCat = Run2016E;
      else if(infilenamevOutput[ifile].Contains("Run2016F")) theCat = Run2016F;
      else if(infilenamevOutput[ifile].Contains("Run2016G")) theCat = Run2016G;
      else if(infilenamevOutput[ifile].Contains("Run2016H")) theCat = Run2016H;
      else {printf("Wrong ERA\n"); return;}
      sumOutput[theCat][0] = sumOutput[theCat][0] + the_input_all->GetEntries();
      sumOutput[theCat][1] = sumOutput[theCat][1] + the_input_tree->GetEntries();
      sumOutput[theCat][2] = sumOutput[theCat][2] + the_SelBit_tree->GetEntries();
      the_input_file->Close();
    }
  }
  
  for(int i=0; i<7; i++){
    printf("Input(%d): %d %d %d\n",i,sumInput[i][0],sumInput[i][1],sumInput[i][2]);
    printf("Output(%d): %d %d %d\n",i,sumOutput[i][0],sumOutput[i][1],sumOutput[i][2]);
  }
}
