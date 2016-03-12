import ROOT

ROOT.gSystem.Load('libMitAnaDataTree.so')

source = ROOT.TFile.Open('/mnt/hadoop/cms/store/user/paus/filefi/043/ggZH_HToInv_ZToLL_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM/F84CC9C4-ACC0-E511-AAB4-0025905A48FC.root')
tree = source.Get('HLT')

triggerTable = getattr(ROOT, 'std::vector<std::string>')()

tree.SetBranchAddress('HLTTriggerTable', triggerTable)
branch = tree.GetBranch('HLTTriggerTable')

branch.GetEntry(0)

for path in triggerTable:
    print path
