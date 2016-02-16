import ROOT

ROOT.gSystem.Load('libMitAnaDataTree.so')

source = ROOT.TFile.Open('/mnt/hadoop/cms/store/user/paus/filefi/042/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM/001A0168-0914-E511-A04A-0025905A48C0.root')
tree = source.Get('HLT')

triggerTable = getattr(ROOT, 'std::vector<std::string>')()

tree.SetBranchAddress('HLTTriggerTable', triggerTable)
branch = tree.GetBranch('HLTTriggerTable')

branch.GetEntry(0)

for path in triggerTable:
    print path
