import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class wzAnalysisProducer(Module):
    def __init__(self):
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("EventMass",  "F");
	#pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        keepIt = True
	eventLeptons = 0
        for lep in muons :
            	if lep.tightId and abs(lep.dz)<0.1 and abs(lep.dxy < 0.02) and lep.pfRelIso04_all < 0.4:
			eventLeptons += 1	
            
        for lep in electrons :
                if lep.cutBased_HLTPreSel == 1:
			eventLeptons += 1
	if eventLeptons < 3:
		keepIt = False
        
        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

wzAnalysisModule = lambda : wzAnalysisProducer() #(jetSelection= lambda j : j.pt > 30) 
 
