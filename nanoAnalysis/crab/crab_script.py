#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from  PhysicsTools.NanoAODTools.postprocessing.wzAnalysis.wzAnalysisModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *


p=PostProcessor(".",inputFiles(),"",modules=[wzAnalysisProducer(),jetmetUncertainties2017(),puWeightProducer("auto",pufile_data2017,"pu_mc","pileup",verbose=False)],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
p.run()

print "DONE"
os.system("ls -lR")
