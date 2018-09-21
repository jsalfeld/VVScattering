from WMCore.Configuration import Configuration
config = Configuration()



config.section_("General")
config.General.requestName = 'nanoTestdxxy'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_scriptdata.sh'
config.JobType.inputFiles = ['./crab_scriptdata.py','../scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
config.Data.inputDataset = '/DoubleMuon/legianni-RunII2017ReReco17Nov17-94X-Nano01Run2017B-17Nov2017-v1-2da645888925dc2add5bb2d3f8d11271/USER'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 2000
config.Data.inputDBS='phys03'
config.Data.outLFNDirBase = '/store/user/jsalfeld/DATAnano2/'
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestPost7'
#config.Data.ignoreLocality = False
config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.whitelist = ['T2_US_*']

#if 'RunII2017ReReco' in config.Data.inputDataset:
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
