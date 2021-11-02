from CRABClient.UserUtilities import config
from WMCore.Configuration import Configuration
config = config()

config.General.requestName = 'W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 300
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 9000
config.JobType.psetName = '/uscms_data/d3/alkaloge/MetStudies/CMSSW_10_6_5/src/Wjets/TM/test/crab_MC/treemaker_cfg_mc2018.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['histo.root']
config.JobType.inputFiles = ['Autumn18_V19_MC.db']

config.Data.inputDataset = '/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/group/lpcsusyhiggs/ntuples/MetStudies/2018/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8'

config.Site.storageSite = 'T3_US_FNALLPC'

