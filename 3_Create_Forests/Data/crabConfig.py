from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.section_('JobType')
config.JobType.psetName = 'step2_DIGI_L1_DIGI2RAW_HLT_L1Reco.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HIMuonHLT_DoubleSingleMu_WZQuarkonia_V48.root']
#config.JobType.maxJobRuntimeMin = 1900
config.section_('Data')
config.Data.inputDataset = '/Hydjet_Quenched_MinBias_5020GeV/HiFall14-START71_V1-v2/GEN-SIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS=100000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/ilaflott/trigger_work/HIPaths_DoubleMu_WZJpsiUpsilon_V42/5020GeVHydjetQuenchedMinBias_Test_6.11.15_allEvents'
config.section_('Site')
config.Site.whitelist = ["T2*"]
config.Site.storageSite = 'T2_US_MIT'
