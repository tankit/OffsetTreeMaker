from CRABClient.UserUtilities import config
config = config()

#config.General.transferLogs = True
config.General.workArea = 'crab/ZeroBias'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_miniaod_ul16.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ["pileup_2016.txt", "pileup_2017.txt", "pileup_2018.txt"]
config.JobType.outputFiles = ["Offset_Data.root"]

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 20
config.Data.unitsPerJob = 200
config.Data.lumiMask = '/afs/cern.ch/user/m/minsuk/public/certificates/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'

config.Data.outLFNDirBase = '/store/user/minsuk/UE'

#config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T2_KR_KNU'
#config.Site.storageSite = 'T2_CH_CERNBOX'
config.Site.storageSite = 'T2_FI_HIP'

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    ##--no lumi to process (a golden json lumimask does not include those runs)
    ##config.General.requestName = 'Run2016B_ver1_v1'
    ##config.Data.inputDataset = '/ZeroBias/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD'
    ##crabCommand('submit', config = config)

    config.General.requestName = 'Run2016B_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016C_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016D_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016E_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016Fe_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016Fl_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016F-21Feb2020_UL2016-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016G_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016G-21Feb2020_UL2016-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2016H_v1'
    config.Data.inputDataset = '/ZeroBias/Run2016H-21Feb2020_UL2016-v1/MINIAOD'
    crabCommand('submit', config = config)
