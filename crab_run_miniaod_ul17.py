from CRABClient.UserUtilities import config
config = config()

#config.General.transferLogs = True
config.General.workArea = 'crab/ZeroBias'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_miniaod_ul17.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ["pileup_2016.txt", "pileup_2017.txt", "pileup_2018.txt"]
config.JobType.outputFiles = ["Offset_Data.root"]

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 20
config.Data.unitsPerJob = 200
config.Data.lumiMask = '/afs/cern.ch/user/m/minsuk/public/certificates/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_mod.txt'

config.Data.outLFNDirBase = '/store/user/minsuk/UE'

#config.Site.storageSite = 'T3_KR_KISTI'
#config.Site.storageSite = 'T2_KR_KNU'
#config.Site.storageSite = 'T2_CH_CERNBOX'
config.Site.storageSite = 'T2_FI_HIP'

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    config.General.requestName = 'Run2017Bv1'
    config.Data.inputDataset = '/ZeroBias/Run2017B-09Aug2019_UL2017-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2017Cv1'
    config.Data.inputDataset = '/ZeroBias/Run2017C-09Aug2019_UL2017-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2017Dv1'
    config.Data.inputDataset = '/ZeroBias/Run2017D-09Aug2019_UL2017-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2017Ev1'
    config.Data.inputDataset = '/ZeroBias/Run2017E-09Aug2019_UL2017-v1/MINIAOD'
    crabCommand('submit', config = config)

    config.General.requestName = 'Run2017Fv1'
    config.Data.inputDataset = '/ZeroBias/Run2017F-09Aug2019_UL2017-v1/MINIAOD'
    crabCommand('submit', config = config)

