# PYTHON configuration file for class: OffsetTreeMaker
# Author: C. Harrington
# Updated by Minsuk Kim for MINIAOD
# Date:  19 - February - 2018

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

readFiles = cms.untracked.vstring()
process.source = cms.Source ("PoolSource", fileNames = readFiles)
readFiles.extend( [
  #'/store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/280000/54C32B50-BAA7-0F4A-AE88-8742F92B46A8.root'
  #'/store/mc/RunIISummer19UL17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/260000/01E437E0-B757-5B45-950F-CF648C57CE32.root'
  #'/store/data/Run2017B/ZeroBias/MINIAOD/09Aug2019_UL2017-v1/260000/0228BC4C-26B0-C04E-8BD3-349FAEC6AABF.root'
  #'/store/data/Run2018A/ZeroBias/MINIAOD/12Nov2019_UL2018-v2/100000/032F9DB5-A22C-5246-A1B3-E0EB0C539A8E.root'
  '/store/data/Run2016B/ZeroBias/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/109E3350-A8F4-6E4F-B4EE-B07695A21179.root' 
] );

isMC = cms.bool(False)
#isMC = cms.bool(True)

if isMC:
  OutputName = "_MC"
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  #process.GlobalTag = GlobalTag( process.GlobalTag, '106X_mc2017_realistic_v6' )
  process.GlobalTag = GlobalTag( process.GlobalTag, '106X_mc2017_realistic_v9' )

else:
  OutputName = "_Data"

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  #process.GlobalTag = GlobalTag( process.GlobalTag, '106X_dataRun2_v20' ) #UL2017 B-F
  #process.GlobalTag = GlobalTag( process.GlobalTag, '106X_dataRun2_v24' ) #UL2018 A-C
  #process.GlobalTag = GlobalTag( process.GlobalTag, '106X_dataRun2_v26' ) #UL2018 D
  process.GlobalTag = GlobalTag( process.GlobalTag, '106X_dataRun2_v27' ) #UL2016 B-H

  # ZeroBias Trigger
  process.HLTZeroBias =cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring('HLT_ZeroBias_part*','HLT_ZeroBias_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw = cms.bool(False)
  )

  #Beam Halo
  process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

  #HCAL HBHE
  process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
  process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
  process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Tight'),
    reverseDecision = cms.bool(False)
  )

process.pf = cms.EDAnalyzer("OffsetTreeMaker",
    numSkip = cms.int32(1),
    RootFileName = cms.string("Offset" + OutputName + ".root"),
    #puFileName = cms.string("lumi-per-bx.root"),
    puFileName = cms.string("pileup_2016.txt"),
#    puFileName = cms.string("pileup_2017.txt"),
#    puFileName = cms.string("pileup_2018.txt"),
    isMC = isMC,
    writeCands = cms.bool(False),
    writeParticles = cms.bool(False),
    #trackTag = cms.InputTag("generalTracks"),
    Generator = cms.InputTag("generator"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    genTag = cms.InputTag("packedGenParticles"),
    pfTag = cms.InputTag("packedPFCandidates"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muTag = cms.InputTag("slimmedAddPileupInfo"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoC0Tag = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    rhoCCTag = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    rhoCentralTag = cms.InputTag("fixedGridRhoFastjetCentral"),
    rhoCentralCaloTag = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    pfJetTag = cms.InputTag("slimmedJets")
)

process.myseq = cms.Sequence( process.pf )

if isMC :
  process.p = cms.Path( process.myseq )
else:
  process.p = cms.Path( process.HLTZeroBias * 
                        process.CSCTightHaloFilter *
                        process.HBHENoiseFilterResultProducer *
                        process.ApplyBaselineHBHENoiseFilter *
                        process.myseq )
