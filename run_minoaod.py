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
  #'/store/data/Run2017B/ZeroBias/MINIAOD/17Nov2017-v1/50000/FAA4C8D3-D4D3-E711-A750-484D7E8DF114.root' 
  '/store/mc/RunIIFall17MiniAODv2/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/FE935FB1-DD44-E811-B398-0CC47AA53D66.root'
] );

#isMC = cms.bool(False)
isMC = cms.bool(True)

if isMC:
  OutputName = "_MC"
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '94X_mc2017_realistic_v14' )

else:
  OutputName = "_Data"

  process.load( "Configuration.Geometry.GeometryIdeal_cff" )
  process.load( "Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
  process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff" )
  from Configuration.AlCa.GlobalTag import GlobalTag
  process.GlobalTag = GlobalTag( process.GlobalTag, '94X_dataRun2_ReReco17_forValidation' )

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
    numSkip = cms.int32(101),
    RootFileName = cms.string("Offset" + OutputName + ".root"),
    puFileName = cms.string("lumi-per-bx.root"),
    isMC = isMC,
    writeCands = cms.bool(True),
    #trackTag = cms.InputTag("generalTracks"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    genTag = cms.InputTag("packedGenParticles"),
    pfTag = cms.InputTag("packedPFCandidates"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muTag = cms.InputTag("slimmedAddPileupInfo"),
    rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoC0Tag = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    rhoCCTag = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
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
