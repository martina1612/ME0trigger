import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/m/mressego/private/ME0trigger/CMSSW_9_3_3/src/samples/step2.root'
    )
)

process.demo = cms.EDAnalyzer('ME0analyzer',
    verbose = cms.bool(True),
    me0DigiToken = cms.InputTag("simMuonME0Digis"),
    me0PadDigiToken = cms.InputTag("simMuonME0PadDigis"),
    me0PadDigiClusterToken = cms.InputTag("simMuonME0PadDigiClusters"),
    genParticles = cms.InputTag("genParticles")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ME0analyzer.root'))

process.p = cms.Path(process.demo)
