import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:../../../samples/step2.root'
    )
)

process.demo = cms.EDAnalyzer('ME0analyzer',
    me0DigiToken = cms.InputTag("simMuonME0Digis"),
    me0PadDigiToken = cms.InputTag("simMuonME0PadDigis"),
    me0PadDigiClusterToken = cms.InputTag("simMuonME0PadDigiClusters"),
    genParticles = cms.InputTag("genParticles")
)


process.p = cms.Path(process.demo)
