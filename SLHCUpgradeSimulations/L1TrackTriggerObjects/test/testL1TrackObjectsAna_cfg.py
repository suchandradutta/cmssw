import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:example.root'
    #'file:example_withEtMiss.root'
    #'file:example_w_Tracks.root'
    )
)




process.ana = cms.EDAnalyzer( 'L1TrackTriggerObjectsAnalyzer' ,
    L1VtxInputTag = cms.InputTag("L1TrackPrimaryVertex"),
    L1EtMissInputTag = cms.InputTag("L1TrackEtMiss","MET"),
    L1TrackElectronsInputTag = cms.InputTag("L1TrackElectrons","NonIsolated"),
    L1TrackPhotonsInputTag = cms.InputTag("L1TrackPhotons","NonIsolated")
)


# root file with histograms produced by the analyzer
filename = "ana.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path( process.ana )




