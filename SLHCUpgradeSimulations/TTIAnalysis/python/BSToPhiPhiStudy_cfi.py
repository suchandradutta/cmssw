import FWCore.ParameterSet.Config as cms
bsStudy = cms.EDAnalyzer("BSToPhiPhiStudy",
     trkSrc      = cms.InputTag("TTTracksFromPixelDigis","Level1TTTracks"),
     trkTruthSrc = cms.InputTag("TTTrackAssociatorFromPixelDigis","Level1TTTracks"),
     pixelTrkSrc = cms.InputTag("L1PixelTrackFit","Level1PixelTracks"),
     DebugFlag   = cms.bool(False)
)
