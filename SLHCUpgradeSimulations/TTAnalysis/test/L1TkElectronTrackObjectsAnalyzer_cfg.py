import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")

file_names = cms.untracked.vstring(
#  'file:/afs/cern.ch/work/d/dutta/public/TkUpgrade/analysis/FullSim/PU140/pid11_prod_newStage2_62X/default/Electron_1.root'
  'file:/afs/cern.ch/user/d/dutta/work/public/TkUpgrade/analysis/FullSim/PU140/pid11_prod_newStage2_62X/degradeZ/WToEnu_1.root'
#   'file:/afs/cern.ch/work/d/dutta/public/TkUpgrade/software/CMSSW_6_2_0_SLHC12/src/SLHCUpgradeSimulations/TTAnalysis/test/WToEnu_1.root'
)     

from SLHCUpgradeSimulations.TTAnalysis.TT_WENeu_PU200_degradeZ_cfi import *
process.source = cms.Source("PoolSource",
                  fileNames = wenuFiles  # the private files are better as the official sample
#                 fileNames = file_names
                  
              )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.TrkAna = cms.EDAnalyzer( 'L1TkElectronTrackObjectsAnalyzer' ,
    L1TkElectronsInputTag = cms.InputTag("L1TkEGamma","EG"),
    L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
    AnalysisOption   = cms.string("Efficiency"),
    EtaCutOff   = cms.double(2.5),
    TrackPtCutOff   = cms.double(10.0),
    GenPtThreshold   = cms.double(20.0),
    EGammaEtThreshold = cms.double(20.0)                              
)
process.IsoTrkAna = cms.EDAnalyzer( 'L1TkElectronTrackObjectsAnalyzer' ,
    L1TkElectronsInputTag = cms.InputTag("L1IsoTkEGamma","EG"),
    L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
    AnalysisOption   = cms.string("Efficiency"),
    EtaCutOff   = cms.double(2.5),
    TrackPtCutOff   = cms.double(10.0),
    GenPtThreshold   = cms.double(20.0),
    EGammaEtThreshold = cms.double(20.0)                              
)


# root file with histograms produced by the analyzer
filename = "efficiency_degradeZ_PU200.root"
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path( process.TrkAna * process.IsoTrkAna )




