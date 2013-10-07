import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("EGamma")

import os

Source_Files = cms.untracked.vstring(
  'root://eoscms//store/mc/Summer13/SingleElectronFlatPt5To50/GEN-SIM-DIGI-RAW/UpgradePhase2BE_2013_DR61SLHCx_PU140Bx25_POSTLS261_V2-v1/10000/00F37080-50C8-E211-B3B3-003048FFCBB0.root'
)     
# Event Source
process.source = cms.Source("PoolSource",
     fileNames = Source_Files
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Basic running parameters (modify to your needs)
index=44

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


# L1 Tracking
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBEReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE_cff')
#   (or 'LB_4PS_2SS' or 'EB'  for the other geometries.)
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")
process.L1Tracks.geometry = cms.untracked.string('BE')
#   (or 'LB_4PS_2SS' or 'EB'  for the other geometries.)
process.pL1Tracks = cms.Path( process.BeamSpotFromSim*process.L1Tracks )
                                        

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V2::All', '')


#MANU:
process.load('Configuration.StandardSequences.L1Reco_cff')

process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.L1CaloTowerProducer.ECALDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")

process.pL1Calo = cms.Path( process.SLHCCaloTrigger)


# -- to produce the l1extra muon particles :
process.l1extraParticles.muonSource = cms.InputTag("simGmtDigis")
process.l1extraParticles.produceCaloParticles = cms.bool(False)
process.L1Reco = cms.Path( process.l1extraParticles )

# Produce L1 EG objects with precise position from Crystal Info
process.L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
   DEBUG = cms.bool(False)
)
process.l1ExtraCrystalProducer = cms.EDProducer("L1ExtraCrystalPosition",
    eGammaSrc = cms.InputTag("SLHCL1ExtraParticles","EGamma"),
    eClusterSrc = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster")
)
process.crystal_producer = cms.Path(process.L1EGammaCrystalsProducer)

process.egcrystal_producer = cms.Path(process.l1ExtraCrystalProducer)

#process.Timing=cms.Service("Timing")

#Add these 3 lines to put back the summary for timing information at the end of the logfile
#(needed for TimeReport report)
#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True)
#    )

# Output root file
filename = 'Electron_EGammaN_BE.root'
print filename
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True))

process.load("SLHCUpgradeSimulations.ElectronAnalysis.TrackTriggerStudy_cfi")
process.trackTriggerStudy.GeometryOption = cms.string("BE")
# default L1EGamma objects       
#process.trackTriggerStudy.egammaSrc            = cms.InputTag("SLHCL1ExtraParticles","EGamma")
#1EGamma objects with Crystal level position 
process.trackTriggerStudy.egammaSrc            = cms.InputTag("l1ExtraCrystalProducer","EGammaCrystal")
process.trackTriggerStudy.StubPtCut            = cms.double(3.0)
process.trackTriggerStudy.StubEGammadPhiCut    = cms.double(0.5)
process.trackTriggerStudy.StubEGammadZCut      = cms.double(35.0)
process.trackTriggerStudy.StubEGammaPhiMissCut = cms.double(0.005)
process.trackTriggerStudy.StubEGammaZMissCut   = cms.double(3.0)

process.EGamAna = cms.Path( process.trackTriggerStudy)


process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path( process.calolocalreco )

# Schedule with default L1EGamma objects
#process.schedule = cms.Schedule(process.pL1Calo,process.L1Reco,process.pL1Tracks,process.EGamAna)
# Schedule with L1EGamma objects with Crystal level position
process.schedule = cms.Schedule(process.raw2digi_step, process.reconstruction_step, process.pL1Calo,process.L1Reco,process.crystal_producer,process.egcrystal_producer,process.pL1Tracks,process.EGamAna)



