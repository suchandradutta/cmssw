import FWCore.ParameterSet.Config as cms
process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
    limit = cms.untracked.int32(10000000)
)     
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# list of files
file_names = cms.untracked.vstring(
  '/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/BsToPhiPhi4K/PU140/Digi_BsToPhiPhi4K_PU140_1.root'
)     
process.source = cms.Source("PoolSource",
            fileNames = file_names
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

############################################################
# import standard configurations
############################################################
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')

# Geometry
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')

## pixel additions
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

# ---------------------------------------------------------------------------
# ---  Run the SLHCCaloSequence  to produce the L1Jets
process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_forTTI_cff")

#process.L1CalibFilterTowerJetProducer.pTCalibrationThreshold = cms.double(40) # applies calibration only to > 40GeV L1 jets

# ---------------------------------------------------------------------------
#
# --- Create the collection of L1 Objects 
#
# L1 Track
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TrackingSequence_cfi")
process.pTracking = cms.Path( process.DefaultTrackingSequence )

# L1 Vertex 
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkPrimaryVertexProducer_cfi")
process.pL1TkPrimaryVertex = cms.Path( process.L1TkPrimaryVertex )

# L1 Muon 
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkMuonSequence_cfi")
process.pMuons = cms.Path( process.L1TkMuons )

# L1 Jet
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkJetProducer_cfi")
process.L1TkJetsL1 = process.L1TkJets.clone()
process.pL1TkJetsL1 = cms.Path( process.L1TkJetsL1 )

BeamSpotFromSim = cms.EDProducer("VtxSmeared")
process.pBeamSpot = cms.Path( process.BeamSpotFromSim )

# pixel stuff
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
process.siPixelRecHits = siPixelRecHits

process.L1PixelTrackFit = cms.EDProducer("L1PixelTrackFit")
process.pixTrk = cms.Path(process.L1PixelTrackFit)

process.myreco = cms.Path(
    process.RawToDigi+
    process.valHcalTriggerPrimitiveDigis+
    process.SLHCCaloTrigger+
    process.siPixelRecHits
)
process.es_prefer_dt = cms.ESPrefer("DTConfigTrivialProducer","L1DTConfig")
# bug fix for missing HCAL TPs in MC RAW
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("valHcalTriggerPrimitiveDigis")

process.raw2digi_step = cms.Path(process.RawToDigi)                                                              

# Tree Producer Code Standalone BstoPhiPhi4k Analysis Code
process.load("SLHCUpgradeSimulations.TTIAnalysis.BSToPhiPhiStudy_cfi")
process.bsStudy.l1JetFlag       = cms.bool(True)
process.bsStudy.l1MuonFlag      = cms.bool(True)
process.bsStudy.recoTrackFlag   = cms.bool(False)
process.bsStudy.genParticleFlag = cms.bool(True)
#process.bsStudy.l1TkMuonSrc     = cms.InputTag("L1TkMuonsMerge", "")
process.BSAna = cms.Path(process.bsStudy)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.eventDump = cms.Path(process.dump)

# Output root file
filename = 'test_L1TkMuon.root'
print filename
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True))

process.schedule = cms.Schedule(process.pBeamSpot, 
                                process.pTracking, 
                                process.pL1TkPrimaryVertex, 
                                process.myreco, 
                                process.pixTrk, 
                                process.pL1TkJetsL1, 
                                process.pMuons, 
                                process.BSAna)

import sys
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI

# call to customisation function cust_2023TTI imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023TTI(process)
