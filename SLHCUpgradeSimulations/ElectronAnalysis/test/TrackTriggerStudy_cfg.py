# Import configurations
import FWCore.ParameterSet.Config as cms

# set up process
process = cms.Process("L1EG")

#from SLHCUpgradeSimulations.L1TrackTriggerObjects.minBiasFiles_p1_cfi import *


Source_Files = cms.untracked.vstring(
  'root://eoscms//store/group/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Electrons/PU140/m1_SingleElectron_E2023TTI_PU140.root'
)     
# Event Source
process.source = cms.Source("PoolSource",
     fileNames = Source_Files
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
        #reportEvery = cms.untracked.int32(500),
        reportEvery = cms.untracked.int32(10),
            limit = cms.untracked.int32(10000000)
        )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

# -- Magnetic Field :
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
# -- Geometry
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

# ---------------------------------------------------------------------------
#
# --- Produces the L1EG objects
#
# To produce L1EG objects corresponding to the "stage-2" algorithms:
# one runs the SLHCCaloTrigger sequence. This produces both the
# "old stage-2" objects (2x2 clustering) and the "new stage-2"
# objects (new clustering from JB Sauvan et al). Note that the
# efficiency of the latter is currently poor at very high PU.

process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')

process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")


# The sequence SLHCCaloTrigger is broken in SLHC8 because of the
# L1Jets stuff (should be fixed in SLHC9).
# Hence, in the "remove" lines below, I hack the sequence such that
# we can run the L1EGamma stuff.

process.SLHCCaloTrigger.remove( process.L1TowerJetProducer )
process.SLHCCaloTrigger.remove( process.L1TowerJetFilter1D )
process.SLHCCaloTrigger.remove( process.L1TowerJetFilter2D )
process.SLHCCaloTrigger.remove( process.L1TowerJetPUEstimator )
process.SLHCCaloTrigger.remove( process.L1TowerJetPUSubtractedProducer )
process.SLHCCaloTrigger.remove( process.L1CalibFilterTowerJetProducer )
process.SLHCCaloTrigger.remove( process.L1CaloJetExpander )

# bug fix for missing HCAL TPs in MC RAW
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis = cms.InputTag("valHcalTriggerPrimitiveDigis")

process.slhccalo = cms.Path( process.RawToDigi + process.valHcalTriggerPrimitiveDigis+process.SLHCCaloTrigger)


# To produce L1EG objects corresponding
# to the Run-1 L1EG algorithm, one just needs to run
# L1Reco. The (Run-1) L1EG algorithm has already been
# run in the DIGI step of the production.
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

# To produce BeamSpot
process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")
process.pBeamSpot = cms.Path( process.BeamSpotFromSim)

# Standalone Electron Analysis Code
process.load("SLHCUpgradeSimulations.ElectronAnalysis.TrackTriggerStudy_cfi")
process.trackTriggerStudy.GeometryOption = cms.string("BE5D")
process.trackTriggerStudy.egammaSrc            = cms.InputTag("SLHCL1ExtraParticles","EGamma")
process.trackTriggerStudy.StubPtCut            = cms.double(2.0)
process.trackTriggerStudy.StubEGammadPhiCut    = cms.double(0.5)
process.trackTriggerStudy.StubEGammadZCut      = cms.double(35.0)
process.trackTriggerStudy.StubEGammaPhiMissCut = cms.double(0.01)
process.trackTriggerStudy.StubEGammaZMissCut   = cms.double(3.0)
process.EGamAna = cms.Path( process.trackTriggerStudy)

# Output root file
filename = 'TrackTriggerStudy.root'
print filename
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True))

process.schedule = cms.Schedule(process.slhccalo, process.L1Reco, process.pBeamSpot,process.EGamAna)
