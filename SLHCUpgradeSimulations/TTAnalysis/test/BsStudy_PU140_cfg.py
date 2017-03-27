import FWCore.ParameterSet.Config as cms
process = cms.Process("Ele")

process.load("FWCore.MessageService.MessageLogger_cfi")

################################################################################
# Example configuratiom file : 
# Here we run the L1EG algorithms (old stage-2 and new clustering),
# we unpack the L1EG objects that were created during the L1 step
# of the central production (i.e. the Run-1 algorithms), and we
# create L1TkEm objects corresponding to the various input
# collections.                                                                            
################################################################################
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
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('Configuration.StandardSequences.RawToDigi_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V3::All', '')

## pixel additions
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')


# ---------------------------------------------------------------------------
#
# --- Create the collection of L1 Objects 
#
# Track
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TrackingSequence_cfi")
process.pTracking = cms.Path( process.DefaultTrackingSequence )
# Muon 
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkMuonSequence_cfi")
process.pMuons = cms.Path( process.L1TkMuons )
# Jet
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkJetProducer_cfi")
process.L1TkJetsL1 = process.L1TkJets.clone()
process.pL1TkJetsL1 = cms.Path( process.L1TkJetsL1 )
# Vertex 
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TkPrimaryVertexProducer_cfi")
process.pL1TkPrimaryVertex = cms.Path( process.L1TkPrimaryVertex )

BeamSpotFromSim =cms.EDProducer("VtxSmeared")
process.pBeamSpot = cms.Path( process.BeamSpotFromSim )

#pixel stuff
from RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi import *
process.siPixelRecHits = siPixelRecHits

process.L1PixelTrackFit = cms.EDProducer("L1PixelTrackFit")
process.pixTrk = cms.Path(process.L1PixelTrackFit)

process.pixRec = cms.Path(
    process.RawToDigi+
    process.siPixelRecHits
)
process.raw2digi_step = cms.Path(process.RawToDigi)
                                                              

# Tree Producer Code Standalone BstoPhiPhi4k Analysis Code
process.load("SLHCUpgradeSimulations.TTAnalysis.BSToPhiPhiStudy_cfi")
process.bsStudy.l1TrackFlag    = cms.bool(True)
process.bsStudy.pixelTrackFlag = cms.bool(True)
process.bsStudy.recoTrackFlag   = cms.bool(False)
process.bsStudy.genParticleFlag = cms.bool(True)

process.BSAna = cms.Path( process.bsStudy)


process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.eventDump = cms.Path(process.dump)
# Output root file
filename = 'Signal_BsStudy_PU140_1.root'
print filename
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename), closeFileFast = cms.untracked.bool(True))

#New Tracking (no pixel)
#process.schedule = cms.Schedule(process.pBeamSpot, process.pTracking, process.pL1TkPrimaryVertex, process.BSAna)
#New Tracking (with pixel)
process.schedule = cms.Schedule(process.pBeamSpot, process.pTracking, process.pL1TkPrimaryVertex, process.pixRec,process.pixTrk, process.pMuons, process.pL1TkJetsL1, process.BSAna)

import sys
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023TTI

#call to customisation function cust_2023TTI imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023TTI(process)

