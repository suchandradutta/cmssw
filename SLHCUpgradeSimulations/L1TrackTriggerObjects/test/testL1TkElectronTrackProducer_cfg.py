
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
  'root://eoscms//store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root',
  'root://eoscms//store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/027029F2-FE38-E311-ACD6-003048678B34.root'
# '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00D6C34E-0339-E311-836A-002618943880.root',
# '/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FEDB1C0F-FF38-E311-A659-0025905938D4.root'
)
# input Events 
process.source = cms.Source("PoolSource",
   fileNames = file_names,
   skipEvents = cms.untracked.uint32(0) 
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

# L1Tracking 
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')

process.BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.BeamSpotFromSim*process.L1Tracks)


# ---------------------------------------------------------------------------
#
# --- Run the L1EG algorithm of Jean-Baptiste (new clustering algorithm).
# --- Note thet the efficiency is poor at very high PU...
#
# --- This also runs the "old" stage-2 algorithm

process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.pSLHCCalo = cms.Path(
    process.RawToDigi+
    process.SLHCCaloTrigger
)
# bug fix for missing HCAL TPs in MC RAW
process.pSLHCCalo.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode    = cms.bool(True)
process.valRctDigis.hcalDigis         = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

# run L1Reco to produce the L1EG objects corresponding
# to the current trigger
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

# ---------------------------------------------------------------------------
#
# --- test L1TrackEmParticle


# --- Run the L1TkEmParticleProducer

# "photons" :
# ---------------------------------------------------------------------------
#
# --- test L1TrackEmParticle


# --- Run the L1TkEmParticleProducer

# "electrons" :

process.L1TkElectrons = cms.EDProducer("L1TkElectronTrackProducer",
	label = cms.string("ElecTrk"),	# labels the collection of L1TkEmParticleProducer that is produced.
                                        # e.g. Elec or IsoElec if all objects are kept, or
                                        # ElecIsoTrk or IsoElecIsoTrk if only the EG or IsoEG
                                        # objects that pass a cut RelIso < RelIsoCut are written
                                        # into the new collection.
        L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","EGamma"),      # input EGamma collection
					# When the standard sequences are used :
					#   - for "old stage-2", use ("l1extraParticles","NonIsolated")
					#     or ("l1extraParticles","Iolated")
					#   - for the new clustering algorithm of Jean-Baptiste et al,
					#     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
					#     ("SLHCL1ExtraParticlesNewClustering","EGamma").                                      
        ETmin = cms.double( 2.0 ),       # Only the L1EG objects that have ET > ETmin in GeV
        TrackEGammaDeltaPhi = cms.double(0.1),  # Delta Phi cutoff to match Track with L1EG objects
        TrackEGammaDeltaR = cms.double(0.06),   # Delta R cutoff to match Track with L1EG objects
        TrackEGammaDeltaEta = cms.double(0.05), # Delta Eta cutoff to match Track with L1EG objects
                                                # are considered. ETmin < 0 means that no cut is applied.
	RelativeIsolation = cms.bool( True ),	# default = True. The isolation variable is relative if True,
						# else absolute.
        IsoCut = cms.double( -0.15 ), 		# Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
						# When RelativeIsolation = False, IsoCut is in GeV.
        # Determination of the isolation w.r.t. L1Tracks :
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	ZMAX = cms.double( 25. ),	# in cm
	CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 5. ),	# in GeV
	DRmin = cms.double( 0.06),
	DRmax = cms.double( 0.5 ),
	DeltaZ = cms.double( 1.0 )    # in cm. Used for tracks to be used isolation calculation
)
process.pElectrons = cms.Path( process.L1TkElectrons )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "L1TrackElectron.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_SLHCL1ExtraParticles_EGamma_*' )
process.Out.outputCommands.append( 'keep *_L1TkElectrons_ElecTrk_*' )
process.Out.outputCommands.append( 'keep SimTracks_g4SimHits_*_*'), 
#process.Out.outputCommands.append('keep *_generator_*_*')
#process.Out.outputCommands.append('keep *')

#process.schedule = cms.Schedule(process.p0,process.L1Reco,process.TT_step,process.pElectrons)
process.FEVToutput_step = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.pSLHCCalo,process.L1Reco,process.TT_step,process.pElectrons,process.FEVToutput_step)




