import FWCore.ParameterSet.Config as cms

# Example configuratiom file that runs over a file that has
# already the L1Tracks.
# (see e.g. runL1Tracks_and_L1PrimaryVertexProducer_cfg.py for how to produce them).

# Here we run the L1EG algorithms (old stage-2 and new clustering),
# we unpack the L1EG objects that were created during the L1 step
# of the central production (i.e. the Run-1 algorithms), and we
# create L1TkEm objects corresponding to the various input
# collections.


process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

from SLHCUpgradeSimulations.L1TrackTriggerObjects.singleElectronFiles_cfi import *

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
   'file:example_w_Tracks_and_vertex.root'
    )
)



# ---- Global Tag and geometry :
#      (needed e.g. when running raw2digi below)

process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

# ---------------------------------------------------------------------------
#
# --- to work only with the Run-1 L1EG algorithm, i.e. what was 
# --- created during the central production when running the L1 step :

        # raw2digi to get the gct digis
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.p0 = cms.Path( process.RawToDigi )
        # run L1Reco
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.L1Reco = cms.Path( process.l1extraParticles )


# ---------------------------------------------------------------------------
#
# --- Run the L1EG algorithm of Jean-Baptiste (new clustering algorithm).
# --- Note thet the efficiency is poor at very high PU...
#
# --- This also runs the "old" stage-2 algorithm

process.load('Configuration/StandardSequences/L1HwVal_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

process.p = cms.Path(
    process.RawToDigi+
    process.SLHCCaloTrigger
    )

# bug fix for missing HCAL TPs in MC RAW
process.p.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
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

	# --- from the Run-1 EG algorithms...

process.L1TkPhotonsRun1EG = cms.EDProducer("L1TkEmParticleProducer",
        label = cms.string("IsoTrk"), # labels the collection of L1TkEmParticleProducer that is produced.
                                                # e.g. EG or IsoEG if all objects are kept, or
                                                # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                                # objects that pass a cut RelIso < RelIsoCut are written
                                                # into the new collection.
        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),      # input L1EG collection
                                                # When the standard sequences are used :
                                                #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
						#   - for the "old stage-2" algo (2x2 clustering), use 
						#     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
        ETmin = cms.double( -1 ),               # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
        RelativeIsolation = cms.bool( True ),   # default = True. The isolation variable is relative if True,
                                                # else absolute.
        IsoCut = cms.double( 0.1 ),             # Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
                                                # When RelativeIsolation = False, IsoCut is in GeV.
           # Determination of the isolation w.r.t. L1Tracks :
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
        ZMAX = cms.double( 25. ),       # in cm
        CHI2MAX = cms.double( 100. ),
        PTMINTRA = cms.double( 2. ),    # in GeV
        DRmin = cms.double( 0.07),
        DRmax = cms.double( 0.30 ),
        PrimaryVtxConstrain = cms.bool( False ),  # default = False
        DeltaZMax = cms.double( 999. ),    # in cm. Used only when PrimaryVtxConstrain = True
        L1VertexInputTag = cms.InputTag("NotUsed"),     # Used only when PrimaryVtxConstrain = True
)
process.pPhotonsRun1EG = cms.Path( process.L1TkPhotonsRun1EG )

process.L1TkPhotonsRunIso1EG = process.L1TkPhotonsRun1EG.clone()
process.L1TkPhotonsRunIso1EG.L1EGammaInputTag = cms.InputTag("l1extraParticles","Isolated")
process.pPhotonsRun1IsoEG = cms.Path( process.L1TkPhotonsRunIso1EG )

	# --- old stage -2 :

process.L1TkPhotonsStage2EG = process.L1TkPhotonsRun1EG.clone()
process.L1TkPhotonsStage2EG.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","EGamma")
process.pPhotonsStage2EG = cms.Path( process.L1TkPhotonsStage2EG )

process.L1TkPhotonsStage2IsoEG = process.L1TkPhotonsRun1EG.clone()
process.L1TkPhotonsStage2IsoEG.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticles","IsoEGamma")
process.pPhotonsStage2IsoEG = cms.Path( process.L1TkPhotonsStage2IsoEG )

	# --- new clustering (Jean Baptiste)

process.L1TkPhotonsNewEG = process.L1TkPhotonsRun1EG.clone()
process.L1TkPhotonsNewEG.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma")
process.pPhotonsNewEG = cms.Path( process.L1TkPhotonsNewEG )

process.L1TkPhotonsNewIsoEG = process.L1TkPhotonsRun1EG.clone()
process.L1TkPhotonsNewIsoEG.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","IsoEGamma")
process.pPhotonsNewIsoEG = cms.Path( process.L1TkPhotonsNewIsoEG )




process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_EGamma.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)


process.Out.outputCommands.append('keep *_L1TkPhotons*_*_*')
process.Out.outputCommands.append('keep *_L1TkEtMiss_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_NonIsolated_*')
process.Out.outputCommands.append('keep *_l1extraParticles_Isolated_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_EGamma_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticles_IsoEGamma_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_EGamma_*')
process.Out.outputCommands.append('keep *_SLHCL1ExtraParticlesNewClustering_IsoEGamma_*')
process.Out.outputCommands.append('keep *_l1extraParticles_MET_*')
process.Out.outputCommands.append('keep *_gen*_*_*')


process.FEVToutput_step = cms.EndPath(process.Out)




