import FWCore.ParameterSet.Config as cms

process = cms.Process("ALL")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

from SLHCUpgradeSimulations.L1TrackTriggerObjects.singleElectronFiles_cfi import *

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/MinBias/BE5D/zmatchingOff/m1_MinBias_BE5D.root',
       '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/MinBias/BE5D/zmatchingOff/m2_MinBias_BE5D.root'
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
# --- Run the L1TrackEtMiss producer :

# needs the L1TkPrimaryVertex
# TrkMET calculated from all tracks that have |z - z_vertex | < DeltaZ,
# pt > PTMINTRA and at least ( >= ) nStubsmin stubs.
# To check the latter condition, the file must contain the stubs.

process.L1TkEtMiss = cms.EDProducer('L1TkEtMissProducer',
     L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
     L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex2"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),
     PTMINTRA = cms.double( 2. ),	# in GeV
     DeltaZ = cms.double( 0.2 ),       # in cm
     nStubsmin = cms.int32( 4 )
)

process.pEtMiss = cms.Path( process.L1TkEtMiss )

#
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#
# --- test L1TrackEmParticle

# --- First, create l1extra objects for L1EGamma 
#     here I use the Run1 trigger. The produced objects dont really
#     make sense, but OK to check the technicalities.

        # raw2digi to get the gct digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )
        # run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )



# --- Now run the L1TkEmParticleProducer 

# "photons" :

process.L1TkPhotons = cms.EDProducer("L1TkEmParticleProducer",
        label = cms.string("EGIsoTrk"), # labels the collection of L1TkEmParticleProducer that is produced.
                                                # e.g. EG or IsoEG if all objects are kept, or
                                                # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                                # objects that pass a cut RelIso < RelIsoCut are written
                                                # into the new collection.
        L1EGammaInputTag = cms.InputTag("l1extraParticles","NonIsolated"),      # input L1EG collection
                                                # When the standard sequences are used :
                                                #   - for "old stage-2", use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
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
process.pPhotons = cms.Path( process.L1TkPhotons )




process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_all.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)


process.Out.outputCommands.append('keep *_L1TkPhotons_*_*')
process.Out.outputCommands.append('keep *_L1TkEtMiss_*_*')
process.Out.outputCommands.append('keep *_l1extraParticles_*_*')


process.FEVToutput_step = cms.EndPath(process.Out)




