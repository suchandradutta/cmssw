import FWCore.ParameterSet.Config as cms

process = cms.Process("ETMISS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:example.root'
    '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/zmatchingOff/m1_TTbar_BE5D.root'
    )
)


# ---- Global Tag :
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')

        # run L1Reco to produce the L1ETM object corresponding
        # to the current trigger
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )

	# now produce the TrkMET :

process.L1TkEtMiss = cms.EDProducer('L1TkEtMissProducer',
     L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
     L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex"),
     ZMAX = cms.double ( 25. ) ,        # in cm
     CHI2MAX = cms.double( 100. ),
     PTMINTRA = cms.double( 2. ),       # in GeV
     DeltaZ = cms.double( 0.2 ),       # in cm
     nStubsmin = cms.int32( 4 )
)

process.p = cms.Path( process.L1TkEtMiss )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_withEtMiss.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append( 'keep *_*_*_ETMISS' )
process.Out.outputCommands.append('keep *_generator_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




