import FWCore.ParameterSet.Config as cms

process = cms.Process("ETMISS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:example.root'
    '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/m1_TTbar_BE5D.root'
    )
)


# ---- Global Tag :
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')



process.L1TrackEtMiss = cms.EDProducer('L1TrackEtMissProducer',
     L1VtxLabel = cms.InputTag("L1TrackPrimaryVertex"),
     ZMAX = cms.double ( 25. ) ,	# in cm
     CHI2MAX = cms.double( 100. ),
     DeltaZ = cms.double( 0.05 ),    	# in cm
     Ptmin = cms.double( 2. ),
     nStubsmin = cms.int32( 4 )
)

process.p = cms.Path( process.L1TrackEtMiss )

process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example_withEtMiss.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)

process.Out.outputCommands.append( 'keep *_*_*_VTX' )
process.Out.outputCommands.append( 'keep *_*_*_ETMISS' )
process.Out.outputCommands.append('keep *_generator_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




