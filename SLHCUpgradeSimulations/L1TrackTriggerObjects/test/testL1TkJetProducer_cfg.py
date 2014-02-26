import FWCore.ParameterSet.Config as cms

process = cms.Process("Jet")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:example_w_Tracks_and_vertex.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/m1_TTbar_BE5D.root'
    #'/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/TTbar_BE5D_97.root'
    )
)


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# --- creates l1extra objects for L1EGamma (here for the Run1 trigger !)

        # raw2digi to get the gct digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )
	# run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path( process.l1extraParticles )


# --- Now run the L1TkJetProducer 

process.L1TkJets = cms.EDProducer("L1TkJetProducer",
        L1CentralJetInputTag = cms.InputTag("l1extraParticles","Central"),      # for Run-1 algos
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
	# cuts on the tracks used to determined the zvertex of the jet (examples) :
	#ZMAX = cms.double( 25. ),	# in cm
	#CHI2MAX = cms.double( 100. ),
        #PTMINTRA = cms.double( 2. ),	# in GeV
)
process.pJets = cms.Path( process.L1TkJets )



process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)


process.Out.outputCommands.append('keep *_L1TkJets_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




