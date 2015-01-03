import FWCore.ParameterSet.Config as cms

process = cms.Process("digiTest")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siPixelRawData'),
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
)
process.source = cms.Source("PoolSource",
    fileNames =  cms.untracked.vstring(
       'file:step2.root'
       )
)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023Muon_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('./DigiTest_pix.root')
)
process.analysis = cms.EDAnalyzer("DigiValidation",
    Verbosity = cms.untracked.bool(False),
#    src = cms.InputTag("SiPixelDigis"),
    src = cms.InputTag("simSiPixelDigis", "Pixel"),
    simG4 = cms.InputTag("g4SimHits"),
    PhiMin = cms.double(10.0),
    PhiMax = cms.double(12.0)                                  
)
#process.digi_step = cms.Sequence(process.siPixelRawData*process.siPixelDigis)
process.p = cms.Path(process.analysis)

# customisation of the process.                                                                                                                              

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms                                                 
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muondev

#call to customisation function cust_2023Muondev imported from SLHCUpgradeSimulations.Configuration.combinedCustoms                                          
process = cust_2023Muondev(process)
