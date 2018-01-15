# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_cfi --conditions auto:phase2_realistic -n 2 --eventcontent FEVTDEBUG -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D17 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D17_cff')


process.load("Configuration.Geometry.GeometryTrackerPhase2TestBeam_cff")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")     





process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring ([
        #"file:/afs/cern.ch/user/g/ghugo/CMSSW_10_0_0_pre1/src/SingleMuPt1000_pythia8_cfi_py_GEN_SIM.root"
        #"file:/afs/cern.ch/user/g/ghugo/CMSSW_10_0_0_pre1/src/MinBias_TuneZ2star_14TeV_pythia6_Phase2_cff_GEN_SIM_2000.root"
        "file:/afs/cern.ch/user/g/ghugo/CMSSW_10_0_0_pre1/src/TrackerPhase2TestBeam_GEN_SIM.root"    
        ])
)

process.options = cms.untracked.PSet(

)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')  # 'auto:phase2_realistic'

process.TFileService = cms.Service ("TFileService",
    fileName = cms.string ("pixelTelescope_simhitmap.root")
)

process.simHitAnalyzer = cms.EDAnalyzer ("SimHitAnalyzer",
    simHitsBarrelHighTof  =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelBarrelHighTof"),
    simHitsBarrelLowTof   =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelBarrelLowTof"),
    simHitsEndcapHighTof  =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelEndcapHighTof"),
    simHitsEndcapLowTof   =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelEndcapLowTof"),
)

process.myPath = cms.Path (process.simHitAnalyzer)
