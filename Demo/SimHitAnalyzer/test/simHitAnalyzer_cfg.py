# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_cfi --conditions auto:upgradePLS3 -n 2 --eventcontent FEVTDEBUG -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023Muon --geometry Extended2023Muon,Extended2023MuonReco --magField 38T_PostLS1 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('ANA')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2023tiltedReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023tilted_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.trackerGeometry.applyAlignment = cms.bool(False)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring ([
        "file:/afs/cern.ch/work/g/ghugo/private/MinBias_TuneZ2star_14TeV_pythia6_Phase2_cff_GEN_SIM_20.root"
        #"file:/afs/cern.ch/user/g/ghugo/CMSSW_9_3_0_pre2/src/MinBias_TuneZ2star_14TeV_pythia6.root"
        #"file:/afs/cern.ch/user/g/ghugo/CMSSW_9_3_0_pre2/src/20000.0_FourMuPt1_200+FourMuPt_1_200_pythia8_2023D17_GenSimHLBeamSpotFull+DigiFullTrigger_2023D17+RecoFullGlobal_2023D17+HARVESTFullGlobal_2023D17/step1.root" 
        ])
)

process.options = cms.untracked.PSet(

)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.TFileService = cms.Service ("TFileService",
    fileName = cms.string ("hsimhit_minbias_Phase2_20.root")
)

process.simHitAnalyzer = cms.EDAnalyzer ("SimHitAnalyzer",
    simHitsBarrelHighTof  =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelBarrelHighTof"),
    simHitsBarrelLowTof   =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelBarrelLowTof"),
    simHitsEndcapHighTof  =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelEndcapHighTof"),
    simHitsEndcapLowTof   =  cms.InputTag  ("g4SimHits",  "TrackerHitsPixelEndcapLowTof"),
)

process.myPath = cms.Path (process.simHitAnalyzer)
