############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")
 

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

<<<<<<< HEAD
=======

>>>>>>> my_dev
############################################################
# input source
############################################################

<<<<<<< HEAD
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

Source_Files = cms.untracked.vstring('root://eoscms//store/mc/UpgFall13d/SingleMuMinusFlatPt0p2To100/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D217598A-0539-E311-8A7E-00261894388D.root')
=======
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

Source_Files = cms.untracked.vstring('root://eoscms//store/group/comm_trigger/L1TrackTrigger/BE5D_612_SLHC6_patch1/singleEle/SingleEle_NoPU.root') 
>>>>>>> my_dev
process.source = cms.Source("PoolSource", fileNames = Source_Files)


############################################################
# track trigger
############################################################

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')


############################################################
# output definition
############################################################

<<<<<<< HEAD
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleMuMinus_BE5D_TrkPerf.root'), closeFileFast = cms.untracked.bool(True))
=======
process.TFileService = cms.Service("TFileService", fileName = cms.string('SingleEl_noPU_BE5D_TrkPerf.root'), closeFileFast = cms.untracked.bool(True))
>>>>>>> my_dev


############################################################
# other statements
############################################################

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# Path definitions & schedule
############################################################

<<<<<<< HEAD
process.BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.BeamSpotFromSim*process.L1Tracks)

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker')
process.ana = cms.Path(process.L1TrackNtuple)

process.schedule = cms.Schedule(process.TT_step,process.ana)
=======
# Remake stubs (running with zMatching=False)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)

process.BeamSpotFromSim = cms.EDProducer("BeamSpotFromSimProducer")
process.TT_step = cms.Path(process.BeamSpotFromSim*process.L1Tracks)


# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(11)
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.TT_step,process.ana)
>>>>>>> my_dev

