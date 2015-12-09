import FWCore.ParameterSet.Config as cms

siStripBadModuleFedErrService = cms.Service("SiStripBadModuleFedErrService",
                                       appendToDataLabel = cms.string(''),
                                       FileName = cms.string('DQM.root'),
                                       BadStripCutoff = cms.double(0.8)
                                       )




