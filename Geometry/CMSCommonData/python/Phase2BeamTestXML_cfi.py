import FWCore.ParameterSet.Config as cms

# This config was generated automatically using generate2023Geometry.py
# If you notice a mistake, please update the generating script, not just this config

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring(
        'Geometry/CMSCommonData/data/materials.xml',
        'Geometry/CMSCommonData/data/rotations.xml',
        'Geometry/CMSCommonData/data/extend/cmsextent.xml',
        'Geometry/CMSCommonData/data/cms.xml',
        'Geometry/CMSCommonData/data/cmsMother.xml',
        
        'Geometry/CMSCommonData/data/Phase2BeamTestConstants.xml',
        'Geometry/CMSCommonData/data/DUT.xml',
        'Geometry/CMSCommonData/data/telescope.xml'
      
    ),
    rootNodeName = cms.string('cms:OCMS')
)
