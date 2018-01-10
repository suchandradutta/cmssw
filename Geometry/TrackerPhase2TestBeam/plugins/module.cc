//<<<<<< INCLUDES                                                       >>>>>>
#include "DetectorDescription/Core/interface/DDAlgorithmFactory.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Geometry/TrackerPhase2TestBeam/plugins/DDTelescopePlanesAlgo.h"


DEFINE_EDM_PLUGIN (DDAlgorithmFactory, DDTelescopePlanesAlgo, "phase2TestBeam:DDTelescopePlanesAlgo");


#include "Geometry/TrackerNumberingBuilder/interface/CmsTrackerStringToEnum.h"
# include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerLevelBuilder.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
