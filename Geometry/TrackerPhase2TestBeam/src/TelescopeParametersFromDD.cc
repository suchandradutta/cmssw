#include "Geometry/TrackerPhase2TestBeam/interface/TelescopeParametersFromDD.h"
#include "CondFormats/GeometryObjects/interface/PTelescopeParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDVectorGetter.h"
#include "DetectorDescription/Core/interface/DDutils.h"


bool TelescopeParametersFromDD::build( const DDCompactView* cvp, PTelescopeParameters& ptp) {

  ptp.detIdShifts = dbl_to_int( DDVectorGetter::get( "telescopeDetIdShifts" ));
  ptp.detIdMasks = dbl_to_int( DDVectorGetter::get( "telescopeDetIdMasks" ));

  return true;
}
