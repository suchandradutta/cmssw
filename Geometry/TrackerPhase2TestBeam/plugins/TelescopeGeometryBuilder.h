#ifndef Geometry_TrackerPhase2TestBeam_TelescopeGeometryBuilder_H
# define Geometry_TrackerPhase2TestBeam_TelescopeGeometryBuilder_H

#include "Geometry/TrackerNumberingBuilder/interface/CmsTrackerStringToEnum.h"
#include "FWCore/ParameterSet/interface/types.h"
#include <string>
#include <vector>

class GeometricDet;
class DDCompactView;

/**
 * High level class to build a tracker. It will only build subdets,
 * then call subdet builders
 */

class TelescopeGeometryBuilder
{
public:
  TelescopeGeometryBuilder( void );
  const GeometricDet* construct( const DDCompactView* cpv, std::vector<int> detidShifts);
  
protected:

  std::string attribute;  
  CmsTrackerStringToEnum theCmsTrackerStringToEnum;
};

#endif
