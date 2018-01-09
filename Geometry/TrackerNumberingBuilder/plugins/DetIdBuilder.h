#ifndef Geometry_TrackerNumberingBuilder_DetIdBuilder_H
# define Geometry_TrackerNumberingBuilder_DetIdBuilder_H

# include "Geometry/TrackerNumberingBuilder/interface/CmsTrackerStringToEnum.h"
# include "FWCore/ParameterSet/interface/types.h"
# include <ostream>
#include <vector>
#include <array>

class GeometricDet;

/**
 * Class to build a geographicalId.
 */

class DetIdBuilder {
 public:
  DetIdBuilder(std::vector<int> detidShifts);
  void build(GeometricDet* telescope);

 private:
  void iterate(uint32_t parentId, GeometricDet* volume, int siblingCounter, int hierarchyLevel);

  std::vector<int> detIdShifts_;
  int numHierarchyLevels_;

  CmsTrackerStringToEnum theCmsTrackerStringToEnum_;
};

#endif
