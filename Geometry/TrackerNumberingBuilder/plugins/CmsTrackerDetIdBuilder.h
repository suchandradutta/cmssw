#ifndef Geometry_TrackerNumberingBuilder_CmsTrackerDetIdBuilder_H
# define Geometry_TrackerNumberingBuilder_CmsTrackerDetIdBuilder_H

# include "FWCore/ParameterSet/interface/types.h"
# include <ostream>
#include <vector>
#include <array>

class GeometricDet;

/**
 * Class to build a geographicalId.
 */

class CmsTrackerDetIdBuilder {
 public:
  CmsTrackerDetIdBuilder(std::vector<int> detidShifts);
  void buildDetIds(GeometricDet* telescope);

 private:
  void iterate(uint32_t parentId, GeometricDet* volume, int siblingCounter, int hierarchyLevel);

  std::vector<int> detIdShifts_;
  int numHierarchyLevels_;
};

#endif
