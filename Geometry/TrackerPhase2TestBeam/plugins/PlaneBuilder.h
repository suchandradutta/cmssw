#ifndef Geometry_TrackerPhase2TestBeam_PlaneBuilder_H
# define Geometry_TrackerPhase2TestBeam_PlaneBuilder_H

# include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerLevelBuilder.h"
# include "FWCore/ParameterSet/interface/types.h"
# include <string>

/**
 * Abstract Class to construct a Level in the hierarchy
 */
class PlaneBuilder : public CmsTrackerLevelBuilder
{
public:
  PlaneBuilder();

private:
  void buildComponent( DDFilteredView& , GeometricDet*, std::string ) override;
  void sortNS( DDFilteredView& , GeometricDet* ) override;
};

#endif
