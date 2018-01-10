#ifndef Geometry_TrackerPhase2TestBeam_DUTHolderOrArmBuilder_H
# define Geometry_TrackerPhase2TestBeam_DUTHolderOrArmBuilder_H

# include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerLevelBuilder.h"
# include "FWCore/ParameterSet/interface/types.h"
# include <string>

/**
 * Abstract Class to construct a Level in the hierarchy
 */
class DUTHolderOrArmBuilder : public CmsTrackerLevelBuilder
{
public:
  DUTHolderOrArmBuilder();

private:
  void buildComponent( DDFilteredView& , GeometricDet*, std::string ) override;
  void sortNS( DDFilteredView& , GeometricDet* ) override;
};

#endif
