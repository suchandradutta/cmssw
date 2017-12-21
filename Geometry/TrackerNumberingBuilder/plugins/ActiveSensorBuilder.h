#ifndef Geometry_TrackerNumberingBuilder_ActiveSensorBuilder_H
# define Geometry_TrackerNumberingBuilder_ActiveSensorBuilder_H

# include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerLevelBuilder.h"
# include "FWCore/ParameterSet/interface/types.h"
# include <string>

/**
 * 
 */
class ActiveSensorBuilder : public CmsTrackerLevelBuilder {
public:
  ActiveSensorBuilder();

private:
  void buildComponent( DDFilteredView& , GeometricDet*, std::string ) override;
};

#endif
