#include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerDetIdBuilder.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <bitset>


CmsTrackerDetIdBuilder::CmsTrackerDetIdBuilder( std::vector<int> detidShifts )
  : detIdShifts_(detidShifts), numHierarchyLevels_(detIdShifts_.size())
{}


void CmsTrackerDetIdBuilder::buildDetIds(GeometricDet* telescope) {
  LogDebug("BuildingTelescopeDetIds") << "Starting to build Telescope DetIds";

  int hierarchyLevel = 0;
  uint32_t siblingCounter = 0;
  uint32_t parentId = 0;

  iterate(parentId, telescope, siblingCounter, hierarchyLevel);


  /*
  //DetId t( DetId::Tracker, 0 );
  telescope->setGeographicalID(DetId(1)); // TO DO: Should create a DetId specific to telescope mother volume (see DataFormats/DetId/interface/DetId.h ).
  // Issue is the space allocated for it is 3 bits, and integers from 1 to 7 are already assigned (cannot used 0).
  */
}


void CmsTrackerDetIdBuilder::iterate(uint32_t parentId, GeometricDet* volume, int siblingCounter, int hierarchyLevel) {

  if (hierarchyLevel >= numHierarchyLevels_) {
    edm::LogError( "TelescopeDetIdBuilder" ) << " ERROR - I reached hierarchyLevel " << hierarchyLevel << ". Level should be strictly lower than " << numHierarchyLevels_;
  }

  else {

    uint32_t indicator = 0;
    if (hierarchyLevel == (numHierarchyLevels_ - 1)) { indicator = volume->geographicalID().rawId(); }
    else { indicator = siblingCounter + 1; }

    uint32_t id = parentId | ( indicator << detIdShifts_.at(hierarchyLevel));
    volume->setGeographicalID(DetId(id));



    std::cout << "GeometricDet = " << volume->type()
	      << ", DetId = " << volume->geographicalID().rawId() 
	      << ", x = " << volume->translation().X() 
	      << ", y = " << volume->translation().Y()
	      << ", z = " << volume->translation().Z()
	      << ", phi = "  << volume->phi() * 180. / M_PI << std::endl;



    hierarchyLevel++;
    parentId = id;
    uint32_t numSiblings = volume->components().size();

    for (uint32_t siblingCounter = 0; siblingCounter < numSiblings; siblingCounter++) {
      GeometricDet* child = volume->component(siblingCounter);
      iterate(parentId, child, siblingCounter, hierarchyLevel);
    }
  }

}
