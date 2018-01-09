#include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerLevelBuilder.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


void CmsTrackerLevelBuilder::build (
				    DDFilteredView& fv, 
				    GeometricDet* parent,
				    std::string attribute){

  LogTrace("GeometricDetBuilding") << std::string(3*fv.history().size(),'-') 
				   << "+ "
				   << ExtractStringFromDDD::getString(attribute,&fv) << " " 
				   << parent->type() << " " 
				   << parent->name() 
				   << std::endl;

  bool doLayers = fv.firstChild(); // descend to the first Layer  

  while (doLayers) {
    buildComponent(fv,parent,attribute);      
    doLayers = fv.nextSibling(); // go to next layer
  }

  fv.parent();

  sortNS(fv,parent);


  // TEST
  /* GeometricDet::ConstGeometricDetContainer& children = parent->components();
  
  for (const GeometricDet* child : children) {
    std::cout << "GeometricDet = " << _CmsTrackerStringToEnum.name(child->type())
	      << ", DetId = " << child->geographicalID().rawId() 
	      << ", x = " << child->translation().X() 
	      << ", y = " << child->translation().Y()
	      << ", z = " << child->translation().Z()
	      << ", phi = "  << child->phi() * 180. / M_PI << std::endl;
	      }*/
  // END TEST


}
