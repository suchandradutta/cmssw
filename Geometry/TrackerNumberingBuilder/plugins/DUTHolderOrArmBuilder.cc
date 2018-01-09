#include "Geometry/TrackerNumberingBuilder/plugins/DUTHolderOrArmBuilder.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TrackerNumberingBuilder/plugins/DUTBuilder.h"
#include "Geometry/TrackerNumberingBuilder/plugins/PlaneBuilder.h"

#include <bitset>

DUTHolderOrArmBuilder::DUTHolderOrArmBuilder() {}

void DUTHolderOrArmBuilder::buildComponent( DDFilteredView& fv, GeometricDet* telescope, std::string attribute ) {

  DUTBuilder myDUTBuilder; // TO DO: why not having the build directly at construction time?
  PlaneBuilder myPlaneBuilder;
  GeometricDet* dutHolderOrArm = new GeometricDet( &fv, theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv )));
  switch( theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv ))) {
    // DUT holder
  case GeometricDet::DUTHolder:
    myDUTBuilder.build(fv, dutHolderOrArm, attribute);      
    break;
    // Telescope arm
  case GeometricDet::Arm:
    // TEST        
    /*std::cout << "arm DetId = " << dutHolderOrArm->geographicalID().rawId() 
	      << ", x = " << dutHolderOrArm->translation().X() 
	      << ", y = " << dutHolderOrArm->translation().Y()
	      << ", z = " << dutHolderOrArm->translation().Z()
	      << ", phi = "  << dutHolderOrArm->phi() * 180. / M_PI << std::endl;*/
    // END TEST
    myPlaneBuilder.build( fv, dutHolderOrArm, attribute);      
    break;
  default:
    edm::LogError( "DUTHolderOrArmBuilder" ) << " ERROR - I was expecting a DUTHolder or an Arm, I got a " << ExtractStringFromDDD::getString( attribute, &fv );
  }
  
  telescope->addComponent(dutHolderOrArm);
}


void DUTHolderOrArmBuilder::sortNS( DDFilteredView& fv, GeometricDet* parent ) {  
  GeometricDet::ConstGeometricDetContainer& myDutHolderOrArms = parent->components();
  std::stable_sort( myDutHolderOrArms.begin(), myDutHolderOrArms.end(), LessZ());
  
  for (uint32_t counter = 1; counter <= myDutHolderOrArms.size(); counter++) {
    //uint32_t id = (parent->geographicalID().rawId() << 2) | counter;
    uint32_t id = counter;
    parent->component(counter-1)->setGeographicalID(DetId(id));
  }
}





