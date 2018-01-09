#include "Geometry/TrackerNumberingBuilder/plugins/DUTBuilder.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TrackerNumberingBuilder/plugins/ActiveSensorBuilder.h"

#include <bitset>

DUTBuilder::DUTBuilder() {}

void DUTBuilder::buildComponent( DDFilteredView& fv, GeometricDet* myDUTHolder, std::string attribute ) {
  ActiveSensorBuilder myActiveSensorBuilder;

  GeometricDet* dut = new GeometricDet( &fv, theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv )));
  switch( theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv ))) {
  case GeometricDet::DUT:
    myActiveSensorBuilder.build( fv, dut, attribute);      
    break;
  default:
    edm::LogError( "DUTBuilder" ) << " ERROR - Could not find a DUT, but found a " << ExtractStringFromDDD::getString( attribute, &fv );
  }
  
  myDUTHolder->addComponent( dut );
}


void DUTBuilder::sortNS( DDFilteredView& fv, GeometricDet* parent ) {  
  GeometricDet::ConstGeometricDetContainer& myDUTs = parent->components();
  std::stable_sort( myDUTs.begin(), myDUTs.end(), LessZ());
  
  for (uint32_t counter = 1; counter <= myDUTs.size(); counter++) {
    //uint32_t id = (parent->geographicalID().rawId() << 5) | counter;
    uint32_t id = counter;
    parent->component(counter-1)->setGeographicalID(DetId(id));
  }
}




