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

/*void
DUTBuilder::sortNS( DDFilteredView& fv, GeometricDet* parent )
{  
  GeometricDet::ConstGeometricDetContainer & children = parent->components();
  std::stable_sort( children.begin(), children.end(), LessZ());
  
  for(auto& child : uint32_t i = 0; i < children.size(); i++ )
  {
    uint32_t temp= children[i]->type();
    det->component(i)->setGeographicalID(temp%100);  // it relies on the fact that the GeometricDet::GDEnumType enumerators used to identify the subdetectors in the upgrade geometries are equal to the ones of the present detector + n*100
  }
  }*/




