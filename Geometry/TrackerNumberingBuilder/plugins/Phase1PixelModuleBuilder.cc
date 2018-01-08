#include "Geometry/TrackerNumberingBuilder/plugins/Phase1PixelModuleBuilder.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TrackerNumberingBuilder/plugins/ActiveSensorBuilder.h"

#include <bitset>

Phase1PixelModuleBuilder::Phase1PixelModuleBuilder() {}

void Phase1PixelModuleBuilder::buildComponent( DDFilteredView& fv, GeometricDet* plane, std::string attribute ) {
  ActiveSensorBuilder myActiveSensorBuilder;

  GeometricDet* phase1PixelModule = new GeometricDet( &fv, theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv )));
  switch( theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv ))) {
  case GeometricDet::Phase1PixelModule:
    myActiveSensorBuilder.build( fv, phase1PixelModule, attribute);      
    break;
  default:
    edm::LogError( "Phase1PixelModuleBuilder" ) << " ERROR - Could not find a Phase1PixelModule, but found a " << ExtractStringFromDDD::getString( attribute, &fv );
  }
  
  plane->addComponent(phase1PixelModule);
}

/*void
Phase1PixelModuleBuilder::sortNS( DDFilteredView& fv, GeometricDet* parent )
{  
  GeometricDet::ConstGeometricDetContainer & children = parent->components();
  std::stable_sort( children.begin(), children.end(), LessZ());
  
  for(auto& child : uint32_t i = 0; i < children.size(); i++ )
  {
    uint32_t temp= children[i]->type();
    det->component(i)->setGeographicalID(temp%100);  // it relies on the fact that the GeometricDet::GDEnumType enumerators used to identify the subdetectors in the upgrade geometries are equal to the ones of the present detector + n*100
  }
  }*/




