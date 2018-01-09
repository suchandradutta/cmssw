#include "Geometry/TrackerNumberingBuilder/plugins/ActiveSensorBuilder.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <bitset>

ActiveSensorBuilder::ActiveSensorBuilder() {}

void ActiveSensorBuilder::buildComponent( DDFilteredView& fv, GeometricDet* parent, std::string attribute ) {

  GeometricDet* activeSensor  = new GeometricDet(&fv, theCmsTrackerStringToEnum.type(ExtractStringFromDDD::getString(attribute,&fv)));
  static const std::string isLower = "DUTInnerSensor";
  static const std::string isUpper = "DUTOuterSensor";
  static const std::string isPixel = "Phase1PixelSensor";

  uint32_t temp = 0;
  if (ExtractStringFromDDD::getString(isLower, &fv) == "true") {
    temp = 1;   
  } 

  else if (ExtractStringFromDDD::getString(isUpper,&fv) == "true") {
    temp = 2;
  } 

  else if (ExtractStringFromDDD::getString(isPixel,&fv) == "true") {
    temp = 0;
  }
  
  else {
    edm::LogError("ActiveSensorBuilder") << "Child is nor DUTInnerSensor nor DUTOuterSensor nor Phase1PixelSensor.";
  }

  //uint32_t id = (parent->geographicalID().rawId() << 2) | temp;    
  uint32_t id = temp;
  activeSensor->setGeographicalID(DetId(id));

  parent->addComponent(activeSensor);
}



