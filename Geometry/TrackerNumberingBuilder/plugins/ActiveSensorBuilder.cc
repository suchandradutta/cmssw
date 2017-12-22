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
  //static const std::string isPixel = "DUTOuterSensor";

  if (ExtractStringFromDDD::getString(isLower, &fv) == "true"){
    uint32_t temp = 1;
    activeSensor->setGeographicalID(DetId(temp));
  } else if (ExtractStringFromDDD::getString(isUpper,&fv) == "true"){
    uint32_t temp = 2;
    activeSensor->setGeographicalID(DetId(temp));
  } else {
    edm::LogError("ActiveSensorBuilder") << "DUT Sensors are not inner nor outer.";
  }
  parent->addComponent(activeSensor);
}



