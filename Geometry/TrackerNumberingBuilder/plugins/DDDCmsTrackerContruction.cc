#include "Geometry/TrackerNumberingBuilder/plugins/DDDCmsTrackerContruction.h"

#include <utility>
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerDetIdBuilder.h"

#include "Geometry/TrackerNumberingBuilder/plugins/DUTBuilder.h"
#include "Geometry/TrackerNumberingBuilder/plugins/PlaneBuilder.h"

using namespace cms;

DDDCmsTrackerContruction::DDDCmsTrackerContruction( void )
{}

const GeometricDet*
DDDCmsTrackerContruction::construct( const DDCompactView* cpv)
{
  attribute = "TelescopeDDDStructure";

  DDSpecificsHasNamedValueFilter filter{ attribute }; 
  DDFilteredView fv( *cpv, filter ); 


  // TELESCOPE VOLUME
  fv.firstChild();
  if (theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString(attribute,&fv)) == GeometricDet::Telescope ) {  
    std::cout << "DDDCmsTrackerContruction::construct I have found the telescope at level 1 " << std::endl;
  }
  else { std::cout << "DDDCmsTrackerContruction::construct I have not found the telescope. " << std::endl; }
  
  GeometricDet* telescope = new GeometricDet( &fv, GeometricDet::Telescope );




  // DUT HOLDER AND ARMS
  bool doLayers = fv.firstChild();

  while (doLayers) {
    //buildComponent(fv,telescope,attribute);      

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
      std::cout << "arm DetId = " << dutHolderOrArm->geographicalID().rawId() 
		<< ", x = " << dutHolderOrArm->translation().X() 
		<< ", y = " << dutHolderOrArm->translation().Y()
		<< ", z = " << dutHolderOrArm->translation().Z()
		<< ", phi = "  << dutHolderOrArm->phi() * 180. / M_PI << std::endl;
      // END TEST
      myPlaneBuilder.build( fv, dutHolderOrArm, attribute);      
      break;
    default:
      edm::LogError( "DDDCmsTrackerContruction" ) << " ERROR - I was expecting a DUTHolder or an Arm, I got a " << ExtractStringFromDDD::getString( attribute, &fv );
    }
  
    telescope->addComponent(dutHolderOrArm);

    doLayers = fv.nextSibling();
  }

  fv.parent(); // come back to telescope volume

  //sortNS(fv,telescope);







  //CmsTrackerBuilder theCmsTrackerBuilder;
  //theCmsTrackerBuilder.build( fv, telescope, attribute );
  
  //CmsTrackerDetIdBuilder theCmsTrackerDetIdBuilder( std::move(detidShifts) );
  //tracker = theCmsTrackerDetIdBuilder.buildId( tracker );
  telescope->setGeographicalID(DetId(500));

  fv.parent(); // come back to world volume
 
  return telescope;
}

