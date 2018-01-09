#include "Geometry/TrackerNumberingBuilder/plugins/DDDCmsTrackerContruction.h"

#include <utility>
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "Geometry/TrackerNumberingBuilder/plugins/CmsTrackerDetIdBuilder.h"

#include "Geometry/TrackerNumberingBuilder/plugins/DUTHolderOrArmBuilder.h"

using namespace cms;

DDDCmsTrackerContruction::DDDCmsTrackerContruction( void )
{}

const GeometricDet* DDDCmsTrackerContruction::construct( const DDCompactView* cpv) {
  attribute = "TelescopeDDDStructure";

  DDSpecificsHasNamedValueFilter filter{ attribute }; 
  DDFilteredView fv( *cpv, filter ); 

  // TELESCOPE VOLUME
  fv.firstChild();


  DUTHolderOrArmBuilder myDUTHolderOrArmBuilder;
  GeometricDet* telescope = new GeometricDet( &fv, theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv )));
  switch( theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv ))) {
  case GeometricDet::Telescope:
     myDUTHolderOrArmBuilder.build(fv, telescope, attribute);      
    break;
  default:
    edm::LogError( "DDDCmsTrackerContruction" ) << " ERROR - I was expecting a Telescope, I got a " << ExtractStringFromDDD::getString( attribute, &fv );
  }


  //CmsTrackerBuilder theCmsTrackerBuilder;
  //theCmsTrackerBuilder.build( fv, telescope, attribute );
  
  //CmsTrackerDetIdBuilder theCmsTrackerDetIdBuilder( std::move(detidShifts) );
  //tracker = theCmsTrackerDetIdBuilder.buildId( tracker );

  telescope->setGeographicalID(DetId(1)); // TO DO: Should create a DetId specific to telescope mother volume (see DataFormats/DetId/interface/DetId.h ).
  // Issue is the space allocated for it is 3 bits, and integers from 1 to 7 are already assigned (cannot used 0).

  fv.parent(); // come back to world volume
 
  return telescope;
}

