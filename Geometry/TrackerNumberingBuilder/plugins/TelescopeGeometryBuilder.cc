#include "Geometry/TrackerNumberingBuilder/plugins/TelescopeGeometryBuilder.h"

#include <utility>
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerNumberingBuilder/plugins/ExtractStringFromDDD.h"
#include "Geometry/TrackerNumberingBuilder/plugins/DetIdBuilder.h"

#include "Geometry/TrackerNumberingBuilder/plugins/DUTHolderOrArmBuilder.h"

using namespace cms;

TelescopeGeometryBuilder::TelescopeGeometryBuilder( void )
{}

const GeometricDet* TelescopeGeometryBuilder::construct( const DDCompactView* cpv, std::vector<int> detidShifts) {
  attribute = "TelescopeDDDStructure";

  DDSpecificsHasNamedValueFilter filter{ attribute }; 
  DDFilteredView fv( *cpv, filter ); 

  // TELESCOPE VOLUME
  fv.firstChild();  // TO DO: Add check that child exist!!


  DUTHolderOrArmBuilder myDUTHolderOrArmBuilder;
  GeometricDet* telescope = new GeometricDet( &fv, theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv )));
  switch( theCmsTrackerStringToEnum.type( ExtractStringFromDDD::getString( attribute, &fv ))) {
  case GeometricDet::Telescope:
     myDUTHolderOrArmBuilder.build(fv, telescope, attribute);      
    break;
  default:
    edm::LogError( "TelescopeGeometryBuilder" ) << " ERROR - I was expecting a Telescope, I got a " << ExtractStringFromDDD::getString( attribute, &fv );
  }


  
  DetIdBuilder myDetIdBuilder( std::move(detidShifts) );
  myDetIdBuilder.build(telescope);

  fv.parent(); // come back to world volume
 
  return telescope;
}

