#include "Geometry/TrackerNumberingBuilder/plugins/TrackerGeometricDetESModule.h"
#include "Geometry/TrackerNumberingBuilder/plugins/TelescopeGeometryBuilder.h"
//#include "Geometry/TrackerNumberingBuilder/plugins/CondDBCmsTrackerConstruction.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDVectorGetter.h"
#include "DetectorDescription/Core/interface/DDutils.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <memory>

using namespace edm;

TrackerGeometricDetESModule::TrackerGeometricDetESModule( const edm::ParameterSet & p ) 
  : fromDDD_( p.getParameter<bool>( "fromDDD" ))
{
  setWhatProduced( this );
}

TrackerGeometricDetESModule::~TrackerGeometricDetESModule( void ) {}


void TrackerGeometricDetESModule::fillDescriptions( edm::ConfigurationDescriptions & descriptions ) {
 
  edm::ParameterSetDescription desc;
  desc.add<bool>( "fromDDD", true );
  descriptions.add( "trackerNumberingGeometry", desc );
}


std::unique_ptr<GeometricDet> TrackerGeometricDetESModule::produce( const IdealGeometryRecord & iRecord ) { 

  std::cout << "TrackerGeometricDetESModule::produce  fromDDD_ = " << fromDDD_ << std::endl;

  edm::ESTransientHandle<DDCompactView> cpv;
  iRecord.get( cpv );

  TelescopeGeometryBuilder theTelescopeGeometryBuilder;
  return std::unique_ptr<GeometricDet> (const_cast<GeometricDet*>( theTelescopeGeometryBuilder.construct(&(*cpv), dbl_to_int( DDVectorGetter::get( "telescopeDetIdShifts" )))));

}

DEFINE_FWK_EVENTSETUP_MODULE( TrackerGeometricDetESModule );
