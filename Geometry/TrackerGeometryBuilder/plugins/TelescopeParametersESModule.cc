#include "Geometry/TrackerPhase2TestBeam/plugins/TelescopeParametersESModule.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/PTelescopeParametersRcd.h"
#include "Geometry/TrackerPhase2TestBeam/interface/TelescopeParametersFromDD.h"
#include "CondFormats/GeometryObjects/interface/PTelescopeParameters.h"

TelescopeParametersESModule::TelescopeParametersESModule( const edm::ParameterSet& )
{
  edm::LogInfo("TELESCOPE") << "TelescopeParametersESModule::TelescopeParametersESModule";

  setWhatProduced(this);
}

TelescopeParametersESModule::~TelescopeParametersESModule()
{ 
}

void
TelescopeParametersESModule::fillDescriptions( edm::ConfigurationDescriptions & descriptions ) 
{
  edm::ParameterSetDescription desc;
  descriptions.add( "telescopeParameters", desc );
}

TelescopeParametersESModule::ReturnType
TelescopeParametersESModule::produce( const PTelescopeParametersRcd& iRecord )
{
  edm::LogInfo("TelescopeParametersESModule") <<  "TelescopeParametersESModule::produce(const PTelescopeParametersRcd& iRecord)" << std::endl;
  edm::ESTransientHandle<DDCompactView> cpv;
  iRecord.getRecord<IdealGeometryRecord>().get( cpv );
    
  PTelescopeParameters* ptp = new PTelescopeParameters();
  TelescopeParametersFromDD builder;
  builder.build( &(*cpv), *ptp );

  return ReturnType( ptp ) ;
}

DEFINE_FWK_EVENTSETUP_MODULE( TelescopeParametersESModule);
