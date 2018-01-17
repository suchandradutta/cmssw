#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/GeometryObjects/interface/PTelescopeParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/TrackerPhase2TestBeam/interface/TelescopeParametersFromDD.h"

class PTelescopeParametersDBBuilder : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
public:
  
  PTelescopeParametersDBBuilder( const edm::ParameterSet& ) {}
  
  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override {}
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {}
};

void
PTelescopeParametersDBBuilder::beginRun( const edm::Run&, edm::EventSetup const& es ) 
{
  PTelescopeParameters* ptp = new PTelescopeParameters;
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  if( !mydbservice.isAvailable())
  {
    edm::LogError( "PTelescopeParametersDBBuilder" ) << "PoolDBOutputService unavailable";
    return;
  }
  edm::ESTransientHandle<DDCompactView> cpv;
  es.get<IdealGeometryRecord>().get( cpv );

  TelescopeParametersFromDD builder;
  builder.build( &(*cpv), *ptp );
  
  if( mydbservice->isNewTagRequest( "PTelescopeParametersRcd" ))
  {
    mydbservice->createNewIOV<PTelescopeParameters>( ptp, mydbservice->beginOfTime(), mydbservice->endOfTime(), "PTelescopeParametersRcd" );
  } else
  {
    edm::LogError( "PTelescopeParametersDBBuilder" ) << "PTelescopeParameters and PTelescopeParametersRcd Tag already present";
  }
}

DEFINE_FWK_MODULE(PTelescopeParametersDBBuilder);
