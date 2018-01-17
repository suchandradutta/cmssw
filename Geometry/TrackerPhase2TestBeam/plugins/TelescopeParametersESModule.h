#ifndef Geometry_TrackerPhase2TestBeam_TelescopeParametersESModule_H
#define Geometry_TrackerPhase2TestBeam_TelescopeParametersESModule_H

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include <memory>

namespace edm {
  class ConfigurationDescriptions;
}
class PTelescopeParameters;
class PTelescopeParametersRcd;

class  TelescopeParametersESModule: public edm::ESProducer
{
 public:
  TelescopeParametersESModule( const edm::ParameterSet & );
  ~TelescopeParametersESModule( void ) override;

  typedef std::shared_ptr<PTelescopeParameters> ReturnType;

  static void fillDescriptions( edm::ConfigurationDescriptions & );
  
  ReturnType produce( const PTelescopeParametersRcd & );
};

#endif
