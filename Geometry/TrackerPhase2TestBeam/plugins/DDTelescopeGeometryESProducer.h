#ifndef Geometry_TrackerPhase2TestBeam_DDTelescopeGeometryESProducer_H
#define Geometry_TrackerPhase2TestBeam_DDTelescopeGeometryESProducer_H

#include "FWCore/Framework/interface/ESProducer.h"

namespace edm {
  class ConfigurationDescriptions;
  class ParameterSet;
}

class GeometricDet;
class IdealGeometryRecord;

class DDTelescopeGeometryESProducer : public edm::ESProducer
{
public:
  DDTelescopeGeometryESProducer( const edm::ParameterSet & p );
  ~DDTelescopeGeometryESProducer( void ) override; 
  std::unique_ptr<GeometricDet> produce( const IdealGeometryRecord & );

  static void fillDescriptions( edm::ConfigurationDescriptions & descriptions );
  
private:
  bool fromDDD_;
};

#endif





