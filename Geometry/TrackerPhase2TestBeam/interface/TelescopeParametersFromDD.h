#ifndef Geometry_TrackerPhase2TestBeam_TelescopeParametersFromDD_h
#define Geometry_TrackerPhase2TestBeam_TelescopeParametersFromDD_h

#include <vector>

class DDCompactView;
class PTelescopeParameters;

class TelescopeParametersFromDD {
 public:
  TelescopeParametersFromDD() {}
  virtual ~TelescopeParametersFromDD() {}

  bool build( const DDCompactView*, PTelescopeParameters& );
};

#endif
