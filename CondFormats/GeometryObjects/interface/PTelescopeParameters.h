#ifndef CondFormats_GeometryObjects_PTelescopeParameters_h
#define CondFormats_GeometryObjects_PTelescopeParameters_h

#include "CondFormats/Serialization/interface/Serializable.h"

class PTelescopeParameters {
 public:
  PTelescopeParameters( void ) { } 
  ~PTelescopeParameters( void ) { }

  std::vector<int> detIdShifts;
  std::vector<int> detIdMasks;
  
  COND_SERIALIZABLE;
};

#endif
