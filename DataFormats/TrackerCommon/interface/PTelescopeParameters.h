#ifndef DataFormats_TrackerCommon_PTelescopeParameters_h
#define DataFormats_TrackerCommon_PTelescopeParameters_h

//#include "CondFormats/Serialization/interface/Serializable.h"

class PTelescopeParameters {
 public:
  PTelescopeParameters( void ) { } 
  ~PTelescopeParameters( void ) { }

  std::vector<int> detIdShifts;
  std::vector<int> detIdMasks;

};

#endif
