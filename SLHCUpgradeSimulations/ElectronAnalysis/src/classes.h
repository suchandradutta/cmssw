#include "SLHCUpgradeSimulations/ElectronAnalysis/interface/AnalObjects.h"
#include <vector>

namespace {
  struct dictionary {
    TTIStudy::Event          rv1;
    TTIStudy::SimTrack       rv2;
    TTIStudy::Track          rv3;
    TTIStudy::Electron       rv4;
    TTIStudy::Tracklet       rv5;
    std::vector<TTIStudy::SimTrack>       vrv1;
    std::vector<TTIStudy::Track>          vrv2;
    std::vector<TTIStudy::Electron>       vrv3;
    std::vector<TTIStudy::Tracklet>       vrv4;
  };
}
