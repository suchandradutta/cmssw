#include "SLHCUpgradeSimulations/TTAnalysis/interface/AnalObjects.h"
#include <vector>

namespace {
  struct dictionary {
    TTStudy::Event          rv1;
    TTStudy::SimTrack       rv2;
    TTStudy::Track          rv3;
    TTStudy::Electron       rv4;
    TTStudy::Tracklet       rv5;
    TTStudy::GenParticle    rv6;
    TTStudy::L1Jet          rv7;
    TTStudy::L1Muon         rv8;
    std::vector<TTStudy::SimTrack>       vrv1;
    std::vector<TTStudy::Track>          vrv2;
    std::vector<TTStudy::Electron>       vrv3;
    std::vector<TTStudy::Tracklet>       vrv4;
    std::vector<TTStudy::GenParticle>    vrv5;
    std::vector<TTStudy::L1Jet>          vrv6;
    std::vector<TTStudy::L1Muon>         vrv7;
  };
}
