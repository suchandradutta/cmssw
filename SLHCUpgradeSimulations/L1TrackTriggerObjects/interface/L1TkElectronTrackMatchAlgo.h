#ifndef L1TkElectronTrackMatchAlgo_HH
#define L1TkElectronTrackMatchAlgo_HH

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

namespace L1TkElectronTrackMatchAlgo {
  typedef L1TkTrack_PixelDigi_   L1TkTrackType;
  typedef std::vector< L1TkTrackType >       L1TkTrackCollectionType;
  void doMatch(l1extra::L1EmParticleCollection::const_iterator egIter, L1TkTrackCollectionType::const_iterator trkIter, double& dph, float&  dr, float& deta);

  float deltaR(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter);
  double deltaPhi(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter);
  float deltaEta(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter);
  GlobalPoint calorimeterPosition(float phi, float eta, float e);

}  
#endif
