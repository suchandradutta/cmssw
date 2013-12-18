#ifndef L1TkElectronTrackMatchAlgo_HH
#define L1TkElectronTrackMatchAlgo_HH

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

namespace L1TkElectronTrackMatchAlgo {
    typedef L1TkTrack_PixelDigi_   L1TkTrackType;
    void doMatch(const edm::Ref< l1extra::L1EmParticleCollection >& egRef, const edm::Ptr< L1TkTrackType >& trkPtr, float& dph, float& dr, float& deta);
   
   float deltaR(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr);
   float deltaPhi(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr);
   float deltaEta(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr);
   GlobalPoint calorimeterPosition(float phi, float eta, float e);
}  
#endif
