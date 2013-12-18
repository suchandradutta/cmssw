// -*- C++ -*-
//
// Package:    L1TrackTriggerObjects
// Class:      L1TkElectronTrackMatchAlgo
// 
/**\class L1TkElectronTrackMatchAlgo 

 Description: Algorithm to match L1EGamma oject with L1Track candidates

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  S. Dutta and A. Modak
//         Created:  Wed Dec 4 12 11:55:35 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronTrackMatchAlgo.h"
namespace L1TkElectronTrackMatchAlgo {
  // ------------ match EGamma and Track
  void doMatch(const edm::Ref< l1extra::L1EmParticleCollection >& egRef, const edm::Ptr< L1TkTrackType >& trkPtr, float& dph, float& dr, float& deta) {
    GlobalPoint egPos = L1TkElectronTrackMatchAlgo::calorimeterPosition(egRef->phi(), egRef->eta(), egRef->energy());
    dph  = L1TkElectronTrackMatchAlgo::deltaPhi(egPos, trkPtr);
    dr   = L1TkElectronTrackMatchAlgo::deltaR(egPos, trkPtr);
    deta = L1TkElectronTrackMatchAlgo::deltaEta(egPos, trkPtr);
  }
  // --------------- calculate deltaR between Track and EGamma object
  float deltaPhi(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr){
  
    float corr = (asin(epos.perp()*trkPtr->getRInv()/(2.0)));
    float trk_phi = trkPtr->getMomentum().phi(); 
    float phi_diff = fabs(reco::deltaPhi(epos.phi(),trk_phi));
    
    double dif1 = phi_diff - corr;
    double dif2 = phi_diff + corr; 
    if (abs(dif1) < abs(dif2)) return dif1;
    else return dif2; 
  }
// --------------- calculate deltaPhi between Track and EGamma object                 
  float deltaR(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr){
    float dPhi = fabs(reco::deltaPhi(epos.phi(), trkPtr->getMomentum().phi()));
    float dEta = (epos.eta() - trkPtr->getMomentum().eta());
    return sqrt(dPhi*dPhi + dEta*dEta);
  }
  // --------------- calculate deltaEta between Track and EGamma object                 
  float deltaEta(GlobalPoint epos, const edm::Ptr< L1TkTrackType >& trkPtr){
    float corr_eta = 999.0;
    float er = epos.perp();
    float ez = epos.z();
    float z0 = trkPtr->getMomentum().z();
    float theta = 0.0;
    if (ez >= 0) theta = atan(er/fabs(ez-z0));
    else theta = M_PI - atan(er/fabs(ez-z0));
    corr_eta = -1.0 * log(tan(theta/2.0));
    float deleta = (corr_eta - trkPtr->getMomentum().eta());
    return deleta;
  }
  // -------------- get Calorimeter position
  GlobalPoint calorimeterPosition(float phi, float eta, float e) {
    float x = 0; 
    float y = 0;
    float z = 0;
    float depth = 0.89*(7.7+ log(e) );
    float theta = 2*atan(exp(-1*eta));
    float r = 0;
    if( fabs(eta) > 1.479 ) 
      { 
	float ecalZ = 315.4*fabs(eta)/eta;
	
	r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
	x = r * cos( phi ) * sin( theta );
	y = r * sin( phi ) * sin( theta );
	z = r * cos( theta );
      }
    else
      {
	float rperp = 129.0;
	float zface =  sqrt( cos( theta ) * cos( theta ) /
			     ( 1 - cos( theta ) * cos( theta ) ) * 
			     rperp * rperp ) * abs( eta ) / eta;  
	r = sqrt( rperp * rperp + zface * zface ) + depth;
	x = r * cos( phi ) * sin( theta );
	y = r * sin( phi ) * sin( theta );
	z = r * cos( theta );
      }
    GlobalPoint pos(x,y,z);
    return pos;
  }
}
