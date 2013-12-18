// -*- C++ -*-
//
//
// Producer of a L1TkElectronParticle, for the algorithm matching a L1Track to the L1EG object.
// This code is just an example, I simply pick up the L1Track that is closest to the
// L1EG algorithm.
// The proper producer will be provided by Suchandra & Atanu. 
// 

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"

#include "DataFormats/Math/interface/LorentzVector.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include <string>
#include "TMath.h"


using namespace l1extra ;

//
// class declaration
//

class L1TkElectronTrackProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TkElectronTrackProducer(const edm::ParameterSet&);
      ~L1TkElectronTrackProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run&, edm::EventSetup const&);
      //virtual void endRun(edm::Run&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
	edm::InputTag L1EGammaInputTag;
	edm::InputTag L1TrackInputTag;
	std::string label;

} ;


//
// constructors and destructor
//
L1TkElectronTrackProducer::L1TkElectronTrackProducer(const edm::ParameterSet& iConfig)
{

   L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
   label = iConfig.getParameter<std::string>("label");
   

   produces<L1TkElectronParticleCollection>(label);
}

L1TkElectronTrackProducer::~L1TkElectronTrackProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TkElectronTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 std::auto_ptr<L1TkElectronParticleCollection> result(new L1TkElectronParticleCollection);

 edm::Handle<L1EmParticleCollection> EGammaHandle;
 iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
 std::vector<L1EmParticle>::const_iterator egIter ;

 edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
 iEvent.getByLabel(L1TrackInputTag, L1TkTrackHandle);
 L1TkTrackCollectionType::const_iterator trackIter;

 int ieg = 0;
 for (egIter = EGammaHandle->begin();  egIter != EGammaHandle->end(); ++egIter) {

    edm::Ref< L1EmParticleCollection > EGammaRef( EGammaHandle, ieg );
    ieg ++;

    int ibx = egIter -> bx();
    if (ibx != 0) continue;

    float eta = egIter -> eta();
    float phi = egIter -> phi();
    
	// match the L1EG object with a L1Track
	// here dummy : I simply take the closest track
	// and require that DR < 0.5

	float drmin = 999;
	int itr = 0;
	int itrack = -1;
	for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {
		float Eta = trackIter->getMomentum().eta();
		float Phi = trackIter->getMomentum().phi();
		float deta = eta - Eta;
		float dphi = phi - Phi;
		if (dphi < 0) dphi = dphi + 2.*TMath::Pi();
		if (dphi > TMath::Pi()) dphi = 2.*TMath::Pi() - dphi;
		float dR = sqrt( deta*deta + dphi*dphi );
		if (dR < drmin) {
		  drmin = dR;
		  itrack = itr;
		}
		itr ++ ;
	}

	if (drmin < 0.5 ) {	// found a L1Track matched to the L1EG object

	    edm::Ptr< L1TkTrackType > L1TrackPtr( L1TkTrackHandle, itrack) ;
	    
 	    float px = L1TrackPtr -> getMomentum().x();
	    float py = L1TrackPtr -> getMomentum().y();
	    float pz = L1TrackPtr -> getMomentum().z();
	    float e = sqrt( px*px + py*py + pz*pz );	// massless particle
            math::XYZTLorentzVector TrackP4(px,py,pz,e);

	    float trkisol = -999; 	// dummy

 	    L1TkElectronParticle trkEm( TrackP4, 
				 EGammaRef,
				 L1TrackPtr, 
			         trkisol );

	    result -> push_back( trkEm );

	}  // endif drmin < 0.5

 }  // end loop over EGamma objects

 iEvent.put( result, label );

}


// ------------ method called once each job just before starting event loop  ------------
void
L1TkElectronTrackProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1TkElectronTrackProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TkElectronTrackProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
L1TkElectronTrackProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TkElectronTrackProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TkElectronTrackProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkElectronTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkElectronTrackProducer);



