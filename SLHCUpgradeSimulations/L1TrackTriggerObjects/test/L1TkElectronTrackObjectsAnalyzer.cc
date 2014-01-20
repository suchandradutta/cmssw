// -*- C++ -*-
//
// Package:    L1TkElectronTrackObjectsAnalyzer
// Class:      L1TkElectronTrackObjectsAnalyzer
// 
/**\class L1TkElectronTrackObjectsAnalyzer L1TkElectronTrackObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TkElectronTrackObjectsAnalyzer/src/L1TkElectronTrackObjectsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"


#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronEtComparator.h"
#include "TH1F.h"


using namespace l1extra;



//
// class declaration
//

class L1TkElectronTrackObjectsAnalyzer : public edm::EDAnalyzer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

  explicit L1TkElectronTrackObjectsAnalyzer(const edm::ParameterSet&);
  ~L1TkElectronTrackObjectsAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void fillIntegralHistos(TH1F* th, float var);  
  void checkEfficiency();
  void checkRate();
  int matchEGWithGenParticle();
  int matchEGWithSimTrack();
  // for L1TrackEmParticles
  edm::InputTag L1EGammaInputTag;
  edm::InputTag L1TkElectronsInputTag;

  int selectedEGTot_;
  int selectedEGTrkTot_;

  TH1F* etaEGamma_;
  TH1F* etEGamma_;
  TH1F* etaEGammaTrk_;
  TH1F* etEGammaTrk_;

  std::string analysisOption_;

  edm::SimTrackContainer simTracks_;
  l1extra::L1EmParticleCollection eGammaCollection;
  edm::Handle<L1TkElectronParticleCollection> L1TrackElectronsHandle;
  edm::Handle<L1EmParticleCollection> EGammaHandle;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  int ievent; 
};

L1TkElectronTrackObjectsAnalyzer::L1TkElectronTrackObjectsAnalyzer(const edm::ParameterSet& iConfig) {

  edm::Service<TFileService> fs;
  L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");

  analysisOption_ = iConfig.getParameter<std::string>("AnalysisOption");
  
}
void L1TkElectronTrackObjectsAnalyzer::beginJob() {
  edm::Service<TFileService> fs;

  etaEGamma_ = fs->make<TH1F>("EGammaEta_Ref","Eta of EGamma", 25, -2.5, 2.5);
  etaEGammaTrk_ = fs->make<TH1F>("EGammaEta_Track","Eta of TrkEGamma", 25, -2.5, 2.5);
    
  if (analysisOption_ == "Efficiency") {
    etEGamma_ = fs->make<TH1F>("EGammaEt_Ref","Et of EGamma", 24, 2.0, 50.0);
    etEGammaTrk_ = fs->make<TH1F>("EGammaEt_Track","Et of TrkEGamma", 24, 2.0, 50.0);
  } else {
    etEGamma_ = fs->make<TH1F>("EGammaEtThresholdEvt_Ref","Et of EGamma (EventEt threshold)", 24, 2.0, 50.0);
    etEGammaTrk_ = fs->make<TH1F>("EGammaEtThresholdEvt_Track","Et of TrkEGamma( Event Et threshold)", 24, 2.0, 50.0);
  }

  selectedEGTot_ = 0;
  selectedEGTrkTot_ = 0;
  ievent = 0;
}


L1TkElectronTrackObjectsAnalyzer::~L1TkElectronTrackObjectsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TkElectronTrackObjectsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  ievent++;  
  
  iEvent.getByLabel(L1TkElectronsInputTag, L1TrackElectronsHandle);
  

  iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
  eGammaCollection = (*EGammaHandle.product());
  sort(eGammaCollection.begin(), eGammaCollection.end(), L1TkElectron::EtComparator());
  
  
  Handle<edm::SimTrackContainer> simTrackHandle;
  iEvent.getByLabel( "g4SimHits", "", simTrackHandle);
  simTracks_ = (*simTrackHandle.product());
   
  iEvent.getByLabel("genParticles", genParticleHandle);

  if (analysisOption_ == "Efficiency") checkEfficiency();
  else checkRate();
}
void L1TkElectronTrackObjectsAnalyzer::endJob() {
  std::cout << " Selected EGammas " << selectedEGTot_ << std::endl;
  std::cout << " Selected Track EGammas " << selectedEGTrkTot_ << std::endl;
}
void L1TkElectronTrackObjectsAnalyzer::checkEfficiency() {

  int nSelectedEG = 0;
  int nSelectedEGTrk =0;
  int igIndx = matchEGWithGenParticle();
  if (igIndx == -1 ) return;
  
  std::cout<< " Event " << ievent << " Selected EGamma # " << igIndx << " Et " << eGammaCollection[igIndx].et() << " Eta " << eGammaCollection[igIndx].eta()<<std::endl;
  nSelectedEG++;
  etaEGamma_->Fill(eGammaCollection[igIndx].eta());
  std::vector<L1TkElectronParticle>::const_iterator egTrkIter ;
  float dRmin = 999.9;
  nSelectedEGTrk = 0;
  float eta_min;
  for (egTrkIter = L1TrackElectronsHandle -> begin(); egTrkIter != L1TrackElectronsHandle->end(); ++egTrkIter) {
    float dPhi = reco::deltaPhi(eGammaCollection[igIndx].phi(), egTrkIter->phi());
    float dEta = (eGammaCollection[igIndx].eta() - egTrkIter->eta());
    float dR =  sqrt(dPhi*dPhi + dEta*dEta);
    if (dR < dRmin) {
      dRmin = dR;
      eta_min = egTrkIter->eta();
    }
  }
  if (dRmin < 999.9) {
    nSelectedEGTrk++;
    if (nSelectedEGTrk == 1) {
      etaEGammaTrk_->Fill(eta_min);
      std::cout<< "Event # " << ievent << " Selected  EGamma matched " << std::endl;
    }
  }

  selectedEGTot_ += nSelectedEG;
  selectedEGTrkTot_ += nSelectedEGTrk;
}
void L1TkElectronTrackObjectsAnalyzer::checkRate() {

  int nSelectedEGTrk = 0;
  int nSelectedEG = 0;
  l1extra::L1EmParticleCollection::const_iterator egIter; 
  for (egIter = eGammaCollection.begin();  egIter != eGammaCollection.end(); ++egIter) {
    
    //    int ibx = egIter -> bx();
    //    if (ibx != 0) continue;
    
    float e_ele   = egIter->energy();
    float eta_ele = egIter->eta(); 
    float phi_ele = egIter->phi(); 
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if (fabs(eta_ele) > 2.3) continue;
    nSelectedEG++;
    if (nSelectedEG == 1) {
      fillIntegralHistos(etEGamma_, et_ele);
    }

    float et_min; 
    float dRmin = 999.9;
    std::vector<L1TkElectronParticle>::const_iterator egTrkIter ;
    for (egTrkIter = L1TrackElectronsHandle -> begin(); egTrkIter != L1TrackElectronsHandle->end(); ++egTrkIter) {
      float dPhi = reco::deltaPhi(phi_ele, egTrkIter->getEGRef()->phi());
      float dEta = (eta_ele - egTrkIter->getEGRef()->eta());
      float dR =  sqrt(dPhi*dPhi + dEta*dEta);
      if (dR < dRmin) {
	dRmin = dR;
	et_min = egTrkIter->pt();
      }
    }
    if (dRmin < 999.9) {
      nSelectedEGTrk++;
      if (nSelectedEGTrk == 1) {
	std::cout << "Selected EGammaTrk objet in the event " << et_ele << " with pt "<< egTrkIter->pt() << std::endl;           
	fillIntegralHistos(etEGammaTrk_, et_min);
      }
    }         
  }
  selectedEGTot_ += nSelectedEG;
  selectedEGTrkTot_ += nSelectedEGTrk;
}
void L1TkElectronTrackObjectsAnalyzer::fillIntegralHistos(TH1F* th, float var){
  int nbin = th->FindBin(var); 
  for (int ibin = 1; ibin < nbin+1; ibin++) th->Fill(th->GetBinCenter(ibin));
}
int L1TkElectronTrackObjectsAnalyzer::matchEGWithSimTrack() {
  int indx;
  float dRmin = 999.9;

  if ( fabs(simTracks_[0].momentum().eta())> 2.3 || simTracks_[0].momentum().pt() <= 20.0) return -1;

  for (unsigned int igam = 0; igam != eGammaCollection.size(); igam++) {
    int ibx = eGammaCollection[igam].bx();
    if (ibx != 0) continue;
    
    float e_ele   = eGammaCollection[igam].energy();
    float eta_ele = eGammaCollection[igam].eta(); 
    float phi_ele = eGammaCollection[igam].phi(); 
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if (fabs(eta_ele) > 2.3) continue;
    if ( et_ele <= 20) continue;
    float dPhi = reco::deltaPhi(simTracks_[0].momentum().phi(), phi_ele);
    float dEta = (simTracks_[0].momentum().eta() - eta_ele);
    float dR =  sqrt(dPhi*dPhi + dEta*dEta);
    if (dR < dRmin) {
      dRmin = dR;
      indx = igam; 
    } 
  } 
  if (dRmin < 0.1) return indx;
  else return -1;
}
int L1TkElectronTrackObjectsAnalyzer::matchEGWithGenParticle() {
  int indx;
  float dRmin = 999.9;
  const reco::Candidate & p = (*genParticleHandle)[0];
  if ( fabs(p.eta()) > 2.3 || p.pt() <= 20.0) return -1;
  for (unsigned int igam = 0; igam != eGammaCollection.size(); igam++) {
    int ibx = eGammaCollection[igam].bx();
    if (ibx != 0) continue;

    float eta_ele = eGammaCollection[igam].eta(); 
    float phi_ele = eGammaCollection[igam].phi(); 
    float e_ele   = eGammaCollection[igam].energy();
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if ( fabs(eta_ele) > 2.3 || et_ele <= 20.0) continue;
    float dPhi = reco::deltaPhi(p.phi(), phi_ele);
    float dEta = (p.eta() - eta_ele);
    float dR =  sqrt(dPhi*dPhi + dEta*dEta);
    if (dR < dRmin) {
      dRmin = dR;
      indx = igam; 
    } 
  } 
  if (dRmin < 0.1) {
    return indx;
  } 
  else return -1;
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1TkElectronTrackObjectsAnalyzer);
