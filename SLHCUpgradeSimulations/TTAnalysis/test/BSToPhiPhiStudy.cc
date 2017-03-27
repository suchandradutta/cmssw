/*********************************/
/*********************************/
/**                             **/
/**   BsToPhiPhi Study          **/
/**                             **/
/**   suchandra.dutta@cern.ch   **/
/**                             **/
/*********************************/
/*********************************/
// system include files
#include <memory>
#include <cmath>
#include <typeinfo>
#include <vector>
#include <algorithm> 
#include <chrono>
#include <ctime>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "L1Trigger/TrackTrigger/interface/TTClusterAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTClusterAlgorithmRecord.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h" 

#include "SLHCUpgradeSimulations/TTAnalysis/interface/AnalObjects.h"
//#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/pTFrom2Stubs.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <sstream>
using namespace l1extra;

namespace {
  typedef edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_ > >, TTStub<Ref_Phase2TrackerDigi_> > stubRef;
  typedef std::vector< stubRef >  stubRefCollection;
  typedef stubRefCollection::iterator stubIter;
  typedef std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > L1TTTrackCollection;
}

//
// class declaration
//
class BSToPhiPhiStudy : public edm::EDAnalyzer {

public:

  explicit BSToPhiPhiStudy(const edm::ParameterSet&);
  ~BSToPhiPhiStudy();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillGenParticleInfo();
  void fillSimTrackInfo();
  void fillL1TrackInfo();
  void fillL1PixelTrackInfo();
  void fillOfflineTrackInfo();
  void fillL1TkJetInfo();
  void fillL1TkMuonInfo();
 
  bool matchGenInfo(const reco::Candidate* a, const reco::Candidate* b);
  int getSimTrkIndex(const unsigned int id);

  // ----------member data ---------------------------
  // variables that fill ntuples:
  int nEvents_;      
  int nstub_;

  int nRecoTrack;
  TTStudy::Event* eventBr_;
  std::vector<TTStudy::GenParticle>* genParBr_;  
  std::vector<TTStudy::SimTrack>* simTracksBr_;  
  std::vector<TTStudy::Track>* tracksBr_;  
  std::vector<TTStudy::Track>* pxTracksBr_;  
  std::vector<TTStudy::Track>* recoTracksBr_;  
  std::vector<TTStudy::L1Jet>* l1JetsBr_;  
  std::vector<TTStudy::L1Muon>* l1MuonsBr_;  

  // ntuples:
  TTree* tree;

  // parameters that can be set in cfg file:
  const edm::InputTag trkSrc_;
  const edm::InputTag trkTruthSrc_;
  const edm::InputTag pxTrkSrc_;
  const edm::InputTag recoTrkSrc_;
  const edm::InputTag beamSpotSrc_;
  const edm::InputTag offVertexSrc_;
  const edm::InputTag l1VertexSrc_;
  const edm::InputTag l1TkJetSrc_;
  const edm::InputTag l1TkMuonSrc_;
  
  const bool debugFlag_;
  const bool l1TrackFlag_;
  const bool simTrackFlag_;
  const bool l1VertexFlag_;
  const bool l1JetFlag_;
  const bool l1MuonFlag_;
  const bool pixelTrackFlag_;
  const bool recoTrackFlag_;
  const bool genParticleFlag_;

  const edm::EDGetTokenT<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackToken_;
  //const edm::EDGetTokenT<L1TkPrimaryVertexCollection> l1VertexToken_;
  //const edm::EDGetTokenT<L1TkJetParticleCollection> l1TkJetToken_;
  //const edm::EDGetTokenT<L1TkMuonParticleCollection> l1TkMuonToken_;
  //  const edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_ > > mcTruthTTTrackToken_;
  const edm::EDGetTokenT<std::vector<TrackingParticle> > tpToken_;
  const edm::EDGetTokenT<reco::TrackCollection> recoTrackToken_;
  const edm::EDGetTokenT<reco::VertexCollection>  offlineVertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot >  beamSpotToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer > simTrackToken_;
  const edm::EDGetTokenT<edm::SimVertexContainer > simVertexToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection >  genParticleToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;

  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<TrackerGeometry> theTrkGeomHandle;

  edm::Handle<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle_;
  /*  edm::Handle<L1TkPrimaryVertexCollection> l1VertexHandle_;
  edm::Handle<L1TkJetParticleCollection> l1TkJetHandle_;
  edm::Handle<L1TkMuonParticleCollection> l1TkMuonHandle_;
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_ > > mcTruthTTTrackHandle_;*/
  edm::Handle<std::vector<TrackingParticle> > tpHandle_;
  edm::Handle<reco::TrackCollection> recoTrackHandle_;
  edm::Handle<reco::VertexCollection> offlineVertexHandle_;
  edm::Handle<reco::BeamSpot> beamSpotHandle_;
  edm::Handle<edm::SimTrackContainer > simTrackHandle_;
  edm::Handle<edm::SimVertexContainer > simVertexHandle_;
  edm::Handle<reco::GenParticleCollection> genParticleHandle_;
  edm::Handle<std::vector<PileupSummaryInfo> > puHandle_;
  
  edm::ESHandle<TrackerTopology> tTopoHandle_;

  reco::GenParticleCollection genParticles_;
  reco::TrackCollection recoTracks_;
  edm::SimTrackContainer simTracks_;
  edm::SimVertexContainer simVertices_;
};

//
// constructors and destructor
//
BSToPhiPhiStudy::BSToPhiPhiStudy(const edm::ParameterSet& iConfig) :
  // Collection src
  trkSrc_(iConfig.getParameter<edm::InputTag>("trkSrc")),
  //  trkTruthSrc_(iConfig.getParameter<edm::InputTag>("trkTruthSrc")),
  pxTrkSrc_(iConfig.getParameter<edm::InputTag>("pixelTrkSrc")),
  recoTrkSrc_(iConfig.getParameter<edm::InputTag>("recoTrkSrc")),
  beamSpotSrc_(iConfig.getParameter<edm::InputTag>("beamspotSrc")),
  offVertexSrc_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  l1VertexSrc_(iConfig.getParameter<edm::InputTag>("l1VertexSrc")),
  l1TkJetSrc_(iConfig.getParameter<edm::InputTag>("l1TkJetSrc")),
  l1TkMuonSrc_(iConfig.getParameter<edm::InputTag>("l1TkMuonSrc")),
  // flags
  debugFlag_       (iConfig.getParameter<bool>("DebugFlag")),
  l1TrackFlag_     (iConfig.getParameter<bool>("l1TrackFlag")),
  simTrackFlag_    (iConfig.getParameter<bool>("simTrackFlag")),
  l1VertexFlag_    (iConfig.getParameter<bool>("l1VertexFlag")),
  l1JetFlag_       (iConfig.getParameter<bool>("l1JetFlag")),
  l1MuonFlag_      (iConfig.getParameter<bool>("l1MuonFlag")),
  pixelTrackFlag_  (iConfig.getParameter<bool>("pixelTrackFlag")),
  recoTrackFlag_   (iConfig.getParameter<bool>("recoTrackFlag")),
  genParticleFlag_ (iConfig.getParameter<bool>("genParticleFlag")),
  // tokens
  TTTrackToken_(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (trkSrc_)),
  // l1VertexToken_(consumes< L1TkPrimaryVertexCollection >(l1VertexSrc_)),
  //l1TkJetToken_(consumes< L1TkJetParticleCollection>(l1TkJetSrc_)),
  //l1TkMuonToken_(consumes<L1TkMuonParticleCollection>(l1TkMuonSrc_)),
  //  mcTruthTTTrackToken_(consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_ > >(trkTruthSrc_)),
  tpToken_(consumes< std::vector< TrackingParticle > >(edm::InputTag("mix", "MergedTrackTruth"))),
  recoTrackToken_(consumes <reco::TrackCollection> (recoTrkSrc_)),
  offlineVertexToken_(consumes< reco::VertexCollection >(offVertexSrc_)),
  beamSpotToken_(consumes< reco::BeamSpot >(beamSpotSrc_)),
  simTrackToken_(consumes < edm::SimTrackContainer >(edm::InputTag ("g4SimHits", ""))),
  simVertexToken_(consumes < edm::SimVertexContainer >(edm::InputTag ("g4SimHits",""))),
  genParticleToken_(consumes < reco::GenParticleCollection >(edm::InputTag ("genParticles"))),
  puToken_(consumes <std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo")))
{  
  // now do whatever initialization is needed
  nEvents_ = 0;
  if (debugFlag_) std::cout << " pixelTrackFlag_ " << pixelTrackFlag_ << std::endl;
}

BSToPhiPhiStudy::~BSToPhiPhiStudy()
{
  delete eventBr_;
  delete tracksBr_; 
  if (simTrackFlag_) delete simTracksBr_; 
  if (genParticleFlag_) delete genParBr_; 
  if (pixelTrackFlag_) delete pxTracksBr_;
  if (recoTrackFlag_) delete recoTracksBr_;
  if (l1JetFlag_) delete l1JetsBr_;
  if (l1MuonFlag_) delete l1MuonsBr_;
}

// ------------ method called to for each event  ------------
void
BSToPhiPhiStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Clear the branches and get ready for the next event
  tracksBr_->clear(); 
  if (simTrackFlag_) simTracksBr_->clear(); 
  if (genParticleFlag_) genParBr_->clear(); 
  if (pixelTrackFlag_) pxTracksBr_->clear(); 
  if (recoTrackFlag_) recoTracksBr_->clear();
  if (l1JetFlag_) l1JetsBr_->clear();
  if (l1MuonFlag_) l1MuonsBr_->clear();

  nEvents_++;
  if (debugFlag_) std::cout << " Event # " << nEvents_ << std::endl;

  eventBr_->event = iEvent.id().event(); 
  //////////////////////////////////////////////////////////
  // Pile Up Vertices
  //////////////////////////////////////////////////////////
  int npu = 0;
  if (!iEvent.isRealData()) {
    iEvent.getByToken(puToken_, puHandle_);
    for (std::vector<PileupSummaryInfo>::const_iterator PVI  = puHandle_->begin(); 
                                                        PVI != puHandle_->end(); ++PVI) {
      // More info about PU is here:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
      if (debugFlag_) std::cout << "bx = "  << PVI->getBunchCrossing() << ", nPU = " << PVI->getPU_NumInteractions() << std::endl;
      if (PVI->getBunchCrossing() == 0) npu = PVI->getPU_NumInteractions(); // in-time PU
      npu++;
    } 
  }
  eventBr_->nPileUp = npu;
 
  //////////////////////////////////////////////////////////
  // Extract Beamspot, Mag Field, etc
  //////////////////////////////////////////////////////////
  /*  iEvent.getByToken(beamSpotToken_, beamSpotHandle_);
  eventBr_->beamSpotX0 = beamSpotHandle_->position().x();
  eventBr_->beamSpotY0 = beamSpotHandle_->position().y();
  eventBr_->beamSpotZ0 = beamSpotHandle_->position().z();*/

  //////////////////////////////////////////////////////////
  // Gun Particle Information from GenParticle
  //////////////////////////////////////////////////////////
  iEvent.getByToken(genParticleToken_, genParticleHandle_);
  if (genParticleHandle_.isValid()) {
    int indx = 0;
    eventBr_->genParticleEt  = (*genParticleHandle_)[indx].pt();
    eventBr_->genParticleEta = (*genParticleHandle_)[indx].eta();
    eventBr_->genParticlePhi = (*genParticleHandle_)[indx].phi();

    if (genParticleFlag_) {
      genParticles_ = (*genParticleHandle_.product());
      fillGenParticleInfo();
    }
  }
  else {
    std::cerr << "analyze: GenParticleCollection for InputTag genParticles not found!" << std::endl; 
  }
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////
  // Tracker Topology 
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle_);

  if (simTrackFlag_) {
    //////////////////////////////////////////////////////////
    // Get Sim Tracks and Sim Vertices
    //////////////////////////////////////////////////////////
    iEvent.getByToken(simTrackToken_, simTrackHandle_);
    if (simTrackHandle_.isValid()) {
      simTracks_ = (*simTrackHandle_.product());
      iEvent.getByToken(simVertexToken_, simVertexHandle_);
      if (simVertexHandle_.isValid()) {
	fillSimTrackInfo();
      } else {
	std::cerr << "analyze: SimVertex Collection for InputTag g4SimHits not found!" << std::endl; 
      }
    } else {
      std::cerr << "analyze: SimTrack Collection for InputTag g4SimHits not found!" << std::endl; 
    }
  }
  ////////////////////////////////////////////////////////////
  // TrackingParticles
  ////////////////////////////////////////////////////////////
  //  iEvent.getByLabel("mix", "MergedTrackTruth", tpHandle_);

  ////////////////////////////////////////////////////////////////
  // L1Tracks From Track Trigger Hits
  ////////////////////////////////////////////////////////////////
  if(l1TrackFlag_){
    iEvent.getByToken(TTTrackToken_, TTTrackHandle_);
    //    iEvent.getByToken(mcTruthTTTrackToken_, mcTruthTTTrackHandle_);
    //    if (TTTrackHandle_.isValid() && mcTruthTTTrackHandle_.isValid()) fillL1TrackInfo();
    if (TTTrackHandle_.isValid() ) fillL1TrackInfo();

    else {
      std::cerr << "analyze: L1 Track Collection for " << trkSrc_ << " or MC Truth Collection for " << trkTruthSrc_ << " not found!" << std::endl; 
    } 
  }

  ////////////////////////////////////////////////////////////////
  // Offline Reconstructed Tracks and Vertex
  ///////////////////////////////////////////////////////////////
  if (recoTrackFlag_) {
    
    iEvent.getByToken(recoTrackToken_, recoTrackHandle_);
    if (recoTrackHandle_.isValid()) {

      iEvent.getByToken(offlineVertexToken_, offlineVertexHandle_);
      if (offlineVertexHandle_.isValid()) {
	const reco::Vertex& vit = offlineVertexHandle_->front();
	eventBr_->nOffPV = (*offlineVertexHandle_.product()).size();
	eventBr_->zOffPV = vit.z();
	eventBr_->sumOffPV = vit.p4().P();

	fillOfflineTrackInfo();
      }
      else {
	std::cerr << "analyze: offlinePrimaryVertex Collection for InputTag " << offVertexSrc_ << " not found!" << std::endl; 
      } 
    }
    else {
      std::cerr << "analyze: offline Track Collection for InputTag " << recoTrkSrc_ << " not found!" << std::endl; 
    } 
  }
  /*
  ////////////////////////////////////////////////////////////
  // L1 Primary Vertices
  ////////////////////////////////////////////////////////////
  if (l1VertexFlag_) {
    iEvent.getByToken(l1VertexToken_, l1VertexHandle_);
    int nvtx = 0;
    if ( l1VertexHandle_.isValid() ) {
      std::vector<L1TkPrimaryVertex>::const_iterator vtxIter = l1VertexHandle_->begin();  // only one algorithm is run in the L1TkPrimaryVertexProducer					                                                        // (in contrast to earlier, under-dev, versions of the code)
      eventBr_->zL1PV   = vtxIter->getZvertex();
      eventBr_->sumL1PV = vtxIter->getSum();
    }
    else {
      std::cerr << "analyze: L1TPrimaryVertex Collection for InputTag " << l1VertexSrc_ << " not found!" << std::endl; 
    } 
    eventBr_->nL1PV = 1;
  }
  ////////////////////////////////////////////////////////////////
  // L1TkJet
  ///////////////////////////////////////////////////////////////
  if (l1JetFlag_) { 
    iEvent.getByToken(l1TkJetToken_, l1TkJetHandle_);
    if ( l1TkJetHandle_.isValid() ) fillL1TkJetInfo();
    else {
      std::cerr << "analyze: L1TkJet Collection for InputTag " << l1TkJetSrc_ << " not found!" << std::endl; 
    }
  }
  ////////////////////////////////////////////////////////////////
  // L1TkMuon
  ///////////////////////////////////////////////////////////////
  if (l1MuonFlag_) {
    iEvent.getByLabel(l1TkMuonToken_, l1TkMuonHandle_);
    if ( l1TkMuonHandle_.isValid() ) fillL1TkMuonInfo();
    else {
      std::cerr << "analyze: L1TkMuon Collection for InputTag " << l1TkMuonSrc_ << " not found!" << std::endl; 
    }
  }
  */
  tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void BSToPhiPhiStudy::beginJob()
{
  ////////////////////////////////////////////////////////////////////////////////
  // Framework handles for the EVENTSETUP tracker geometry, L1 stack geometry, etc...
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("ElectronTriggerInfo","Electron Trigger Information");
  
  eventBr_ = new TTStudy::Event();
  tree->Branch("Event", "TTStudy::Event", &eventBr_, 32000, 2);

  if (l1TrackFlag_) {
    tracksBr_ = new std::vector<TTStudy::Track>(); 
    tree->Branch("L1Track", "std::vector<TTStudy::Track>", &tracksBr_);
  }
  if (genParticleFlag_) {
    genParBr_ = new std::vector<TTStudy::GenParticle>(); 
    tree->Branch("GenParticle", "std::vector<TTStudy::GenParticle>", &genParBr_);
  }

  // Sim Tracks
  if (simTrackFlag_) {
    simTracksBr_ = new std::vector<TTStudy::SimTrack>(); 
    tree->Branch("SimTrack", "std::vector<TTStudy::SimTrack>", &simTracksBr_);
  }
  if (pixelTrackFlag_) {
    pxTracksBr_ = new std::vector<TTStudy::Track>(); 
    tree->Branch("PixelTrack", "std::vector<TTStudy::Track>", &pxTracksBr_);
  }

  if (recoTrackFlag_) {
    recoTracksBr_ = new std::vector<TTStudy::Track>(); 
    tree->Branch("RecoTrack", "std::vector<TTStudy::Track>", &recoTracksBr_);
  }
  // L1Jet
  if (l1JetFlag_) {
    l1JetsBr_ = new std::vector<TTStudy::L1Jet>(); 
    tree->Branch("L1Jet", "std::vector<TTStudy::L1Jet>", &l1JetsBr_);
  }

  // L1Muon
  if (l1MuonFlag_) {
    l1MuonsBr_ = new std::vector<TTStudy::L1Muon>(); 
    tree->Branch("L1Muon", "std::vector<TTStudy::L1Muon>", &l1MuonsBr_);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void BSToPhiPhiStudy::endJob() {
}

int BSToPhiPhiStudy::getSimTrkIndex(unsigned int id) {
  int index = -5;
  for(unsigned int i = 0; i<simTracks_.size(); i++) {
    unsigned int simTrkId = simTracks_[i].trackId();
    if (simTrkId == id) {
      index = i;
      break;
    }
  }
  return index;
}

//
// -- Fill GenParticle Information
//
void BSToPhiPhiStudy::fillGenParticleInfo() {
  reco::GenParticleCollection::const_iterator gbeg = genParticles_.begin(); 
  if (debugFlag_ > 0) {
    std::cout << std::setiosflags(std::ios::fixed);
    std::cout << std::setprecision(2);
    std::cout << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
	      << std::endl;
  }
  unsigned int i = 0;
  for (reco::GenParticleCollection::const_iterator it = genParticles_.begin(); it != genParticles_.end(); ++it) {
    TTStudy::GenParticle genPar;
    genPar.eta = it->eta();
    genPar.phi = it->phi();
    genPar.p = it->p();
    genPar.px = it->px();
    genPar.py = it->py();
    genPar.pz = it->pz();
    genPar.pt = it->pt();
    genPar.energy = it->energy();
    genPar.pdgId = it->pdgId();
    genPar.vx = it->vx();
    genPar.vy = it->vy();
    genPar.vz = it->vz();
    genPar.status = it->status();
    genPar.charge = it->charge();

    // First mother
    int idm = -1;
    const reco::Candidate* m = it->mother();
    if (m != nullptr) {
      for (reco::GenParticleCollection::const_iterator mit = genParticles_.begin(); mit != genParticles_.end(); ++mit) {
        const reco::Candidate* ap = &(*mit);
        if (matchGenInfo(m, ap)) {
	  idm = std::distance(gbeg, mit);
  	  break;
        }
      }
    }
    genPar.motherIndex = idm;

    std::ostringstream dID;
    for (size_t j = 0; j < it->numberOfDaughters(); ++j) {
      const reco::Candidate* d = it->daughter(j);
      for (reco::GenParticleCollection::const_iterator dit = genParticles_.begin(); dit != genParticles_.end(); ++dit) {
        const reco::Candidate* ap = &(*dit);  
	if (matchGenInfo(d, ap)) {
	  int idd = std::distance(gbeg, dit);
	  genPar.daughterIndices.push_back(idd);
	  dID << " " << idd;
	  break;
	}
      }
    }
    if (debugFlag_ > 0) {
      std::string ds = dID.str();
      if (!ds.length()) ds = " -";
      std::cout << std::setw(4)  << i++
		<< std::setw(8)  << it->status()
		<< std::setw(10) << it->pdgId()
		<< std::setw(8)  << it->charge()
		<< std::setw(10) << it->eta()
		<< std::setw(9)  << it->phi()
		<< std::setw(9)  << it->pt()
		<< std::setw(9)  << it->energy()
		<< std::setw(16) << idm
		<< ds
		<< std::endl;
    }
    genParBr_->push_back(genPar);
  }
  std::cout << std::resetiosflags(std::ios::fixed);
}
//
// -- Fill SimTrack Information
//
void BSToPhiPhiStudy::fillSimTrackInfo() {
  for (edm::SimTrackContainer::const_iterator it = simTracks_.begin(); it != simTracks_.end(); ++it) {
    TTStudy::SimTrack simTk;
  
    simTk.pt = it->momentum().pt();
    simTk.eta = it->momentum().eta();
    simTk.phi = it->momentum().phi();
    int vertIndex = it->vertIndex();
    simTk.vtxIndx = vertIndex; 
    if (static_cast<int>(simVertices_.size()) > vertIndex) {
      simTk.vx = simVertices_[vertIndex].position().x();
      simTk.vy = simVertices_[vertIndex].position().y();
      simTk.vz = simVertices_[vertIndex].position().z();
    }
    simTk.type = it->type(); 
    simTracksBr_->push_back(simTk);
  }
}
//
// -- Fill Offline Track Information
//
void BSToPhiPhiStudy::fillOfflineTrackInfo() {
  reco::TrackCollection recoTracks = *(recoTrackHandle_.product());

  for (reco::TrackCollection::const_iterator trk = recoTracks.begin(); trk != recoTracks.end(); ++trk) {
    TTStudy::Track offlineTk;
    offlineTk.pt = trk->pt();
    offlineTk.eta = trk->eta();
    offlineTk.phi = trk->phi();
    offlineTk.curvature = trk->qoverp();
    offlineTk.chiSquare = trk->chi2();
    offlineTk.vertexX = trk->vx();
    offlineTk.vertexY = trk->vy();
    offlineTk.vertexZ = trk->vz();
    double trkd0 = trk->d0();
    double trkdz = trk->dz();
    if (beamSpotHandle_.isValid()) {
      trkd0 = -(trk->dxy(beamSpotHandle_->position()));
      trkdz = trk->dz(beamSpotHandle_->position());
    }
    offlineTk.d0 = trkd0;
    offlineTk.z0 = trkdz;

    double dxyWrtPV = 999;
    double dzWrtPV = 999;
    if (offlineVertexHandle_.isValid()) {
      const reco::Vertex& vit = offlineVertexHandle_->front();
      dxyWrtPV = trk->dxy(vit.position());
      dzWrtPV = trk->dz(vit.position());
    }
    offlineTk.d0PV = dxyWrtPV;
    offlineTk.z0PV = dzWrtPV;
    
    recoTracksBr_->push_back(offlineTk);  
  }
}
bool BSToPhiPhiStudy::matchGenInfo(const reco::Candidate* a, const reco::Candidate* b) {
  bool result = false;
  if ( a->pdgId()  == b->pdgId()  &&
       a->status() == b->status() &&        
       a->pt()     == b->pt()     &&        
       a->eta()    == b->eta()    &&       
       a->phi()    == b->phi() ) result = true;
  return result;
}
//
// -- Fill L1 Track Information
//
void BSToPhiPhiStudy::fillL1TrackInfo() {
  L1TTTrackCollection l1TTTracks_ = *(TTTrackHandle_.product());
  if (debugFlag_) std::cout << "Found " << l1TTTracks_.size() << " L1 Tracks" << std::endl;

  const TrackerTopology* tTopo = tTopoHandle_.product();

  int itrk = 0;
  for (L1TTTrackCollection::const_iterator trk = l1TTTracks_.begin(); trk != l1TTTracks_.end() ; ++trk) {
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle_, itrk);
    itrk++;
    int PDGId = -9999;
    int VtxId = -1;
    /*    edm::Ptr< TrackingParticle > tpPtr = mcTruthTTTrackHandle_->findTrackingParticlePtr(l1track_ptr);
    if (!tpPtr.isNull()) {
      PDGId = tpPtr->pdgId();
      if (tpPtr->g4Tracks().size() > 0) VtxId = tpPtr->g4Tracks()[0].vertIndex();
      }*/
    TTStudy::Track recTk;
    // Track Pt = 0.3* B * radius
    recTk.pt        = trk->getMomentum(5).perp();
    recTk.eta       = trk->getMomentum(5).eta();
    recTk.phi       = trk->getMomentum(5).phi();
    recTk.chiSquare = trk->getChi2(5);
    recTk.chiSquareRed = trk->getChi2Red(5);
    recTk.curvature = trk->getRInv(5);
    //    recTk.curvature = 0.3 * magnetStrength /trk->getMomentum().perp();
    recTk.vertexX   = trk->getPOCA(5).x();
    recTk.vertexY   = trk->getPOCA(5).y();
    recTk.vertexZ   = trk->getPOCA(5).z();
    
    recTk.d0        = trk->getPOCA(5).x()*sin(trk->getMomentum(5).phi()) + 
                      trk->getPOCA(5).y()*cos(trk->getMomentum(5).phi());
    recTk.z0        = trk->getPOCA(5).z();
    
    recTk.pt_p4        = trk->getMomentum(4).perp();
    recTk.eta_p4       = trk->getMomentum(4).eta();
    recTk.phi_p4       = trk->getMomentum(4).phi();
    recTk.chiSquare_p4 = trk->getChi2(4);
    recTk.chiSquareRed_p4 = trk->getChi2Red(4);
    recTk.curvature_p4 = trk->getRInv(4);
    recTk.vertexX_p4   = trk->getPOCA(4).x();
    recTk.vertexY_p4   = trk->getPOCA(4).y();
    recTk.vertexZ_p4   = trk->getPOCA(4).z();
    
    stubRefCollection stubs = trk->getStubRefs();
    recTk.nStub = stubs.size();
    
    recTk.nStub_PS = 0;
    recTk.nStub_SS = 0;
    for (stubIter it = stubs.begin(); it != stubs.end(); ++it) {
      DetId detid = (*it)->getDetId();
      if (detid.det() != DetId::Detector::Tracker) continue;
      if (detid.subdetId() == StripSubdetector::TOB) {
	(tTopo->tobLayer(detid) <= 3) ? recTk.nStub_PS++ : recTk.nStub_SS++;
      }	else if (detid.subdetId() == StripSubdetector::TID) {
	(tTopo->tidRing(detid) <= 9) ? recTk.nStub_PS++ : recTk.nStub_SS++;
      }
    }
    recTk.pdgId     = PDGId;
    recTk.vertexId  = VtxId;
    //    recTk.matchedSimTrack = mcTruthTTTrackHandle_->isGenuine(l1track_ptr);
    
    tracksBr_->push_back(recTk);
  }
}
/*void BSToPhiPhiStudy::fillL1TkJetInfo() {
  for (std::vector<L1TkJetParticle>::const_iterator jetIter = l1TkJetHandle_->begin(); jetIter != l1TkJetHandle_->end(); ++jetIter) {
    TTStudy::L1Jet l1Jet;
    l1Jet.pt = jetIter->pt();
    l1Jet.phi = jetIter->phi();
    l1Jet.eta = jetIter->eta();
    l1Jet.zvtx = jetIter->getJetVtx(); // Z position
    const edm::Ref<L1JetParticleCollection> jetRef = jetIter->getJetRef();
    if (jetRef.isNonnull()) l1Jet.et = jetRef->et();
    
    const std::vector<edm::Ptr<L1TkJetParticle::L1TkTrackType > > trkPtrs = jetIter->getTrkPtrs() ;
    //    L1JetParticle::JetType type = Jetref -> type();
    
    float sum_pt = 0;
    float sum_wt_vtx = 0;
    for (size_t it = 0; it < trkPtrs.size(); it++) {
      edm::Ptr<L1TkJetParticle::L1TkTrackType> aTrack = trkPtrs.at(it);
      sum_pt += aTrack -> getMomentum().perp();
      sum_wt_vtx += aTrack->getPOCA(5).z() * aTrack->getMomentum().perp();
    }
    if (sum_pt) sum_wt_vtx /= sum_pt;
    l1Jet.zvtx_tk = sum_wt_vtx;  
    l1Jet.nTk = trkPtrs.size();
    
    l1JetsBr_->push_back(l1Jet);  
  }
}
void BSToPhiPhiStudy::fillL1TkMuonInfo() {
  if (debugFlag_) std::cout << " -----  Accessing L1TkMuonPaticle objects ---- " << std::endl;
  for (std::vector<L1TkMuonParticle>::const_iterator muIter = l1TkMuonHandle_->begin(); muIter != l1TkMuonHandle_->end(); ++muIter) {
    TTStudy::L1Muon l1Muon;
    
    l1Muon.pt = muIter->pt();
    l1Muon.eta = muIter->eta();
    l1Muon.phi = muIter->phi();
    if (muIter->getMuRef().isNonnull()) l1Muon.isMip = muIter->getMuRef()->isMip() ? 1: 0;
    l1Muon.isolation = muIter->getTrkIsol();
    l1Muon.z0 = muIter->getTrkzVtx();
    
    const edm::Ptr<L1TkMuonParticle::L1TkTrackType> & trk = muIter->getTrkPtr() ;
    l1Muon.pt_tk  = trk->getMomentum(5).perp();
    l1Muon.eta_tk = trk->getMomentum(5).eta();
    l1Muon.phi_tk = trk->getMomentum(5).phi();
    l1Muon.chiSquare = trk->getChi2(5);
    l1Muon.chiSquareRed = trk->getChi2Red(5);
    l1Muon.curvature = trk->getRInv(5);
    l1Muon.vertexX   = trk->getPOCA(5).x();
    l1Muon.vertexY   = trk->getPOCA(5).y();
    l1Muon.d0        = trk->getPOCA(5).x()*sin(trk->getMomentum(5).phi()) + 
                       trk->getPOCA(5).y()*cos(trk->getMomentum(5).phi());
    
    l1MuonsBr_->push_back(l1Muon);  
    }
    }*/
//define this as a plug-in
DEFINE_FWK_MODULE(BSToPhiPhiStudy);
