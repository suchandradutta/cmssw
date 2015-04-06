/*********************************/
/*********************************/
/**                             **/
/** Stacked Tracker Simulations **/
/**        Laura  Fields        **/
/**             2009            **/
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
#include "SimDataFormats/SLHC/interface/L1CaloCluster.h"
#include "SimDataFormats/SLHC/interface/L1CaloTower.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTPixelTrack.h"

#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h" 

#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithm.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithm_a.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithmRecord.h"
#include "SLHCUpgradeSimulations/TTIAnalysis/interface/AnalObjects.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/pTFrom2Stubs.h"

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
#include "SimDataFormats/SLHC/interface/L1EGCrystalCluster.h"
#include "Utilities/Timing/interface/TimingReport.h"

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <sstream>



namespace {
  typedef edm::Ref<edmNew::DetSetVector<TTStub<Ref_PixelDigi_ > >, TTStub<Ref_PixelDigi_> > stubRef;
  typedef std::vector< stubRef >  stubRefCollection;
  //  typedef std::vector<stubRef>::const_iterator stubIter;
  typedef stubRefCollection::iterator stubIter;
  typedef std::vector< TTTrack< Ref_PixelDigi_ > > L1TTTrackCollection;
  typedef std::vector<TTPixelTrack> L1TTPixelTrackCollection;
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
  void fillRecTrackInfo();
 
  bool matchGenInfo(const reco::Candidate* a, const reco::Candidate* b);
  int getSimTrkIndex(const unsigned int id);

  // ----------member data ---------------------------
  // variables that fill ntuples:
  int nEvents_;      
  int npu_;
  int nstub_;

  int nRecoTrack;

  TTIStudy::Event* eventBr_;
  std::vector<TTIStudy::GenParticle>* genParBr_;  
  std::vector<TTIStudy::SimTrack>* simTracksBr_;  
  std::vector<TTIStudy::Track>* tracksBr_;  
  std::vector<TTIStudy::Track>* pxTracksBr_;  

  // ntuples:

  TTree* tree;


  // parameters that can be set in cfg file:
  edm::InputTag trkSrc_;
  edm::InputTag trkTruthSrc_;
  edm::InputTag pxTrkSrc_;
  bool debugFlag_;

  edm::Handle<reco::BeamSpot> beamSpotHandle;

  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
  const StackedTrackerGeometry*         theStackedTracker;
  edm::ESHandle<TrackerGeometry> theTrkGeomHandle;


  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  edm::Handle< std::vector< TTPixelTrack > > TTPixelTrackHandle;
  edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > mcTruthTTTrackHandle;

  edm::Handle<std::vector<TrackingParticle> > tpHandle;
  edm::SimTrackContainer simTracks_;
  edm::SimVertexContainer simVertices_;
  reco::GenParticleCollection genParticles_;
  L1TTTrackCollection l1TTTracks_;
  L1TTPixelTrackCollection l1TTPixelTracks_;
};

//
// constructors and destructor
//
BSToPhiPhiStudy::BSToPhiPhiStudy(const edm::ParameterSet& iConfig)
{
  
  trkSrc_           = iConfig.getParameter<edm::InputTag>("trkSrc");
  trkTruthSrc_      = iConfig.getParameter<edm::InputTag>("trkTruthSrc");
  pxTrkSrc_           = iConfig.getParameter<edm::InputTag>("pixelTrkSrc");
  debugFlag_        = iConfig.getParameter<bool>("DebugFlag");

  //now do what ever initialization is needed

  nEvents_ = 0;

}

BSToPhiPhiStudy::~BSToPhiPhiStudy()
{
  delete eventBr_;
  delete genParBr_; 
  delete simTracksBr_; 
  delete tracksBr_; 
  delete pxTracksBr_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BSToPhiPhiStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  genParBr_->clear(); 
  simTracksBr_->clear(); 
  tracksBr_->clear(); 
  pxTracksBr_->clear(); 

  nEvents_++;
  std::cout << " Event # " << nEvents_ << std::endl;

  eventBr_->event = iEvent.id().event(); 
  //////////////////////////////////////////////////////////
  // Pile Up Vertices
  //////////////////////////////////////////////////////////
  int npu = 0;
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    for (std::vector<PileupSummaryInfo>::const_iterator PVI  = PupInfo->begin(); 
                                                        PVI != PupInfo->end(); ++PVI) {
      // More info about PU is here:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
       if (PVI->getBunchCrossing() == 0) npu_ = PVI->getPU_NumInteractions();
    } 
  }
  eventBr_->nPileUp = npu;
  //////////////////////////////////////////////////////////
  // Extract Beamspot, Mag Field, etc
  //////////////////////////////////////////////////////////
  iEvent.getByLabel("BeamSpotFromSim", "BeamSpot", beamSpotHandle);
  eventBr_->beamSpotX0 = beamSpotHandle->x0();
  eventBr_->beamSpotY0 = beamSpotHandle->y0();
  eventBr_->beamSpotZ0 = beamSpotHandle->z0();
 
  //////////////////////////////////////////////////////////
  // Gun Particle Information from GenParticle
  //////////////////////////////////////////////////////////
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  iEvent.getByLabel("genParticles", genParticleHandle);
  genParticles_ = (*genParticleHandle.product());
  int indx = 0;
  eventBr_->genParticleEt   = (*genParticleHandle)[indx].pt();
  eventBr_->genParticleEta  = (*genParticleHandle)[indx].eta();
  eventBr_->genParticlePhi  = (*genParticleHandle)[indx].phi();
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////

  iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
  theStackedTracker = stackedGeometryHandle.product(); /// Note this is different 

  //////////////////////////////////////////////////////////
  // Get Sim Tracks
  //////////////////////////////////////////////////////////
  Handle<edm::SimTrackContainer> simTrackHandle;
  iEvent.getByLabel( "g4SimHits", "", simTrackHandle);
  simTracks_ = (*simTrackHandle.product());
  ////////////////////////////////////////////////////////////
  // Get Sim Vertices
  ////////////////////////////////////////////////////////////
  Handle<edm::SimVertexContainer > simVertexHandle;
  iEvent.getByLabel( "g4SimHits","", simVertexHandle);
  simVertices_ =  (*simVertexHandle.product());
    
  // TrackingParticles
  iEvent.getByLabel("mix", "MergedTrackTruth", tpHandle);

  ////////////////////////////////////////////////////////////////
  // L1Tracks From Track Trigger Hits
  ////////////////////////////////////////////////////////////////

  iEvent.getByLabel(trkSrc_, TTTrackHandle);
  l1TTTracks_ = (*TTTrackHandle.product());
  std::cout << "Found " << l1TTTracks_.size() << " L1 Tracks" << std::endl;

  iEvent.getByLabel(pxTrkSrc_, TTPixelTrackHandle);
  l1TTPixelTracks_ = (*TTPixelTrackHandle.product());
  std::cout << "Found " << l1TTPixelTracks_.size() << " L1 Pixel Tracks" << std::endl;
  //  std::cout << "Found " << TTPixelTrackHandle->size() << " L1 Pixel Tracks" << std::endl;
  

  iEvent.getByLabel(trkTruthSrc_, mcTruthTTTrackHandle);
  fillGenParticleInfo();
  fillSimTrackInfo();
  fillRecTrackInfo();

  tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
BSToPhiPhiStudy::beginJob()
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Framework handles for the EVENTSETUP tracker geometry, L1 stack geometry, etc...
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("ElectronTriggerInfo","Electron Trigger Information");
  
  eventBr_ = new TTIStudy::Event();
  tree->Branch("Event", "TTIStudy::Event", &eventBr_, 32000, 2);

  genParBr_ = new std::vector<TTIStudy::GenParticle>(); 
  tree->Branch("GenParticle", "std::vector<TTIStudy::GenParticle>", &genParBr_);

  simTracksBr_ = new std::vector<TTIStudy::SimTrack>(); 
  tree->Branch("SimTrack", "std::vector<TTIStudy::SimTrack>", &simTracksBr_);
  

  tracksBr_ = new std::vector<TTIStudy::Track>(); 
  tree->Branch("Track", "std::vector<TTIStudy::Track>", &tracksBr_);

  pxTracksBr_ = new std::vector<TTIStudy::Track>(); 
  tree->Branch("PixelTrack", "std::vector<TTIStudy::Track>", &pxTracksBr_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BSToPhiPhiStudy::endJob() {
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
  int i = 0;
  reco::GenParticleCollection::const_iterator gbeg = genParticles_.begin(); 
  if (debugFlag_ > 0) {
    std::cout << std::setiosflags(std::ios::fixed);
    std::cout << std::setprecision(2);
    std::cout << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
	      << std::endl;
  }
  for (reco::GenParticleCollection::const_iterator it = genParticles_.begin(); it != genParticles_.end(); ++it) {
    TTIStudy::GenParticle genPar;
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
    TTIStudy::SimTrack simTk;
  
    simTk.pt = it->momentum().pt();
    simTk.eta = it->momentum().eta();
    simTk.phi = it->momentum().phi();
    int vertIndex = it->vertIndex();
    simTk.vtxIndx = vertIndex; 
    if((int)simVertices_.size()>vertIndex) {
      simTk.vx = simVertices_[vertIndex].position().x();
      simTk.vy = simVertices_[vertIndex].position().y();
      simTk.vz = simVertices_[vertIndex].position().z();
    }
    simTk.type = it->type(); 
    simTracksBr_->push_back(simTk);
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
void BSToPhiPhiStudy::fillRecTrackInfo() {

  // Store Pixel Tracks
  int itrk = 0;
  for(L1TTTrackCollection::const_iterator trk = l1TTTracks_.begin(); trk != l1TTTracks_.end() ; ++trk){

    edm::Ptr< TTTrack< Ref_PixelDigi_ > > l1track_ptr(TTTrackHandle, itrk);
    itrk++;
    int PDGId = -9999;
    int VtxId = -1;
    edm::Ptr< TrackingParticle > tpPtr = mcTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);
    if (!tpPtr.isNull()) {
      PDGId = tpPtr->pdgId();
      if (tpPtr->g4Tracks().size() > 0) VtxId = tpPtr->g4Tracks()[0].vertIndex();
    }
    TTIStudy::Track recTk;
    // Track Pt = 0.3* B * radius
    recTk.pt        = trk->getMomentum().perp();
    recTk.ptFromStub= pTFrom2Stubs::pTFrom2( trk, theStackedTracker);
    recTk.eta       = trk->getMomentum().eta();
    recTk.phi       = trk->getMomentum().phi();
    recTk.chiSquare = trk->getChi2();
    recTk.curvature = trk->getRInv();
    //    recTk.curvature = 0.3 * magnetStrength /trk->getMomentum().perp();
    recTk.vertexX   = trk->getPOCA(5).x();
    recTk.vertexY   = trk->getPOCA(5).y();
    recTk.vertexZ   = trk->getPOCA(5).z();
    recTk.nStub     = trk->getStubRefs().size();
    recTk.pdgId     = PDGId;
    recTk.vertexId  = VtxId;
    recTk.matchedSimTrack = mcTruthTTTrackHandle->isGenuine(l1track_ptr);
    tracksBr_->push_back(recTk);
  }

  //Store Pixel Tracks
  L1TTPixelTrackCollection::const_iterator iterL1PixelTrack;
  for ( L1TTPixelTrackCollection::const_iterator pxTrk = l1TTPixelTracks_.begin(); pxTrk != l1TTPixelTracks_.end(); pxTrk++ ) {
    TTIStudy::Track recPxTk;
    recPxTk.pt        = pxTrk->getMomentum().perp();
    recPxTk.eta       = pxTrk->getMomentum().eta();
    recPxTk.phi       = pxTrk->getMomentum().phi();
    recPxTk.chiSquare = pxTrk->getChi2();
    recPxTk.curvature = pxTrk->getRInv();
    recPxTk.vertexX   = pxTrk->getPOCA().x();
    recPxTk.vertexY   = pxTrk->getPOCA().y();
    recPxTk.vertexZ   = pxTrk->getPOCA().z();
    pxTracksBr_->push_back(recPxTk);
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(BSToPhiPhiStudy);
