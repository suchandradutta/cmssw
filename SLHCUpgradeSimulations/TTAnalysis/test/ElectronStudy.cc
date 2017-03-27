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
#include "SLHCUpgradeSimulations/TTAnalysis/interface/AnalObjects.h"
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

using namespace std;
using namespace edm;

const double twoPi=6.283185307;
const double caloThresh = 2.0;

namespace {
  typedef edm::Ref<edmNew::DetSetVector<TTStub<Ref_PixelDigi_ > >, TTStub<Ref_PixelDigi_> > stubRef;
  typedef std::vector< stubRef >  stubRefCollection;
  //  typedef std::vector<stubRef>::const_iterator stubIter;
  typedef stubRefCollection::iterator stubIter;
  typedef std::vector< TTTrack< Ref_PixelDigi_ > > L1TTTrackCollection;
  typedef std::vector<TTPixelTrack> L1TTPixelTrackCollection;
}
namespace TTStudy {
  bool compareStubLayer(const stubRef& s1, const stubRef & s2) {
    unsigned int l1, l2;
    StackedTrackerDetId stackedDetId1(s1->getDetId());
    StackedTrackerDetId stackedDetId2(s2->getDetId());
    if (stackedDetId1.isBarrel()) {
      l1 = stackedDetId1.iLayer();
    } else if (stackedDetId1.isEndcap()) {
      l1 = stackedDetId1.iSide() * 10 + stackedDetId1.iDisk();
    }
    if (stackedDetId2.isBarrel()) {
      l2 = stackedDetId2.iLayer();
    } else if (stackedDetId2.isEndcap()) {
      l2 = stackedDetId2.iSide() * 10 + stackedDetId2.iDisk();
    }
    return l1 < l2;
  }
  class EtComparator {
  public:
    bool operator()(const l1extra::L1EmParticle& a, const l1extra::L1EmParticle& b) const {
      double et_a = 0.0;
      double et_b = 0.0;    
      if (cosh(a.eta()) > 0.0) et_a = a.energy()/cosh(a.eta());
      if (cosh(b.eta()) > 0.0) et_b = b.energy()/cosh(b.eta());
      
      return et_a > et_b;
    }
  };
  class L1EGCrystalETComparator {
  public:
    bool operator()(const l1slhc::L1EGCrystalCluster a, const l1slhc::L1EGCrystalCluster b) const {
      double et_a = 0.0;
      double et_b = 0.0;
      if (cosh(a.eta()) > 0.0) et_a = a.energy()/cosh(a.eta());
      if (cosh(b.eta()) > 0.0) et_b = b.energy()/cosh(b.eta());

      return et_a > et_b;
    }
  };
}  


//
// class declaration
//

class ElectronStudy : public edm::EDAnalyzer {

public:

  explicit ElectronStudy(const edm::ParameterSet&);
  ~ElectronStudy();
  
      
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillSimTrackInfo();
  void fillRecTrackInfo();
  void fillElectronInfo(l1extra::L1EmParticleCollection::const_iterator igam);
  void matchStubEGamma(GlobalPoint epos, double ee, TTStudy::Electron& electron);

 
  int getSimTrkIndex(stubIter s); 
  int getSimTrkIndex(const unsigned int id);

  double getZ0(stubIter s1, stubIter s2);
  double getAngle(stubIter s1,stubIter s2,double z0);

  GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
  template<typename T> void removeDuplicates(std::vector<T>& vec);
  double getTwoPointPt(GlobalPoint a, GlobalPoint b);
  bool goodTwoPointZ(double innerR, double outerR, double innerZ, double outerZ );
  bool goodTwoPointPhi(double innerR, double outerR, double innerPhi, double outerPhi, double m_strength);
  double getDPhi(GlobalPoint epos, double eet, double r, double phi, double m_strength);
  double getZIntercept(GlobalPoint epos, double r, double z);
  double getPhiMiss(double eet, GlobalPoint spos1, GlobalPoint spos2);
  double getZMiss(GlobalPoint epos, double r1, double r2, double z1, double z2, bool bar);
  double getScaledZInterceptCut(unsigned int layer, double cut, double cfac, double eta);
  double getScaledZMissCut(int layer1, int layer2, double cut, double cfac, double eta);

  double getTrkEGamdP(GlobalPoint epos, L1TTTrackCollection::const_iterator trk);
  double getTrkEGamOriginZ(GlobalPoint epo, L1TTTrackCollection::const_iterator trk);
  double getTrkEGamDeltaR(GlobalPoint epos, L1TTTrackCollection::const_iterator trk);
  double getTrkEGamDeltaRCorr(GlobalPoint epos, stubIter s1, stubIter s2, L1TTTrackCollection::const_iterator trk);

  //  void addTrackletCandidate(L1TkStubIters& stubs,  L1TkStubIter& a_stub);
  int getMatchedGenIndex(GlobalPoint epos, double eet);
  unsigned int getLayerId(StackedTrackerDetId id);

  int getMotherParticleId(const reco::Candidate& gp);

  // helpers for EGamma Seed finding - From Avi
  double normalPhihalf(double phi) const {
    while (phi > M_PI) { phi -= M_PI; }
    while (phi < 0) { phi += M_PI; }
    return phi;
  }
  double normalPhi(double phi) const {
    while (phi > 2.* M_PI) { phi -= 2.*M_PI; }
    while (phi < 0) { phi += 2.*M_PI; }
    return phi;
  }
  double phiDiff(double phi1, double phi2){
    double result = normalPhi(phi1) - normalPhi(phi2);
    if(result > M_PI) result -= 2*M_PI;
    if(result < -M_PI) result += 2*M_PI;
    return result;
  }
  double unwrapPhi(double phi) const {
    while (phi > M_PI) { phi -= 2.*M_PI; }
    while (phi < -M_PI) { phi += 2.*M_PI; }
    return phi;
  }
  
  float getBreamStregth(float eta, float phi);
  // ----------member data ---------------------------
  double magnetStrength;
  unsigned int nLay_;  // number of layers
  double rcal_;

  // variables that fill ntuples:
  int nEvents_;      
  int npu_;
  int nstub_;

  int nElectron;
  int nGenTrack;
  int nRecoTrack;

  TTStudy::Event* eventBr_;
  std::vector<TTStudy::Electron>* electronsBr_; 
  std::vector<TTStudy::SimTrack>* simTracksBr_;  
  std::vector<TTStudy::Track>* tracksBr_;  
  std::vector<TTStudy::Track>* pxTracksBr_;  

  // ntuples:

  TTree* tree;


  // parameters that can be set in cfg file:
  unsigned int mcType_ ;
  edm::InputTag egammaSrc_;
  edm::InputTag trkSrc_;
  edm::InputTag pxTrkSrc_;
  edm::InputTag stubSrc_;
  edm::InputTag trkTruthSrc_;
  edm::InputTag stubTruthSrc_;
  bool debugFlag_;
  std::vector<double> radii_;
  std::vector<double> discs_;
  double egETThreshold_;
  std::string geometryOption_;
  double stubPtCut_;
  double trkEGdPhiCut_;
  double stubEGdPhiCut_;
  double stubEGdZCut_;
  double stubEGPhiMissCut_;
  double stubEGZMissCut_;

  std::vector<double> trkEGdRCut_;
  std::vector<double> trkEGEtoPCut_;

  std::vector< PCaloHit > theEcalSimHits;

  edm::Handle<reco::BeamSpot> beamSpotHandle;

  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
  const StackedTrackerGeometry*         theStackedTracker;
  ESHandle<TrackerGeometry> theTrkGeomHandle;

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > >    stubHandle;
  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > > mcTruthTTStubHandle;


  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  edm::Handle< std::vector< TTPixelTrack > > TTPixelTrackHandle;
  edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > mcTruthTTTrackHandle;

  edm::Handle<std::vector<TrackingParticle> > tpHandle;
  edm::SimTrackContainer simTracks_;
  edm::SimVertexContainer simVertices_;
  L1TTTrackCollection l1TTTracks_;
  l1slhc::L1EGCrystalClusterCollection clusColl_;
  L1TTPixelTrackCollection l1TTPixelTracks_;

};

//
// constructors and destructor
//
ElectronStudy::ElectronStudy(const edm::ParameterSet& iConfig)
{
  
  mcType_           = iConfig.getParameter<unsigned int>("mcType");
  egammaSrc_        = iConfig.getParameter<edm::InputTag>("egammaSrc");
  trkSrc_           = iConfig.getParameter<edm::InputTag>("trkSrc");
  trkTruthSrc_      = iConfig.getParameter<edm::InputTag>("trkTruthSrc");
  pxTrkSrc_           = iConfig.getParameter<edm::InputTag>("pixelTrkSrc");
  stubSrc_          = iConfig.getParameter<edm::InputTag>("stubSrc");
  stubTruthSrc_     = iConfig.getParameter<edm::InputTag>("stubTruthSrc");
  debugFlag_        = iConfig.getParameter<bool>("DebugFlag");
  egETThreshold_    = iConfig.getParameter<double>("EGammaETThreshold");
  geometryOption_   = iConfig.getParameter<std::string>("GeometryOption");
  stubPtCut_        = iConfig.getParameter<double>("StubPtCut");
  trkEGdPhiCut_     = iConfig.getParameter<double>("TrackEGammadPhiCut");
  trkEGdRCut_       = iConfig.getParameter< std::vector<double> >("TrackEGammadRCut");
  trkEGEtoPCut_     = iConfig.getParameter< std::vector<double> >("TrackEGammaEtoPCut");
  stubEGdPhiCut_    = iConfig.getParameter<double>("StubEGammadPhiCut");
  stubEGdZCut_      = iConfig.getParameter<double>("StubEGammadZCut");
  stubEGPhiMissCut_ = iConfig.getParameter<double>("StubEGammaPhiMissCut");
  stubEGZMissCut_   = iConfig.getParameter<double>("StubEGammaZMissCut");

  //now do what ever initialization is needed

  nEvents_ = 0;

  nLay_=10;

  if (geometryOption_ == "LB_6PS") radii_     = iConfig.getParameter< std::vector <double> > ("LongBarrelRadii");
  else if (geometryOption_ == "BE5D")radii_     = iConfig.getParameter< std::vector <double> > ("BarrelEndCap5DRadii");
  else radii_     = iConfig.getParameter< std::vector <double> > ("BarrelEndCapRadii");
  discs_     = iConfig.getParameter< std::vector <double> > ("EndCapDiscPositions");  
  rcal_= 129;	
  for (std::vector<double>::iterator it = radii_.begin(); it != radii_.end(); it++) {
    std::cout<< " Layer radius " << (*it) << std::endl; 
  }
}

ElectronStudy::~ElectronStudy()
{
  delete eventBr_;
  delete simTracksBr_; 
  delete tracksBr_; 
  delete electronsBr_; 
  delete pxTracksBr_;  

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ElectronStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  simTracksBr_->clear(); 
  tracksBr_->clear(); 
  electronsBr_->clear(); 
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
 
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  magnetStrength = theMagField.product()->inTesla(GlobalPoint(0,0,0)).z();
  //////////////////////////////////////////////////////////
  // Gun Particle Information from GenParticle
  //////////////////////////////////////////////////////////
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  iEvent.getByLabel("genParticles", genParticleHandle);
  int indx = 0;
  //  for (size_t i = 0; i < (*genParticleHandle.product()).size(); ++i) {
  //    const reco::Candidate & p = (*genParticleHandle)[i];
  //    if (abs(p.pdgId()) == 11 && p.status() == 1) {
  //      if (abs(getMotherParticleId(p)) == 24) indx = i;
  //      break;
  //    }
  //  }
  eventBr_->genParticleEt   = (*genParticleHandle)[indx].pt();
  eventBr_->genParticleEta  = (*genParticleHandle)[indx].eta();
  eventBr_->genParticlePhi = (*genParticleHandle)[indx].phi();
  //////////////////////////////////////////////////////////
  // Geometry 
  //////////////////////////////////////////////////////////

  iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
  theStackedTracker = stackedGeometryHandle.product(); /// Note this is different 

  //////////////////////////////////////////////////////////
  // Get Sim Tracks
  //////////////////////////////////////////////////////////
  Handle<edm::SimTrackContainer> simTrackHandle;
  if(mcType_ ==1 ) 
    iEvent.getByLabel( "g4SimHits", "", simTrackHandle);
  else
    iEvent.getByLabel( "famosSimHits", "", simTrackHandle);
  simTracks_ = (*simTrackHandle.product());
  ////////////////////////////////////////////////////////////
  // Get Sim Vertices
  ////////////////////////////////////////////////////////////
  Handle<edm::SimVertexContainer > simVertexHandle;

  if(mcType_ ==1 ) 
    iEvent.getByLabel( "g4SimHits","", simVertexHandle);
  else
    iEvent.getByLabel( "famosSimHits" , "" ,  simVertexHandle);
  simVertices_ =  (*simVertexHandle.product());
    
  ////////////////////////////////////////////////////////////////
  // Stubs From Track Trigger Hits
  ////////////////////////////////////////////////////////////////

  // TTStubs
  iEvent.getByLabel(stubSrc_, stubHandle);
  iEvent.getByLabel(stubTruthSrc_, mcTruthTTStubHandle );

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

  iEvent.getByLabel(trkTruthSrc_, mcTruthTTTrackHandle);
  fillSimTrackInfo();
  fillRecTrackInfo();

  Handle<l1extra::L1EmParticleCollection> CandHandle;
  iEvent.getByLabel(egammaSrc_,CandHandle);
  l1extra::L1EmParticleCollection eGammaCollection = (*CandHandle.product());
  sort(eGammaCollection.begin(), eGammaCollection.end(), TTStudy::EtComparator());

  std::cout << " EGammas " << eGammaCollection.size() << std::endl;


  edm::Handle<l1slhc::L1EGCrystalClusterCollection> L1EGCrystalClusters;
  iEvent.getByLabel("L1EGammaCrystalsProducer","EGCrystalCluster", L1EGCrystalClusters);
  clusColl_ = (*L1EGCrystalClusters.product());   
  std::sort(clusColl_.begin(), clusColl_.end(), TTStudy::L1EGCrystalETComparator());

  for ( l1extra::L1EmParticleCollection::const_iterator EGammaIter = eGammaCollection.begin(); 
	EGammaIter != eGammaCollection.end();  ++EGammaIter) {

    fillElectronInfo(EGammaIter);
   
  }
  tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
ElectronStudy::beginJob()
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Framework handles for the EVENTSETUP tracker geometry, L1 stack geometry, etc...
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("ElectronTriggerInfo","Electron Trigger Information");
  
  eventBr_ = new TTStudy::Event();
  tree->Branch("Event", "TTStudy::Event", &eventBr_, 32000, 2);

  electronsBr_ = new std::vector<TTStudy::Electron>(); 
  tree->Branch("Electron", "std::vector<TTStudy::Electron>", &electronsBr_);

  simTracksBr_ = new std::vector<TTStudy::SimTrack>(); 
  tree->Branch("SimTrack", "std::vector<TTStudy::SimTrack>", &simTracksBr_);
  
  tracksBr_ = new std::vector<TTStudy::Track>(); 
  tree->Branch("Track", "std::vector<TTStudy::Track>", &tracksBr_);

  pxTracksBr_ = new std::vector<TTStudy::Track>(); 
  tree->Branch("PixelTrack", "std::vector<TTStudy::Track>", &pxTracksBr_);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronStudy::endJob() {
}
int ElectronStudy::getSimTrkIndex(stubIter s) {
  int index = -5;
  
  if (mcTruthTTStubHandle->isGenuine((*s))) {
    const std::vector<SimTrack>& simTks = mcTruthTTStubHandle->findTrackingParticlePtr((*s))->g4Tracks();
    if (simTks.size() > 0) {
      index = getSimTrkIndex(simTks[0].trackId());
    }
  }
  return index;
}
int ElectronStudy::getSimTrkIndex(unsigned int id) {
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



double ElectronStudy::getZ0(stubIter s1, stubIter s2) {
  GlobalPoint hitpos1 = theStackedTracker->findGlobalPosition(&(*(*s1)));
  GlobalPoint hitpos2 = theStackedTracker->findGlobalPosition(&(*(*s2)));
  
  double r1x = hitpos1.perp();
  double r2x = hitpos2.perp();
  double z1x = hitpos1.z();
  double z2x = hitpos2.z();
  
  return (z1x - r1x*(z2x-z1x)/(r2x-r1x));
}

double ElectronStudy::getAngle(stubIter s1, stubIter s2, double z0) {
  
  GlobalPoint primary(theStackedTracker->findGlobalPosition(&(*(*s1))).x(),
		      theStackedTracker->findGlobalPosition(&(*(*s1))).y(),
		      theStackedTracker->findGlobalPosition(&(*(*s1))).z() - z0);
  GlobalPoint secondary(theStackedTracker->findGlobalPosition(&(*(*s2))).x(),
			theStackedTracker->findGlobalPosition(&(*(*s2))).y(),
			theStackedTracker->findGlobalPosition(&(*(*s2))).z() - z0);
  
  double angle = (primary.x()*secondary.x()+
		  primary.y()*secondary.y()+
		  primary.z()*secondary.z())/
    primary.mag()/secondary.mag();
  
  return angle;
}


GlobalPoint ElectronStudy::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0; 
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 ) 
    { 
      double ecalZ = 314*fabs(eta)/eta;
      
      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {
      double rperp = 129.0;
      double zface =  sqrt( cos( theta ) * cos( theta ) /
			    ( 1 - cos( theta ) * cos( theta ) ) * 
			    rperp * rperp ) * fabs( eta ) / eta;  
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}

int ElectronStudy::getMatchedGenIndex(GlobalPoint epos, double eet){
  double dRMin = 9999.9; 
  int indx = -1;
  for(int it = 0; it < (int)simTracks_.size(); it++) {
    float dR = deltaR(epos.eta(), epos.phi(), simTracks_[it].momentum().eta(), simTracks_[it].momentum().phi());
    if (dR < dRMin) {
      dRMin = dR;
      indx = it;
    }
  }  
  return indx;
}

template<typename T>
void ElectronStudy::removeDuplicates(std::vector<T>& vec)
{
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


double ElectronStudy::getTwoPointPt(GlobalPoint a, GlobalPoint b)
{
  double phi_a = a.phi();
  double phi_b = b.phi();
  double r_a = a.perp();
  double r_b = b.perp();
  
  return fabs((r_b-r_a)/(phi_b-phi_a))*(2.998e8*3.8)/(2e11);     
}

bool ElectronStudy::goodTwoPointZ(double innerR, double outerR, double innerZ, double outerZ ) {
    
  double mIPWidth = 200.0;
  double positiveZBoundary = (mIPWidth - outerZ) * (outerR - innerR);
  double negativeZBoundary = -(mIPWidth + outerZ) * (outerR - innerR);
  double multipliedLocation = (innerZ - outerZ) * outerR;
    
    
  if( multipliedLocation < positiveZBoundary &&
      multipliedLocation > negativeZBoundary )
    return true;
  return false;
}
bool ElectronStudy::goodTwoPointPhi(double innerR, double outerR, double innerPhi, double outerPhi, double m_strength) {
  
  // Rebase the angles in terms of 0-2PI, should
  if ( innerPhi < 0.0 ) innerPhi += 2.0 * TMath::Pi();
  if ( outerPhi < 0.0 ) outerPhi += 2.0 * TMath::Pi();
    
  // Check for seed compatibility given a pt cut
  // Threshold computed from radial location of hits
    double mCompatibilityScalingFactor =
      (100.0 * 2.0e+9 * 2.0) / (TMath::C() * m_strength);
    
    mCompatibilityScalingFactor = 1.0 / mCompatibilityScalingFactor;
    
    double deltaPhiThreshold =
      (outerR - innerR) * mCompatibilityScalingFactor;
    
    // Delta phi computed from hit phi locations
    double deltaPhi = outerPhi - innerPhi;
    if(deltaPhi<0) deltaPhi = -deltaPhi;
    
    if(deltaPhi<deltaPhiThreshold) return true;
    else return false;
}
double ElectronStudy::getDPhi(GlobalPoint epos, double eet, double r, double phi, double m_strength) {
    
  double er = epos.perp();
    
  double phiVsRSlope = -3.00e-3 * m_strength / eet / 2.;
    
  // preselecton variable
  double psi = reco::deltaPhi(phi,epos.phi());
  double deltaPsi = psi - (er-r)*phiVsRSlope;
  double antiDeltaPsi = psi - (r-er)*phiVsRSlope;
  double dP;
  if (fabs(deltaPsi)<fabs(antiDeltaPsi)){
    dP = deltaPsi;
  }else{
    dP = antiDeltaPsi;
  }
  return dP;
}
double ElectronStudy::getZIntercept(GlobalPoint epos, double r, double z) {

  double er = epos.perp();
  double ez = epos.z();
        
  double zint = (er*z - r*ez)/(er-r);
  return zint;
}

// Extrapolate a L1Track from its vertex to the ECal Cluster position::method
double ElectronStudy::getTrkEGamdP(GlobalPoint epos, L1TTTrackCollection::const_iterator trk) {
  double er = epos.perp();

  // Using track fit curvature
  //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
  // double curv = 0.003 * magnetStrength / trk->getMomentum().perp(); 
  double curv = trk->getRInv();
  double x1 = (asin(er*curv/(2.0)));
  double phi1 = phiDiff(trk->getMomentum().phi(), epos.phi());

  double dif1 = phi1 - x1;
  double dif2 = phi1 + x1; 
  if (fabs(dif1) < fabs(dif2)) return dif1;
  else return dif2; 
}

double ElectronStudy::getTrkEGamOriginZ(GlobalPoint epos, L1TTTrackCollection::const_iterator trk) {
	  
  double er = epos.perp();
  double ez = epos.z();
  
  //  TVector3 trkvtx;
  //  trkvtx.SetXYZ(trk->vx(), trk->vy(), trk->vz());
  double r1 = trk->getPOCA().perp();
  double z1 = trk->getPOCA().z();

  double originZ1 = (er*z1 - r1*ez)/(er-r1);
  return originZ1;
}

double ElectronStudy::getTrkEGamDeltaR(GlobalPoint epos, L1TTTrackCollection::const_iterator trk){
  double dPhi = phiDiff(epos.phi(), trk->getMomentum().phi());
  //double dPhi = (ephi - trk->getMomentum().phi());
  //if (fabs(dPhi) > 3.14285)
    //dPhi = (trk->getMomentum().phi() - ephi);
  double dEta = (epos.eta() - trk->getMomentum().eta());
  double deltaR = sqrt(dPhi*dPhi + dEta*dEta);
  return deltaR;
}
double ElectronStudy::getTrkEGamDeltaRCorr(GlobalPoint epos, stubIter s1, stubIter s2, L1TTTrackCollection::const_iterator trk){
  double er = epos.perp();
  double ez = epos.z();
  double z0 = getZ0(s1,s2);
  double theta = 0;
  if (ez >= 0)
    theta = atan(er/fabs(ez-z0));
  else theta = 3.14 - atan(er/fabs(ez-z0));
  double corr_eta = -log(tan(theta/2.0));

  double dPhi = phiDiff(epos.phi(), trk->getMomentum().phi());
  double dEta = (corr_eta - trk->getMomentum().eta());
  double deltaR = sqrt(dPhi*dPhi + dEta*dEta);
  return deltaR;
}

double ElectronStudy::getPhiMiss(double eet, GlobalPoint spos1, GlobalPoint spos2) {

  double pT = eet;
  double curv = pT*100*.877;
  if (curv == 0) return 999.9;
    
    
  double r1 = spos1.perp();
  double r2 = spos2.perp();
    
  double phi1 = spos1.phi();
  double phi2 = spos2.phi();
    
  //Predict phi of hit 2
  double a = (r2-r1)/(2*curv);
  double b = reco::deltaPhi(phi2,phi1);
    
  double phiMiss = 0;
  if(fabs(b - a)<fabs(b + a)) phiMiss = b - a;
  if(fabs(b - a)>fabs(b + a)) phiMiss = b + a;

  return phiMiss;

}

// ZMiss
double ElectronStudy::getZMiss(GlobalPoint epos, double r1, double r2, double z1, double z2, bool bar) {

  double er = epos.perp();
  double ez = epos.z();
    
    
  double missVal ;
  if (bar) {
    missVal = z2 - (r2*(ez-z1)-r1*ez +
		    er*z1)/(er-r1);
  } else {
    missVal = r2 - er - (er-r1)*(z2-ez)/(ez-z1);
  }
  return missVal;
}
double ElectronStudy::getScaledZInterceptCut(unsigned int layer, double cut, double cfac, double eta) {
  double mult = (rcal_-radii_[layer])/(rcal_-radii_[1]);
  return (mult*(cut+cfac*(1.0/(1.0-cos(2*atan(exp(-1.0*fabs(eta))))))));
}
// Z-issCut
double ElectronStudy::getScaledZMissCut(int layer1, int layer2, double cut, double cfac, double eta) {
  double mult = ( (radii_[layer2] - radii_[layer1])/(rcal_-radii_[layer1]) )*( (rcal_ -radii_[1])/(radii_[2]-radii_[1]));
  return (mult*(cut+cfac*(1.0/(1.0-cos(2*atan(exp(-1.0*fabs(eta))))))));
}
unsigned int ElectronStudy::getLayerId(StackedTrackerDetId id) {
  unsigned int layer = 999999;
  if (id.isBarrel()) {
    layer = id.iLayer();
  } else if (id.isEndcap()) {
    layer = id.iSide() * 10 + id.iDisk();
  } else {
    std::cout << "Neither Barrel nor Endcap " << layer << std::endl;
  }
  if (layer > 100) std::cout << " Wrong Layer " << layer << std::endl;
  return layer;
}
//
// -- Fill SimTrack Information
//
void ElectronStudy::fillSimTrackInfo() {
  for (edm::SimTrackContainer::const_iterator it = simTracks_.begin(); it != simTracks_.end(); ++it) {
    TTStudy::SimTrack simTk;
  
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
void ElectronStudy::fillRecTrackInfo() {
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
    TTStudy::Track recTk;
    // Track Pt = 0.3* B * radius
    recTk.pt        = trk->getMomentum().perp();
    recTk.ptFromStub= pTFrom2Stubs::pTFrom2( trk, theStackedTracker);
    recTk.eta       = trk->getMomentum().eta();
    recTk.phi       = trk->getMomentum().phi();
    recTk.chiSquare = trk->getChi2();
    recTk.curvature = trk->getRInv();
    //    recTk.curvature = 0.3 * magnetStrength /trk->getMomentum().perp();
    recTk.vertexX   = trk->getPOCA().x();
    recTk.vertexY   = trk->getPOCA().y();
    recTk.vertexZ   = trk->getPOCA().z();
    recTk.nStub     = trk->getStubRefs().size();
    recTk.pdgId     = PDGId;
    recTk.vertexId  = VtxId;
    recTk.matchedSimTrack = mcTruthTTTrackHandle->isGenuine(l1track_ptr);
    tracksBr_->push_back(recTk);
  }
  //Store Pixel Tracks
  L1TTPixelTrackCollection::const_iterator iterL1PixelTrack;
  for ( L1TTPixelTrackCollection::const_iterator pxTrk = l1TTPixelTracks_.begin(); pxTrk != l1TTPixelTracks_.end(); pxTrk++ ) {
    TTStudy::Track recPxTk;
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
//
// -- Fill EGamma Properties
//
void ElectronStudy::fillElectronInfo(l1extra::L1EmParticleCollection::const_iterator igam) {
  bool selected = true;
  float e_ele   = igam->energy();
  float eta_ele = igam->eta();
  float et_ele = 0;
  if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
  else et_ele = -1.0;     

  if (et_ele <= egETThreshold_) selected = false; 
    
  if (!selected) return;
  
  TTStudy::Electron elec; 
  
  elec.phi = igam->phi();
  elec.eta = eta_ele;
  elec.e   = e_ele;

  GlobalPoint epos = getCalorimeterPosition(igam->phi(), igam->eta(), igam->energy());

  elec.x = epos.x();
  elec.y = epos.y();
  elec.z = epos.z();
  elec.r = epos.perp();
  elec.et = et_ele;
  

  elec.simTkIndx = getMatchedGenIndex(epos, et_ele);    
  elec.bStrength = getBreamStregth(eta_ele,igam->phi());  

  if (debugFlag_) std::cout << " Et, E, Eta and Phi " << et_ele << " " << e_ele << " " <<  eta_ele << " " << epos.phi() << std::endl;
  
  matchStubEGamma(epos, e_ele, elec);    
  electronsBr_->push_back(elec);
}  
//
// -- Match Track and EGamma Clusters
//
void ElectronStudy::matchStubEGamma(GlobalPoint epos, double ee, TTStudy::Electron& electron) {

  double et_ele = 0;
  if (cosh(epos.eta()) > 0.0) et_ele = ee/cosh(epos.eta());
  else et_ele = -1.0;      
  
  // Get stubs consistent with EGamma
  stubRefCollection preSelectedStubs; 
  //  std::chrono::time_point<std::chrono::system_clock> start, end;
  //  start = std::chrono::system_clock::now();
  for (edmNew::DetSetVector< TTStub<Ref_PixelDigi_> >::const_iterator it  = stubHandle->begin();
       it != stubHandle->end();++it) {
    for (edmNew::DetSet<TTStub<Ref_PixelDigi_> >::const_iterator jt  = it->begin();
	 jt != it->end(); ++jt) {
      /// Make the reference 
      stubRef stub_ref = edmNew::makeRefTo(stubHandle, jt);

      StackedTrackerDetId stackedDetId(stub_ref->getDetId());

      // Store Track information in maps, skip if the Cluster is not good
      //      bool isGenuine = mcTruthTTStubHandle->isGenuine(stub_ref);
      //      if (!isGenuine) continue;
      float stub_pt = theStackedTracker->findRoughPt(magnetStrength,&(*stub_ref)); 

      if (stubPtCut_ > 0.0 && stub_pt <= stubPtCut_) continue;
      unsigned int ilayer = getLayerId(stackedDetId);
      double originZMultiplier = 1.0;
      double r   = stackedGeometryHandle->findGlobalPosition(&(*stub_ref)).perp();
      double phi = stackedGeometryHandle->findGlobalPosition(&(*stub_ref)).phi();
      double z   = stackedGeometryHandle->findGlobalPosition(&(*stub_ref)).z();
      
      double dPhi = getDPhi(epos, et_ele, r, phi, magnetStrength);
      double zIntercept = getZIntercept(epos, r, z);
      double scaledZIntCut;
      double scaledDPhiCut;
      
      if (ilayer < 10) {
	originZMultiplier = radii_[ilayer]/radii_[1];
	scaledZIntCut = originZMultiplier*(stubEGdZCut_+0.75*(1.0/(1.0-cos(2*atan(exp(-1.0*TMath::Abs(epos.eta()))))))); 
	scaledDPhiCut = stubEGdPhiCut_;
      } else{
	scaledDPhiCut = 1.6*stubEGdPhiCut_;
	scaledZIntCut = stubEGdZCut_;
      }
      if ( (fabs(dPhi) < scaledDPhiCut) && 
	   (fabs(zIntercept)<scaledZIntCut) ) preSelectedStubs.push_back(stub_ref);
    }
  }
  //  end = std::chrono::system_clock::now();

  //  std::chrono::duration<double> elapsed_seconds = end-start;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  //  std::cout << "finished computation at " << std::ctime(&end_time)
  //	    << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
  std::sort(preSelectedStubs.begin(), preSelectedStubs.end(), TTStudy::compareStubLayer);
  // Get two-point tracklets consistent with this EGamma
  // and loop over them
  for (stubIter istub1 = preSelectedStubs.begin(); istub1 != preSelectedStubs.end(); istub1++) {
    for (stubIter istub2 = istub1+1; istub2 != preSelectedStubs.end(); istub2++) {
      
      StackedTrackerDetId stackedDetId1((*istub1)->getDetId());
      StackedTrackerDetId stackedDetId2((*istub2)->getDetId());
      
      unsigned layer1 = getLayerId(stackedDetId1);

      unsigned layer2 = getLayerId(stackedDetId2);
      if (layer1 >= layer2) continue;
      bool barrel = true;
      if (layer2 > 10) barrel = false; 

      double innerZ = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).z();
      double outerZ = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).z();
      double innerR = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).perp();
      double outerR = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).perp();
      double innerPhi = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).phi();
      double outerPhi = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).phi();
      
      const GlobalPoint beamSpot(beamSpotHandle->x0(),beamSpotHandle->y0(),beamSpotHandle->z0());
      GlobalPoint s1pos(stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).x()-beamSpot.x(),
			stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).y()-beamSpot.y(),
			stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).z()-beamSpot.z());
      
      GlobalPoint s2pos(stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).x()-beamSpot.x(),
			stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).y()-beamSpot.y(),
			stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).z()-beamSpot.z());
      
      
      //    if (debugFlag_) cout<<"Found "<<twoPointCands.size()<<" matched 2-point tracklets."<<endl;
      //    for(vector <L1TkStubIters>::const_iterator ip = twoPointCands.begin(); ip != twoPointCands.end(); ip++) {
      
      if (layer1 > 100 || layer2 > 100) std::cout << " Wrong Layers " << layer1 << " " << layer2 << std::endl;
      
      //      int trk_index1 = getSimTrkIndex(istub1);
      //      int trk_index2 = getSimTrkIndex(istub2);
      int trk_index1 = -5;
      int trk_index2 = -5;
    
      if(!goodTwoPointPhi(innerR, outerR, innerPhi, outerPhi, magnetStrength)) continue;
      if(!goodTwoPointZ(innerR, outerR, innerZ, outerZ)) continue;
      
      double zMiss = getZMiss(epos, innerR, outerR, innerZ, outerZ, barrel);
      double phiMiss = getPhiMiss(et_ele, s1pos, s2pos);
      /*      double phiMissScaledCut = stubEGPhiMissCut_;
      if (fabs(epos.eta()) >= 1.1) {
	if (layer1 <= 3 && layer2 <= 3) phiMissScaledCut *= 1.4;
	else phiMissScaledCut *= 1.8;
      }
      double zMissScaledCut;
      if (barrel) {
	zMissScaledCut = getScaledZMissCut(layer1, layer2,stubEGPhiMissCut_, 0.04, epos.eta());
      } else {
	zMissScaledCut = 2.0;
      }
      if(fabs(phiMiss)< phiMissScaledCut && fabs(zMiss) < zMissScaledCut) {*/
      TTStudy::Tracklet tracklet;
      tracklet.deltaPhi1       = getDPhi(epos, et_ele, innerR, innerPhi, magnetStrength);
      tracklet.zIntercept1     = getZIntercept(epos, innerR, innerZ);
      tracklet.layer1          = layer1;
      tracklet.phi1            = theStackedTracker->findGlobalPosition(&(*(*istub1))).phi();
      tracklet.eta1            = theStackedTracker->findGlobalPosition(&(*(*istub1))).eta();
      tracklet.x1              = theStackedTracker->findGlobalPosition(&(*(*istub1))).x();
      tracklet.y1              = theStackedTracker->findGlobalPosition(&(*(*istub1))).y();
      tracklet.z1              = theStackedTracker->findGlobalPosition(&(*(*istub1))).z();
      tracklet.r1              = theStackedTracker->findGlobalPosition(&(*(*istub1))).perp();
      tracklet.pt1             = theStackedTracker->findRoughPt(magnetStrength,&(*(*istub1)));
      tracklet.trackIndex1     = trk_index1;
      if(trk_index1 != -5)  { // not combinatorial 
	tracklet.particleId1 = simTracks_[trk_index1].type();
	tracklet.truePt1     = simTracks_[trk_index1].momentum().pt();
      } 
      tracklet.deltaPhi2       = getDPhi(epos, et_ele, outerR, outerPhi, magnetStrength);
      tracklet.zIntercept2     = getZIntercept(epos, outerR, outerZ);
      tracklet.layer2          = layer2;
      tracklet.phi2            = theStackedTracker->findGlobalPosition(&(*(*istub2))).phi();
      tracklet.eta2            = theStackedTracker->findGlobalPosition(&(*(*istub2))).eta();
      tracklet.x2              = theStackedTracker->findGlobalPosition(&(*(*istub2))).x();
      tracklet.y2              = theStackedTracker->findGlobalPosition(&(*(*istub2))).y();
      tracklet.z2              = theStackedTracker->findGlobalPosition(&(*(*istub2))).z();
      tracklet.r2              = theStackedTracker->findGlobalPosition(&(*(*istub2))).perp();
      tracklet.pt2             = theStackedTracker->findRoughPt(magnetStrength,&(*(*istub2)));
      tracklet.trackIndex2     = trk_index2;
      if(trk_index2 != -5)  { // not combinatorial 
	tracklet.particleId2   = simTracks_[trk_index2].type();
	tracklet.truePt2       = simTracks_[trk_index2].momentum().pt();
      }
      tracklet.phiMiss = phiMiss;
      tracklet.zMiss   = zMiss;
      GlobalPoint a = theStackedTracker->findGlobalPosition(&(*(*istub1)));
      GlobalPoint b = theStackedTracker->findGlobalPosition(&(*(*istub2)));
      tracklet.twoPointPt = getTwoPointPt(a,b);
      tracklet.twoPointZIntercept = getZ0(istub1, istub2);
      electron.matchedTracklets.push_back(tracklet);
      //    }
    }
  }
} 
int ElectronStudy::getMotherParticleId(const reco::Candidate& gp) {
  int mid = -1;
  if (gp.numberOfMothers() == 0) return mid;
  const reco::Candidate* m0 = gp.mother(0);
  if (!m0)  return mid;

  mid = m0->pdgId();
  while (gp.pdgId() == mid) {

    const reco::Candidate* m = m0->mother(0);
    if (!m) {
      mid = -1;
      break;
    }   
    mid = m->pdgId();
    m0 = m;
  }
  return mid;
}
float ElectronStudy::getBreamStregth(float eta, float phi) {
  float dr_min = 999.9;
  size_t indx_min = 9999;
  float bs = -1.0;
  for (size_t i = 0; i < clusColl_.size(); i++) {
    float dr = reco::deltaR(clusColl_[i].eta(), clusColl_[i].phi(), eta, phi);
    if (dr < dr_min) {
      indx_min = i;
      dr_min = dr;
    }
  }
  if (indx_min != 999.9) bs = clusColl_[indx_min].bremStrength(); 
  return bs;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronStudy);
