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

// system include files
#include <memory>
#include <typeinfo>
#include <vector>
#include <algorithm> 

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

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h" 

#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithm.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithm_a.h"
#include "SLHCUpgradeSimulations/L1TrackTrigger/interface/ClusteringAlgorithmRecord.h"
#include "SLHCUpgradeSimulations/ElectronAnalysis/interface/AnalObjects.h"

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

#include <iostream>
#include <math.h>
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"

using namespace std;
using namespace edm;

const double twoPi=6.283185307;
const double caloThresh = 2.0;
namespace TTIStudy {
  bool compareStubLayer( L1TkStub_PixelDigi_Collection::const_iterator s1, L1TkStub_PixelDigi_Collection::const_iterator s2) {
    return StackedTrackerDetId(s1->getDetId()).iLayer() < StackedTrackerDetId(s2->getDetId()).iLayer();
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
}  


//
// class declaration
//

class TrackTriggerStudy : public edm::EDAnalyzer {

public:
  typedef L1TkTrack_PixelDigi_                                       L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollection;
  typedef L1TkStub_PixelDigi_Collection::const_iterator              L1TkStubIter;
  typedef vector< L1TkStubIter >                                     L1TkStubIters;
  typedef pair < L1TkStubIter, L1TkStubIter >                        L1TkStubPair;
  typedef vector< L1TkStubPair >                                     L1TkStubPairs;

  explicit TrackTriggerStudy(const edm::ParameterSet&);
  ~TrackTriggerStudy();
  
      
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillSimTrackInfo();
  void fillRecTrackInfo();
  void fillElectronInfo(l1extra::L1EmParticleCollection::const_iterator igam);
  void matchTrackEGamma(GlobalPoint epos, double ee, TTIStudy::Electron& electron);
  void matchStubEGamma(GlobalPoint epos, double ee, TTIStudy::Electron& electron);

 
  int getSimTrkId(L1TkStubIter s); 
  int getSimTrkIndex(L1TkStubIter s); 
  int getSimTrkIndex(const unsigned int id);

  double getZ0(L1TkStubIter s1,
	       L1TkStubIter s2);
  double getAngle(L1TkStubIter s1,
		  L1TkStubIter s2,
		  double z0);
  GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
  void getTwoPointTracklets(GlobalPoint epos, double et, map< unsigned int, L1TkStubIters >& matched_stubs, vector< L1TkStubIters >& theCands);
  void getThreePointTracklets(vector< L1TkStubIters >& twoPointCands, vector< L1TkStubIters >& theCands);
  void getFourPointTracklets(vector< L1TkStubIters >& twoPointCands, vector< L1TkStubIters >& threePointCands, vector< L1TkStubIters >& theCands);
  template<typename T> void removeDuplicates(std::vector<T>& vec);
  double getTwoPointPt(GlobalPoint a, GlobalPoint b);
  bool goodTwoPointPhi(L1TkStubIter s1, L1TkStubIter s2) ;
  bool goodTwoPointZ(L1TkStubIter s1, L1TkStubIter s2) ;
  double getdP(GlobalPoint epos, double et, L1TkStubIter s1);
  double getOriginZ(GlobalPoint epos, L1TkStubIter s1);
  double getPhiMiss(double et, L1TkStubIter s1, L1TkStubIter s2) ;
  double getPhiMissUncorrected(double et, L1TkStubIter s1, L1TkStubIter s2) ;
  double getPhiMissTruth(L1TkStubIter s1, L1TkStubIter s2) ;
  double getZMiss(GlobalPoint epos, L1TkStubIter s1, L1TkStubIter s2) ;
  L1TkStubIters  getStubsMatchedTo(int simTrackId, unsigned int layer);
  double getTrkEGamdP(GlobalPoint epos, L1TkTrackCollection::const_iterator trk);
  double getTrkEGamOriginZ(GlobalPoint epo, L1TkTrackCollection::const_iterator trk);
  double getTrkEGamDeltaR(GlobalPoint epos, L1TkTrackCollection::const_iterator trk);
  double getTrkEGamDeltaRCorr(GlobalPoint epos, L1TkStubIter s1, L1TkStubIter s2, L1TkTrackCollection::const_iterator trk);
  bool matchedSimTrack(L1TkTrackCollection::const_iterator trk);
  void findMatchedPairs(GlobalPoint epos, double et, L1TkStubIters& stub_vec1, L1TkStubIters& stub_vec2, vector<L1TkStubIters>& candidates); 
  unsigned int findCompatibleStubs(L1TkStubIters & stub_vec1, L1TkStubIters& stub_vec2, L1TkStubIters& candidate);
  unsigned int checkMatchedCandidates( L1TkStubIters & cand1, L1TkStubIters& cand2);
  bool isUniqueCandidate( L1TkStubIters& cand1, L1TkStubIters& cand2);
  void addTrackletCandidate(L1TkStubIters& stubs,  L1TkStubIter& a_stub);
  int getMatchedGenIndex(GlobalPoint epos, double eet);
 
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

  TTIStudy::Event* eventBr_;
  std::vector<TTIStudy::Electron>* electronsBr_; 
  std::vector<TTIStudy::SimTrack>* simTracksBr_;  
  std::vector<TTIStudy::Track>* tracksBr_;  

  // ntuples:

  TTree* tree;


  // parameters that can be set in cfg file:
  unsigned int mcType_ ;
  edm::InputTag egammaSrc_;
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

  edm::Handle<reco::BeamSpot> theBeamSpot;
  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
  const StackedTrackerGeometry*         theStackedTracker;
  ESHandle<TrackerGeometry> theTrkGeomHandle;

  Handle< L1TkStub_PixelDigi_Collection > theGStubs, theGStubs2;
  std::vector< PCaloHit > theEcalSimHits;

  edm::SimTrackContainer simTracks_;
  edm::SimVertexContainer simVertices_;
  L1TkTrackCollection l1TkTracks_;
};

//
// constructors and destructor
//
TrackTriggerStudy::TrackTriggerStudy(const edm::ParameterSet& iConfig)
{
  
  mcType_    = iConfig.getParameter<unsigned int>("mcType");
  egammaSrc_ = iConfig.getParameter<edm::InputTag>("egammaSrc");
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
    std::cout << " Layer radius " << (*it) << std::endl; 
  }
}

TrackTriggerStudy::~TrackTriggerStudy()
{
  delete eventBr_;
  delete simTracksBr_; 
  delete tracksBr_; 
  delete electronsBr_; 
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackTriggerStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  simTracksBr_->clear(); 
  tracksBr_->clear(); 
  electronsBr_->clear(); 

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
  iEvent.getByLabel("BeamSpotFromSim", "BeamSpot", theBeamSpot);
  const GlobalPoint beamSpot(theBeamSpot->x0(),theBeamSpot->y0(),theBeamSpot->z0());
  //  const GlobalPoint beamSpot(0,0,0);

  iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  magnetStrength = theMagField.product()->inTesla(GlobalPoint(0,0,0)).z();
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
    
  iEvent.getByLabel("L1TkStubsFromPixelDigis","StubsPass", theGStubs);
  std::cout << "Found " << theGStubs->size() << " global stubs from Pixel Digis." << std::endl;
  nstub_ = theGStubs->size();
  ////////////////////////////////////////////////////////////////
  // L1Tracks From Track Trigger Hits
  ////////////////////////////////////////////////////////////////

  edm::Handle<L1TkTrackCollection> l1TkTrackHandle;
  iEvent.getByLabel("L1Tracks", "Level1TkTracks", l1TkTrackHandle);
  l1TkTracks_ = (*l1TkTrackHandle.product());
  std::cout << "Found " << l1TkTracks_.size() << " L1 Tracks" << std::endl;

  fillSimTrackInfo();
  fillRecTrackInfo();

  Handle<l1extra::L1EmParticleCollection> CandHandle;
  iEvent.getByLabel(egammaSrc_,CandHandle);
  l1extra::L1EmParticleCollection eGammaCollection = (*CandHandle.product());
  sort(eGammaCollection.begin(), eGammaCollection.end(), TTIStudy::EtComparator());

  std::cout << " EGammas " << eGammaCollection.size() << std::endl;

  for ( l1extra::L1EmParticleCollection::const_iterator EGammaIter = eGammaCollection.begin(); 
	EGammaIter != eGammaCollection.end();  ++EGammaIter) {

    fillElectronInfo(EGammaIter);
   
  }
  tree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void 
TrackTriggerStudy::beginJob()
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Framework handles for the EVENTSETUP tracker geometry, L1 stack geometry, etc...
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("ElectronTriggerInfo","Electron Trigger Information");
  
  eventBr_ = new TTIStudy::Event();
  tree->Branch("Event", "TTIStudy::Event", &eventBr_, 32000, 2);

  electronsBr_ = new std::vector<TTIStudy::Electron>(); 
  tree->Branch("Electron", "std::vector<TTIStudy::Electron>", &electronsBr_);

  simTracksBr_ = new std::vector<TTIStudy::SimTrack>(); 
  tree->Branch("SimTrack", "std::vector<TTIStudy::SimTrack>", &simTracksBr_);
  
  tracksBr_ = new std::vector<TTIStudy::Track>(); 
  tree->Branch("Track", "std::vector<TTIStudy::Track>", &tracksBr_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackTriggerStudy::endJob() {
}
int TrackTriggerStudy::getSimTrkId(L1TkStub_PixelDigi_Collection::const_iterator s) {
  
  int id = -5;
  if (s->isGenuine()) {
    edm::Ptr< SimTrack > iSimtk = s->getSimTrackPtr();
    if (!iSimtk.isNull()) id = iSimtk->trackId();
  }
  return id;
}
int TrackTriggerStudy::getSimTrkIndex(L1TkStub_PixelDigi_Collection::const_iterator s) {
  int index = -5;
  if (s->isGenuine()) {
    edm::Ptr< SimTrack > iSimtk = s->getSimTrackPtr();
    if (!iSimtk.isNull()) index = getSimTrkIndex(iSimtk->trackId());
  }
  return index;
}



double TrackTriggerStudy::getZ0(L1TkStubIter s1,
				L1TkStubIter s2) {
  GlobalPoint hitpos1 = theStackedTracker->findGlobalPosition(&(*s1));
  GlobalPoint hitpos2 = theStackedTracker->findGlobalPosition(&(*s2));
  
  double r1x = hitpos1.perp();
  double r2x = hitpos2.perp();
  double z1x = hitpos1.z();
  double z2x = hitpos2.z();
  
  return (z1x - r1x*(z2x-z1x)/(r2x-r1x));
}

double TrackTriggerStudy::getAngle(L1TkStubIter s1,
				   L1TkStubIter s2,
				   double z0) {

  GlobalPoint primary(theStackedTracker->findGlobalPosition(&(*s1)).x(),
		      theStackedTracker->findGlobalPosition(&(*s1)).y(),
		      theStackedTracker->findGlobalPosition(&(*s1)).z() - z0);
  GlobalPoint secondary(theStackedTracker->findGlobalPosition(&(*s2)).x(),
			theStackedTracker->findGlobalPosition(&(*s2)).y(),
			theStackedTracker->findGlobalPosition(&(*s2)).z() - z0);

  double angle = (primary.x()*secondary.x()+
		  primary.y()*secondary.y()+
		  primary.z()*secondary.z())/
    primary.mag()/secondary.mag();

  return angle;
}

TrackTriggerStudy::L1TkStubIters TrackTriggerStudy::getStubsMatchedTo(int simTrackId, unsigned int layer) {

  TrackTriggerStudy::L1TkStubIters matchedStubs;

  L1TkStubIter stubIter;
  for (  stubIter = theGStubs->begin(); stubIter != theGStubs->end() ; ++stubIter ) {
    if(StackedTrackerDetId(stubIter->getDetId()).iLayer()+1 == layer && getSimTrkId(stubIter) == simTrackId)
      matchedStubs.push_back(stubIter);    
  }

  return matchedStubs;
}

GlobalPoint TrackTriggerStudy::getCalorimeterPosition(double phi, double eta, double e) {
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
			    rperp * rperp ) * abs( eta ) / eta;  
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}

int TrackTriggerStudy::getMatchedGenIndex(GlobalPoint epos, double eet){
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
void TrackTriggerStudy::getTwoPointTracklets(GlobalPoint epos, double et, map< unsigned int, L1TkStubIters >& matched_stubs, 
				       vector<L1TkStubIters>& theCands) {
  for (map< unsigned int, L1TkStubIters >::const_iterator it1 = matched_stubs.begin(); it1 != matched_stubs.end(); it1++) {

    L1TkStubIters it_stub1 = it1->second;
    unsigned iLayer = it1->first;
     
    for (map< unsigned int, L1TkStubIters >::const_iterator it2 = matched_stubs.begin(); it2 != matched_stubs.end(); it2++) {
      unsigned jLayer = it2->first;
      if (jLayer <= iLayer) continue;
      L1TkStubIters it_stub2 = it2->second;
      findMatchedPairs(epos, et, it_stub1, it_stub2, theCands); 
    }
  }
  removeDuplicates(theCands); 
  if (debugFlag_) {
    int ival = 0;
    for (vector<L1TkStubIters>::const_iterator ic = theCands.begin(); ic != theCands.end(); ic++) {
      ival++;
      //      cout << ival << " Layer1 " << StackedTrackerDetId((*ic)[0]->getDetId()).iLayer() << "  Layer2 " << StackedTrackerDetId((*ic)[1]->getDetId()).iLayer() << std::endl;
    }
  }
}
void TrackTriggerStudy::findMatchedPairs(GlobalPoint epos, double et, L1TkStubIters & stub_vec1, L1TkStubIters& stub_vec2, vector<L1TkStubIters >& candidates) {

  for (L1TkStubIters::const_iterator stub1 = stub_vec1.begin(); stub1 != stub_vec1.end(); stub1++) {
    for (L1TkStubIters::const_iterator stub2 = stub_vec2.begin(); stub2 != stub_vec2.end(); stub2++) {
      unsigned int layer1 = StackedTrackerDetId((*stub1)->getDetId()).iLayer(); 
      unsigned int layer2 = StackedTrackerDetId((*stub2)->getDetId()).iLayer(); 
      if (layer1 == layer2) continue;
      if(!goodTwoPointPhi((*stub1),(*stub2))) continue;
      if(!goodTwoPointZ((*stub1),(*stub2))) continue;
      double zMiss = getZMiss(epos, (*stub1), (*stub2));
      double phiMiss = getPhiMiss(et, (*stub1), (*stub2));
      double zMissMultiplier = 1.0;
      if (geometryOption_ == "LB_6PS") 
	zMissMultiplier = ((radii_[1])/(radii_[layer2]))*(radii_[layer2]-radii_[layer1])/(rcal_-radii_[layer1])*(rcal_-radii_[1])/(radii_[2]-radii_[1]);
      else if (geometryOption_ == "BE" || geometryOption_ == "BE5D") {
	if (TMath::Abs(epos.eta()) > 1.0) zMissMultiplier = 1.0;
	else  zMissMultiplier = (radii_[layer2]-radii_[layer1])/(rcal_-radii_[layer1])*(rcal_-radii_[1])/(radii_[2]-radii_[1]);
      }
      double zMissCut = zMissMultiplier*(stubEGZMissCut_+0.07*
					 (1.0/(1.0-cos(2*atan(exp(-1.0*TMath::Abs(epos.eta())))))));
      L1TkStubIters aCand;
      if (fabs(phiMiss)< stubEGPhiMissCut_) {
	//	std::cout << " TrackTriggerStudy::findMatchedPairs : phiMiss and ZMiss "
	//		  << fabs(phiMiss) << "  " 
	//		  << fabs(zMiss) <<"  " 
	//		  << " Zmiss cut " << zMissCut << " Layers " << layer1 << " " << layer2 << std::endl;
      }
      if(fabs(phiMiss)< stubEGPhiMissCut_ && fabs(zMiss) < zMissCut) {
        aCand.push_back((*stub1)); 
	aCand.push_back((*stub2)); 
        sort( aCand.begin(), aCand.end(), TTIStudy::compareStubLayer); 
	if (aCand.size() > 0 && layer1 != layer2){
	  sort( aCand.begin(), aCand.end(), TTIStudy::compareStubLayer);                                                       
	  candidates.push_back(aCand);                                                                                          
	}
      } 
    } 
  }
}
void TrackTriggerStudy::getThreePointTracklets( vector< L1TkStubIters >& twoPointCands, vector< L1TkStubIters >& theCands) {
  
  for(vector< L1TkStubIters >::const_iterator it1 = twoPointCands.begin(); it1 != twoPointCands.end(); it1++) {
    L1TkStubIters stub_vec1 = (*it1);
    for(vector< L1TkStubIters >::const_iterator it2 = it1+1; it2 != twoPointCands.end(); it2++) {
      L1TkStubIters stub_vec2 = (*it2);
      if (it1 == it2) continue; 
      L1TkStubIters aCand;
      unsigned int ncomp  = findCompatibleStubs(stub_vec1, stub_vec2, aCand); 
      if (ncomp == 3) {
	vector< L1TkStubIters >::iterator ipos = std::find(theCands.begin(),theCands.end(),aCand);
        if (ipos == theCands.end()) theCands.push_back(aCand);                                                                                                                       
	/*	std::cout << " Size of theCands " << theCands.size() << std::endl; */
      }
    }
  }
  removeDuplicates(theCands);
  if (debugFlag_) {
    int ival = 0; 
    for (vector<L1TkStubIters>::const_iterator ic = theCands.begin(); ic != theCands.end(); ic++) {
      ival++;
      cout << ival << " Layer1 " << StackedTrackerDetId((*ic)[0]->getDetId()).iLayer() << "  Layer2 " << StackedTrackerDetId((*ic)[1]->getDetId()).iLayer() << " layer3 " << StackedTrackerDetId((*ic)[2]->getDetId()).iLayer() << std::endl;
    }
  }
}

void TrackTriggerStudy::getFourPointTracklets( vector< L1TkStubIters> & twoPointCands, vector< L1TkStubIters >& threePointCands, vector< L1TkStubIters >& theCands) {
  // look for overlaps between 2 and 3 point candidates
  for(vector< L1TkStubIters>::const_iterator it2 = twoPointCands.begin(); it2 != twoPointCands.end(); it2++) {
    L1TkStubIters stub_doublet = (*it2);
    for(vector< L1TkStubIters>::const_iterator it3 = threePointCands.begin(); it3 != threePointCands.end(); it3++) {
      L1TkStubIters stub_triplet = (*it3);
      L1TkStubIters aCand;

      unsigned int ncomp  = findCompatibleStubs(stub_doublet, stub_triplet, aCand); 
      if (ncomp == 4) {
	vector< L1TkStubIters >::iterator ipos = std::find(theCands.begin(),theCands.end(),aCand);
	if (ipos == theCands.end()) theCands.push_back(aCand);
      }
    }
  }
  removeDuplicates(theCands);

  if (debugFlag_) {
    int ival = 0; 
    for (vector<L1TkStubIters>::const_iterator ic = theCands.begin(); ic != theCands.end(); ic++) {
      ival++;
      cout << ival << " Layer1 " << StackedTrackerDetId((*ic)[0]->getDetId()).iLayer() << " Layer2 " << StackedTrackerDetId((*ic)[1]->getDetId()).iLayer() << " layer3 " << StackedTrackerDetId((*ic)[2]->getDetId()).iLayer() << " Layer4 " << StackedTrackerDetId((*ic)[3]->getDetId()).iLayer() <<  std::endl;
    }
  }
}

template<typename T>
void TrackTriggerStudy::removeDuplicates(std::vector<T>& vec)
{
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


double TrackTriggerStudy::getTwoPointPt(GlobalPoint a, GlobalPoint b)
{
  double phi_a = a.phi();
  double phi_b = b.phi();
  double r_a = a.perp();
  double r_b = b.perp();
  
  return fabs((r_b-r_a)/(phi_b-phi_a))*(2.998e8*3.8)/(2e11);     
}

bool TrackTriggerStudy::goodTwoPointZ(L1TkStubIter s1,
		     L1TkStubIter s2) {

  double innerPointZ = theStackedTracker->findGlobalPosition(&(*s1)).z();
  double outerPointZ = theStackedTracker->findGlobalPosition(&(*s2)).z();
  double innerPointRadius = theStackedTracker->findGlobalPosition(&(*s1)).perp();
  double outerPointRadius = theStackedTracker->findGlobalPosition(&(*s2)).perp();

  double mIPWidth = 200.0;
  double positiveZBoundary = (mIPWidth - outerPointZ) * (outerPointRadius - innerPointRadius);
  double negativeZBoundary = -(mIPWidth + outerPointZ) * (outerPointRadius - innerPointRadius);
  double multipliedLocation = (innerPointZ - outerPointZ) * outerPointRadius;
  
  
  if( multipliedLocation < positiveZBoundary &&
      multipliedLocation > negativeZBoundary )
    return true;
  return false;
}

bool TrackTriggerStudy::goodTwoPointPhi(L1TkStubIter s1, L1TkStubIter s2) {
  
  double innerPointPhi = theStackedTracker->findGlobalPosition(&(*s1)).phi();
  double outerPointPhi = theStackedTracker->findGlobalPosition(&(*s2)).phi();
  
  double innerPointRadius = theStackedTracker->findGlobalPosition(&(*s1)).perp();
  double outerPointRadius = theStackedTracker->findGlobalPosition(&(*s2)).perp();

  // Rebase the angles in terms of 0-2PI, should
  if ( innerPointPhi < 0.0 ) innerPointPhi += 2.0 * TMath::Pi();
  if ( outerPointPhi < 0.0 ) outerPointPhi += 2.0 * TMath::Pi();

  // Check for seed compatibility given a pt cut
  // Threshold computed from radial location of hits
  double mCompatibilityScalingFactor = 
    (100.0 * 2.0e+9 * 2.0) / (TMath::C() * magnetStrength);
  
  mCompatibilityScalingFactor = 
    1.0 / mCompatibilityScalingFactor;
  
  double deltaPhiThreshold = 
    (outerPointRadius - innerPointRadius) * 
    mCompatibilityScalingFactor;  
  
  // Delta phi computed from hit phi locations
  double deltaPhi = outerPointPhi - innerPointPhi;
  if(deltaPhi<0) deltaPhi = -deltaPhi;
		   
  if(deltaPhi<deltaPhiThreshold)
    return true;
  return false;
}	

double TrackTriggerStudy::getdP(GlobalPoint epos, double eet, L1TkStubIter s1) {
  
  double er = epos.perp();

  double magneticField = 3.8;  
  double phiVsRSlope = -3.00e-3 * magneticField / eet / 2.;
  
  double r1 = theStackedTracker->findGlobalPosition(&(*s1)).perp();
  double phi1 = theStackedTracker->findGlobalPosition(&(*s1)).phi();

  // preselecton variable
  double psi1 = phiDiff(phi1,epos.phi());
  double deltaPsi1 = psi1 - (er-r1)*phiVsRSlope;
  double antiDeltaPsi1 = psi1 - (r1-er)*phiVsRSlope;
  double dP;
  if (fabs(deltaPsi1)<fabs(antiDeltaPsi1)){
    dP = deltaPsi1;
  }else{
    dP = antiDeltaPsi1;
  }
  return dP;
}

double TrackTriggerStudy::getOriginZ(GlobalPoint epos, L1TkStubIter s1) {
	  
  double er = epos.perp();
  double ez = epos.z();
  
  double r1 = theStackedTracker->findGlobalPosition(&(*s1)).perp();
  double z1 = theStackedTracker->findGlobalPosition(&(*s1)).z();

  double originZ1 = (er*z1 - r1*ez)/(er-r1);
  return originZ1;
}

// Extrapolate a L1Track from its vertex to the ECal Cluster position::method
double TrackTriggerStudy::getTrkEGamdP(GlobalPoint epos, L1TkTrackCollection::const_iterator trk) {
  double er = epos.perp();

  // Using track fit curvature
  //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
  double curv = trk->getRInv();
  double x1 = (asin(er*curv/(2.0)));
  double phi1 = phiDiff(trk->getMomentum().phi(), epos.phi());

  double dif1 = phi1 - x1;
  double dif2 = phi1 + x1; 
  if (abs(dif1) < abs(dif2)) return dif1;
  else return dif2; 
}

double TrackTriggerStudy::getTrkEGamOriginZ(GlobalPoint epos, L1TkTrackCollection::const_iterator trk) {
	  
  double er = epos.perp();
  double ez = epos.z();
  
  //  TVector3 trkvtx;
  //  trkvtx.SetXYZ(trk->vx(), trk->vy(), trk->vz());
  double r1 = trk->getVertex().perp();
  double z1 = trk->getVertex().z();

  double originZ1 = (er*z1 - r1*ez)/(er-r1);
  return originZ1;
}

double TrackTriggerStudy::getTrkEGamDeltaR(GlobalPoint epos, L1TkTrackCollection::const_iterator trk){
  double dPhi = phiDiff(epos.phi(), trk->getMomentum().phi());
  //double dPhi = (ephi - trk->getMomentum().phi());
  //if (fabs(dPhi) > 3.14285)
    //dPhi = (trk->getMomentum().phi() - ephi);
  double dEta = (epos.eta() - trk->getMomentum().eta());
  double deltaR = sqrt(dPhi*dPhi + dEta*dEta);
  return deltaR;
}
double TrackTriggerStudy::getTrkEGamDeltaRCorr(GlobalPoint epos, L1TkStubIter s1, L1TkStubIter s2, L1TkTrackCollection::const_iterator trk){
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

double TrackTriggerStudy::getPhiMiss(double eet, L1TkStubIter s1, L1TkStubIter s2) {

  
  double pT = eet;
  double curv = pT*100*.877;
  if (curv == 0) return 999.9;
 	    
  GlobalPoint s1pos(theStackedTracker->findGlobalPosition(&(*s1)).x()-theBeamSpot->x0(),
		    theStackedTracker->findGlobalPosition(&(*s1)).y()-theBeamSpot->y0(),
		    theStackedTracker->findGlobalPosition(&(*s1)).z()-theBeamSpot->z0());

  GlobalPoint s2pos(theStackedTracker->findGlobalPosition(&(*s2)).x()-theBeamSpot->x0(),
		    theStackedTracker->findGlobalPosition(&(*s2)).y()-theBeamSpot->y0(),
		    theStackedTracker->findGlobalPosition(&(*s2)).z()-theBeamSpot->z0());

  double r1 = s1pos.perp();
  double r2 = s2pos.perp();

  double phi1 = s1pos.phi();
  double phi2 = s2pos.phi();

  //Predict phi of hit 2
  double a = (r2-r1)/(2*curv);
  double b = phiDiff(phi2,phi1);
  
  double phiMiss = 0;
  if(fabs(b - a)<fabs(b + a)) phiMiss = b - a;
  if(fabs(b - a)>fabs(b + a)) phiMiss = b + a;
	    
  return phiMiss;

}

// calculate phimiss uncorrected for beamspot offset
double TrackTriggerStudy::getPhiMissUncorrected(double eet, L1TkStubIter s1, L1TkStubIter s2) {

  
  double pT = eet;
  double curv = pT*100*.877;
  if (curv == 0) return 999.9;
	    
  GlobalPoint s1pos(theStackedTracker->findGlobalPosition(&(*s1)).x(),
		    theStackedTracker->findGlobalPosition(&(*s1)).y(),
		    theStackedTracker->findGlobalPosition(&(*s1)).z());

  GlobalPoint s2pos(theStackedTracker->findGlobalPosition(&(*s2)).x(),
		    theStackedTracker->findGlobalPosition(&(*s2)).y(),
		    theStackedTracker->findGlobalPosition(&(*s2)).z());

  double r1 = s1pos.perp();
  double r2 = s2pos.perp();

  double phi1 = s1pos.phi();
  double phi2 = s2pos.phi();

  //Predict phi of hit 2
  double a = (r2-r1)/(2*curv);
  double b = phiDiff(phi2,phi1);
  
  double phiMiss = 0;
  if(fabs(b - a)<fabs(b + a)) phiMiss = b - a;
  if(fabs(b - a)>fabs(b + a)) phiMiss = b + a;
	    
  return phiMiss;

}

double TrackTriggerStudy::getPhiMissTruth(L1TkStubIter s1, L1TkStubIter s2) {
  int trkindx = getSimTrkIndex(s1);
  int trkindx2 = getSimTrkIndex(s2);
  double pT = 0;
  if(trkindx != -5 && trkindx==trkindx2) {       
    pT = simTracks_[trkindx].momentum().pt();
    
    double curv = pT*100*.877;
    if (curv == 0) return 999.9;
    
    GlobalPoint s1pos(theStackedTracker->findGlobalPosition(&(*s1)).x()-theBeamSpot->x0(),
		      theStackedTracker->findGlobalPosition(&(*s1)).y()-theBeamSpot->y0(),
		      theStackedTracker->findGlobalPosition(&(*s1)).z()-theBeamSpot->z0());
    
    GlobalPoint s2pos(theStackedTracker->findGlobalPosition(&(*s2)).x()-theBeamSpot->x0(),
		      theStackedTracker->findGlobalPosition(&(*s2)).y()-theBeamSpot->y0(),
		      theStackedTracker->findGlobalPosition(&(*s2)).z()-theBeamSpot->z0());
    
    double r1 = s1pos.perp();
    double r2 = s2pos.perp();
    
    double phi1 = s1pos.phi();
    double phi2 = s2pos.phi();
    
    //Predict phi of hit 2
    double a = (r2-r1)/(2*curv);
    double b = phiDiff(phi2,phi1);
    
    double phiMiss = 0;
    if(fabs(b - a)<fabs(b + a)) phiMiss = b - a;
    if(fabs(b - a)>fabs(b + a)) phiMiss = b + a;
    
    return phiMiss;
  }
  return -1;

}

double TrackTriggerStudy::getZMiss(GlobalPoint epos, L1TkStubIter s1, L1TkStubIter s2) {

  double er = epos.perp();
  double ez = epos.z();

  double r1 = theStackedTracker->findGlobalPosition(&(*s1)).perp();
  double r2 = theStackedTracker->findGlobalPosition(&(*s2)).perp();

  double z1 = theStackedTracker->findGlobalPosition(&(*s1)).z();
  double z2 = theStackedTracker->findGlobalPosition(&(*s2)).z();
 
  double zMiss = z2 - (r2*(ez-z1)-r1*ez + 
		       er*z1)/(er-r1);
  return zMiss;
}

int TrackTriggerStudy::getSimTrkIndex(unsigned int id) {
  int index = -5;
  for(unsigned int i = 0; i<simTracks_.size(); i++) {
    unsigned int simTrkId = simTracks_[i].trackId();
    if (simTrkId == id) {
      index = i;
    }
  }
  return index;
}

bool TrackTriggerStudy::matchedSimTrack(L1TkTrackCollection::const_iterator trk){
  bool matched = false;
  edm::Ptr< SimTrack > simTrackPtr = trk->getSimTrackPtr();
  if (simTrackPtr.isNull()) return matched; 
  unsigned int simtrackid = simTrackPtr->trackId();
  for(unsigned int i = 0; i<simTracks_.size(); i++) {
    unsigned int Id = simTracks_[i].trackId();
    if (Id == simtrackid) {
      matched = true;
      break;
    }
  }

  return matched;
}
unsigned int TrackTriggerStudy::findCompatibleStubs(L1TkStubIters & stub_vec1, L1TkStubIters& stub_vec2,L1TkStubIters& candidate) {
  /*  std::cout << " Layers in input Vectors "; 
  for (L1TkStubIters::const_iterator ic = stub_vec1.begin(); ic != stub_vec1.end(); ic++) {
    std::cout << "  " <<(*ic)->getDetId().layer() << " " ;
  }
  for (L1TkStubIters::const_iterator ic = stub_vec2.begin(); ic != stub_vec2.end(); ic++) {
    std::cout <<  "  " << (*ic)->getDetId().layer() << " ";
    }*/
  unsigned int nCommon = checkMatchedCandidates(stub_vec1, stub_vec2);
  /*  std::cout << " common " << nCommon ;*/
  if (nCommon == 1) {

   for (L1TkStubIters::const_iterator it1 = stub_vec1.begin(); it1 != stub_vec1.end(); it1++) {                                    
     L1TkStubIter stub1 = (*it1);
     addTrackletCandidate(candidate, stub1);
     for (L1TkStubIters::const_iterator it2 = stub_vec2.begin(); it2 != stub_vec2.end(); it2++) {  
	L1TkStubIter stub2 = (*it2); 
	if (stub1 != stub2) addTrackletCandidate(candidate, stub2);
      }
   }
   sort( candidate.begin(), candidate.end(), TTIStudy::compareStubLayer);
   /*   std::cout << " ==> layers in candidate " ; 
   for (L1TkStubIters::const_iterator ic = candidate.begin(); ic != candidate.end(); ic++) {
     std::cout << (*ic)->getDetId().layer() << " ";
     }*/

  }
  /*  std::cout  << std::endl; */
  return candidate.size();
}
unsigned int TrackTriggerStudy::checkMatchedCandidates( L1TkStubIters & cand1, L1TkStubIters& cand2) {
  //  unsigned int jmatched = 999;
  int nMatched = 0; 
  for(unsigned int i = 0; i < cand1.size(); i++) {
    L1TkStubIter element1 = cand1[i];
    unsigned int L1 = StackedTrackerDetId(element1->getDetId()).iLayer(); 
    for(unsigned int j = 0; j < cand2.size(); j++) {
      L1TkStubIter element2 = cand2[j];
      unsigned int L2 = StackedTrackerDetId(element2->getDetId()).iLayer(); 
      if (element1 == element2 && L1 == L2) {
	//        jmatched = j;
        nMatched++;
        break;
      }
    }
  }
  return nMatched;
}
bool TrackTriggerStudy::isUniqueCandidate(L1TkStubIters& cand1, L1TkStubIters& cand2) {
  if (cand1.size() == cand2.size()) {
    unsigned int nS = 0;     
    for (unsigned int ii = 0; ii < cand1.size(); ii++) {
      if (cand1[ii] == cand2[ii] && StackedTrackerDetId(cand1[ii]->getDetId()).iLayer() == StackedTrackerDetId(cand2[ii]->getDetId()).iLayer()) nS++; 
    }
    if (nS == cand1.size()) return false;
  } else std::cout << " Trying to match wrong candidates " << std::endl;
  return true;
}
void TrackTriggerStudy::addTrackletCandidate(L1TkStubIters& stubs,  L1TkStubIter& a_stub) {

  bool new_cand = true;
  unsigned int Layer_a = StackedTrackerDetId(a_stub->getDetId()).iLayer();
  for (L1TkStubIters::const_iterator it = stubs.begin(); it != stubs.end(); it++) {
    L1TkStubIter i_stub = (*it);
    unsigned int Layer_i = StackedTrackerDetId(i_stub->getDetId()).iLayer();
    if (a_stub == i_stub || Layer_a == Layer_i) {
      new_cand = false;
      break;
    }
  }
  if (new_cand) stubs.push_back(a_stub);
}
//
// -- Fill SimTrack Information
//
void TrackTriggerStudy::fillSimTrackInfo() {
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
void TrackTriggerStudy::fillRecTrackInfo() {
  for(L1TkTrackCollection::const_iterator trk = l1TkTracks_.begin(); trk != l1TkTracks_.end() ; ++trk){

    unsigned int PDGId = 0;
    int VtxId = -1;

    edm::Ptr< SimTrack > simTrackPtr = trk->getSimTrackPtr();
    if (!simTrackPtr.isNull()) {
      PDGId = simTrackPtr->type();
      VtxId = simTrackPtr->vertIndex();
    }
    TTIStudy::Track recTk;
    // Track Pt = 0.3* B * radius
    recTk.pt        = trk->getMomentum().perp();
    recTk.eta       = trk->getMomentum().eta();
    recTk.phi       = trk->getMomentum().phi();
    recTk.chiSquare = trk->getChi2Red();
    //    float Curvature = 0.3 * magnetStrength * trk->getCharge()/trk->getMomentum().perp();
    double Curvature = trk->getRInv();
    recTk.curvature = Curvature;
    recTk.vertexZ   = trk->getVertex().z();
    recTk.vertexEta = trk->getVertex().eta();
    recTk.vertexPhi = trk->getVertex().phi();
    recTk.nStub     = trk->getStubPtrs().size();
    recTk.pdgId     = PDGId;
    recTk.vertexId  = VtxId;
    recTk.matchedSimTrack = matchedSimTrack(trk);

    tracksBr_->push_back(recTk);
  } 
}
//
// -- Fill EGamma Properties
//
void TrackTriggerStudy::fillElectronInfo(l1extra::L1EmParticleCollection::const_iterator igam) {

  bool selected = true;
  float e_ele   = igam->energy();
  float eta_ele = igam->eta();
  float et_ele = 0;
  if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
  else et_ele = -1.0;     

  if (et_ele <= egETThreshold_) selected = false; 
    
  /*  if (geometryOption_ == "LB_6PS") {
    if (et_ele <= egETThreshold_) selected = false;
  }  else {
    if(et_ele <= egETThreshold_ || abs(eta_ele) >=1.0) selected = false;
    }*/
  if (!selected) return;
  
  TTIStudy::Electron elec; 
  
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
  
  if (debugFlag_) std::cout << " Et, E, Eta and Phi " << et_ele << " " << e_ele << " " <<  eta_ele << " " << epos.phi() << std::endl;

  matchStubEGamma(epos, e_ele, elec);    
  matchTrackEGamma(epos, e_ele, elec);        
  
  electronsBr_->push_back(elec);
}  
//
// -- Match Track and EGamma Clusters
//
void TrackTriggerStudy::matchTrackEGamma(GlobalPoint epos, double ee, TTIStudy::Electron& elec) {
  
  int ntrackB_ = 0; // Track with Pt > 5 GeV
  int ntrackC_ = 0; // calo+track, dP + deltaR_0.05 + No E/P
  int ntrackD_ = 0; // calo+track, dP + deltaR_0.1 + No E/P
  int ntrackE_ = 0; // calo+track, dP + deltaR_0.15 + No E/P
  int ntrackF_ = 0; // calo+track, dP + deltaR_0.05 + E/P
  int ntrackG_ = 0; // calo+track, dP + deltaR_0.1 + E/P
  int ntrackH_ = 0; // calo+track, dP + deltaR_0.15 + E/P

  for(L1TkTrackCollection::const_iterator trk = l1TkTracks_.begin(); trk != l1TkTracks_.end() ; ++trk){
    if (trk->getMomentum().perp() > 5 && trk->getStubPtrs().size() >= 4 ){
      ntrackB_++;
      float trkEGamDp = getTrkEGamdP(epos, trk);
      float trkEGamDr = getTrkEGamDeltaR(epos, trk);
      float trkEGamEtoP = ee/trk->getMomentum().mag(); 
      //Here we apply the cuts to match Track and EGamma
      if (abs(trkEGamDp) < trkEGdPhiCut_) {
	if (trkEGamDr < trkEGdRCut_[0])  ntrackC_++;
	
	if (trkEGamDr < trkEGdRCut_[1])  ntrackD_++;          
	  
	if (trkEGamDr < trkEGdRCut_[2])  ntrackE_++;

	if ( (trkEGamDr < trkEGdRCut_[0]) && 
	     (trkEGamEtoP > trkEGEtoPCut_[0] && trkEGamEtoP < trkEGEtoPCut_[1]) ) ntrackF_++;
	if ( (trkEGamDr < trkEGdRCut_[1]) && 
	     (trkEGamEtoP > trkEGEtoPCut_[0] && trkEGamEtoP < trkEGEtoPCut_[1]) ) ntrackG_++;
	if ( (trkEGamDr < trkEGdRCut_[2]) && 
	     (trkEGamEtoP > trkEGEtoPCut_[0] && trkEGamEtoP < trkEGEtoPCut_[1]) ) ntrackH_++;
      }
    }
  }
  elec.matchedTrackCounts.push_back(ntrackB_);
  elec.matchedTrackCounts.push_back(ntrackC_);
  elec.matchedTrackCounts.push_back(ntrackD_);
  elec.matchedTrackCounts.push_back(ntrackE_);
  elec.matchedTrackCounts.push_back(ntrackF_);
  elec.matchedTrackCounts.push_back(ntrackG_);
  elec.matchedTrackCounts.push_back(ntrackH_);
}
//
// -- Match Track and EGamma Clusters
//
void TrackTriggerStudy::matchStubEGamma(GlobalPoint epos, double ee, TTIStudy::Electron& electron) {

  int ncandB_ = 0; // 1st+2nd layer (2/2, 4cm sep)
  int ncandC_ = 0; // 2nd+3rd layer (2/2, 12cm sep)
  int ncandD_ = 0; // 2/3 
  int ncandE_ = 0; // 3/3 
  int ncandF_ = 0; // 3/4 
  int ncandG_ = 0; // 4/4 
  int ncandH_ = 0; // 4/6 
  
  int nisoA_ = 0;
  int nisoB_ = 0;
  int nisoC_ = 0;
  int nisoD_ = 0;
  int nisoE_ = 0;


  double et_ele = 0;
  if (cosh(epos.eta()) > 0.0) et_ele = ee/cosh(epos.eta());
  else et_ele = -1.0;      

  // Get stubs consistent with EGamma
  std::map< unsigned int, L1TkStubIters > matchedStubMap; 
  for (L1TkStubIter istub = theGStubs->begin(); istub != theGStubs->end(); ++istub)  {
    float stub_pt = theStackedTracker->findRoughPt(magnetStrength,&(*istub)); 
    if (stubPtCut_ > 0.0 && stub_pt <= stubPtCut_) continue;
    unsigned int ilayer = 999999;
    double originZMultiplier = 1.0;
    if (StackedTrackerDetId(istub->getDetId()).isBarrel()) {
      ilayer = StackedTrackerDetId(istub->getDetId()).iLayer();
      originZMultiplier = radii_[ilayer]/radii_[1];
    } else if (StackedTrackerDetId(istub->getDetId()).isEndcap()) {
      if (StackedTrackerDetId(istub->getDetId()).iSide() == 1) ilayer=StackedTrackerDetId(istub->getDetId()).iDisk()+20;
      else if (StackedTrackerDetId(istub->getDetId()).iSide() == 2) ilayer=StackedTrackerDetId(istub->getDetId()).iDisk()+30;
      else ilayer = 999999;
    }
      
    //    if (ilayer > 10) { // making it 10 for BE
    //      //	std::cout << " WRONG DETECTOR in WRONG LAYER " << ilayer << std::endl;
    //      continue;
    //    }

    double dPi = getdP(epos, et_ele, istub);
    double originZi = getOriginZ(epos, istub);
    if ( (fabs(dPi) < stubEGdPhiCut_) && 
	 (fabs(originZi)<originZMultiplier*(stubEGdZCut_+0.75*(1.0/(1.0-cos(2*atan(exp(-1.0*TMath::Abs(epos.eta())))))))) ) {
      std::map< unsigned int, L1TkStubIters >::iterator iPos = matchedStubMap.find(ilayer);
      if (iPos != matchedStubMap.end()) iPos->second.push_back(istub);
      else {
	vector < L1TkStubIter > stub_vec;
	stub_vec.push_back(istub);
	matchedStubMap.insert(make_pair(ilayer, stub_vec));
      }         
    }
  }
  if (debugFlag_) {
    for (std::map< unsigned int, L1TkStubIters >::const_iterator im = matchedStubMap.begin(); im != matchedStubMap.end(); im++) {
      std::cout << " Layer " << im->first << " matched stubs " << im->second.size() << std::endl;
    }   
  }
  // Get two-point tracklets consistent with this EGamma
  // and loop over them
  vector <L1TkStubIters>  twoPointCands;
  getTwoPointTracklets(epos, et_ele, matchedStubMap, twoPointCands);
  cout<<"Found "<<twoPointCands.size()<<" matched 2-point tracklets."<<endl;
  for(vector <L1TkStubIters>::const_iterator ip = twoPointCands.begin(); ip != twoPointCands.end(); ip++) {

    L1TkStubIter stub1 = (*ip)[0];
    L1TkStubIter stub2 = (*ip)[1]; 
    unsigned int layer1 = StackedTrackerDetId(stub1->getDetId()).iLayer();
    unsigned int layer2 = StackedTrackerDetId(stub2->getDetId()).iLayer();
    int trk_index1 = getSimTrkIndex(stub1);
    int trk_index2 = getSimTrkIndex(stub2);
    
    TTIStudy::Tracklet tracklet;
    tracklet.deltaPhi1       = getdP(epos, et_ele,stub1);
    tracklet.zIntercept1     = getOriginZ(epos,stub1);
    tracklet.layer1          = layer1;
    tracklet.phi1            = theStackedTracker->findGlobalPosition(&(*stub1)).phi();
    tracklet.eta1            = theStackedTracker->findGlobalPosition(&(*stub1)).eta();
    tracklet.x1              = theStackedTracker->findGlobalPosition(&(*stub1)).x();
    tracklet.y1              = theStackedTracker->findGlobalPosition(&(*stub1)).y();
    tracklet.z1              = theStackedTracker->findGlobalPosition(&(*stub1)).z();
    tracklet.r1              = theStackedTracker->findGlobalPosition(&(*stub1)).perp();
    tracklet.pt1             = theStackedTracker->findRoughPt(magnetStrength,&(*stub1));
    tracklet.trackIndex1     = trk_index1;
    if(trk_index1 != -5)  { // not combinatorial 
      tracklet.particleId1 = simTracks_[trk_index1].type();
      tracklet.truePt1     = simTracks_[trk_index1].momentum().pt();
    } 
    tracklet.deltaPhi2       = getdP(epos, et_ele,stub2);
    tracklet.zIntercept2     = getOriginZ(epos,stub2);
    tracklet.layer2          = layer2;
    tracklet.phi2            = theStackedTracker->findGlobalPosition(&(*stub2)).phi();
    tracklet.eta2            = theStackedTracker->findGlobalPosition(&(*stub2)).eta();
    tracklet.x2              = theStackedTracker->findGlobalPosition(&(*stub2)).x();
    tracklet.y2              = theStackedTracker->findGlobalPosition(&(*stub2)).y();
    tracklet.z2              = theStackedTracker->findGlobalPosition(&(*stub2)).z();
    tracklet.r2              = theStackedTracker->findGlobalPosition(&(*stub2)).perp();
    tracklet.pt2             = theStackedTracker->findRoughPt(magnetStrength,&(*stub2));
    tracklet.trackIndex2     = trk_index2;
    if(trk_index2 != -5)  { // not combinatorial 
      tracklet.particleId2   = simTracks_[trk_index2].type();
      tracklet.truePt2       = simTracks_[trk_index2].momentum().pt();
    }
    tracklet.phiMiss = getPhiMiss(et_ele, stub1, stub2);
    tracklet.zMiss   = getZMiss(epos,stub1,stub2);
    GlobalPoint a = theStackedTracker->findGlobalPosition(&(*stub1));
    GlobalPoint b = theStackedTracker->findGlobalPosition(&(*stub2));
    tracklet.twoPointPt = getTwoPointPt(a,b);
    tracklet.twoPointZIntercept = getZ0(stub1, stub2);
    
    if(layer1 == 1 && layer2 == 2) ncandB_++;
    
    if(layer1 == 2 && layer2 == 3) ncandC_++;
    
    if(layer1 < 4  && layer2 < 4) ncandD_++;

    //Here we are doing the tracker-based isolation study
    if(layer1 < 4  && layer2 < 4){
      double zpoint = getZ0(stub1, stub2);
      double iso0 = 0;
      double iso1 = 0;
      double iso2 = 0;
      double iso3 = 0;
      double iso4 = 0;
      for(L1TkTrackCollection::const_iterator trk = l1TkTracks_.begin(); trk != l1TkTracks_.end() ; ++trk){
         double DR = getTrkEGamDeltaRCorr (epos, stub1, stub2, trk);
         if (fabs(zpoint - trk->getVertex().z()) < 0.4 && DR > 0.08 && DR < 0.5) iso0 += trk->getMomentum().perp()/et_ele;
         if (fabs(zpoint - trk->getVertex().z()) < 0.6 && DR > 0.08 && DR < 0.5) iso1 += trk->getMomentum().perp()/et_ele;
         if (fabs(zpoint - trk->getVertex().z()) < 0.8 && DR > 0.08 && DR < 0.5) iso2 += trk->getMomentum().perp()/et_ele;
         if (fabs(zpoint - trk->getVertex().z()) < 1.0 && DR > 0.08 && DR < 0.5) iso3 += trk->getMomentum().perp()/et_ele;
         if (fabs(zpoint - trk->getVertex().z()) < 1.2 && DR > 0.08 && DR < 0.5) iso4 += trk->getMomentum().perp()/et_ele;
      }
      //Here we take the minimum isolation value considering all two point candidiates       
      if (iso0 < 0.15) nisoA_ += 1; 
      if (iso1 < 0.15) nisoB_ += 1; 
      if (iso2 < 0.15) nisoC_ += 1; 
      if (iso3 < 0.15) nisoD_ += 1; 
      if (iso4 < 0.15) nisoE_ += 1; 
    }
    electron.matchedTracklets.push_back(tracklet);
  }
  // get three point tracklets
  vector< L1TkStubIters > threePointCands;
  getThreePointTracklets( twoPointCands, threePointCands );
  cout<<"Found "<<threePointCands.size()<<" matched 3-point tracklets."<<endl;
  for(int threePointIndex = 0; threePointIndex<(int)threePointCands.size(); threePointIndex++) {
    unsigned int layer0 = StackedTrackerDetId(threePointCands[threePointIndex][0]->getDetId()).iLayer();
    unsigned int layer1 = StackedTrackerDetId(threePointCands[threePointIndex][1]->getDetId()).iLayer();
    unsigned int layer2 = StackedTrackerDetId(threePointCands[threePointIndex][2]->getDetId()).iLayer();
    if(layer0 == 1 &&  layer1 == 2 && layer2 == 3)  ncandE_++;
    
    if(layer0 < 5 && layer1 < 5 && layer2 < 5)   ncandF_++;
  }
  
  vector< L1TkStubIters > fourPointCands;
  std::cout << " Size of : Two Point Tracklet " << twoPointCands.size() << " Three Point Tracklet " << threePointCands.size() << std::endl;
  getFourPointTracklets( twoPointCands, threePointCands, fourPointCands);
  cout<<"Found "<<fourPointCands.size()<<" matched 4-point tracklets."<<endl;
  for(int fourPointIndex = 0; fourPointIndex<(int)fourPointCands.size(); fourPointIndex++) {
    unsigned int layer0 = StackedTrackerDetId(fourPointCands[fourPointIndex][0]->getDetId()).iLayer();
    unsigned int layer1 = StackedTrackerDetId(fourPointCands[fourPointIndex][1]->getDetId()).iLayer();
    unsigned int layer2 = StackedTrackerDetId(fourPointCands[fourPointIndex][2]->getDetId()).iLayer();
    unsigned int layer3 = StackedTrackerDetId(fourPointCands[fourPointIndex][3]->getDetId()).iLayer();
    
    if(layer0 == 1 && layer1 == 2 && layer2 == 3 && layer3 == 4) ncandG_++;
    ncandH_++;
  }
  electron.matchedTrackletCounts.push_back(ncandB_);   
  electron.matchedTrackletCounts.push_back(ncandC_);   
  electron.matchedTrackletCounts.push_back(ncandD_);
  electron.matchedTrackletCounts.push_back(ncandE_);      
  electron.matchedTrackletCounts.push_back(ncandF_);   
  electron.matchedTrackletCounts.push_back(ncandG_);          
  electron.matchedTrackletCounts.push_back(ncandH_);          
  std::cout << " size   of matchedTrackletCounts " << electron.matchedTrackletCounts.size() << std::endl;
  electron.matchedIsoTrackletCounts.push_back(nisoA_); 
  electron.matchedIsoTrackletCounts.push_back(nisoB_); 
  electron.matchedIsoTrackletCounts.push_back(nisoC_); 
  electron.matchedIsoTrackletCounts.push_back(nisoD_); 
  electron.matchedIsoTrackletCounts.push_back(nisoE_); 

  std::cout << "B, C, D, E, F, G, H " <<  ncandB_ << " " << ncandC_ << " " << ncandD_ << " " << ncandE_ << " " << ncandF_ << " " << ncandG_ << " " << ncandH_ << std::endl;
} 
//define this as a plug-in
DEFINE_FWK_MODULE(TrackTriggerStudy);
