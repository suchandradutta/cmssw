//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

<<<<<<< HEAD
  
=======

>>>>>>> my_dev
////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/Common/interface/DetSetVector.h"

#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

<<<<<<< HEAD
=======
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>

>>>>>>> my_dev
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"


////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackNtupleMaker : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackNtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~L1TrackNtupleMaker();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
protected:
  
private:
  
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
<<<<<<< HEAD
=======
  int MyProcess;
>>>>>>> my_dev

  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;

<<<<<<< HEAD
  // basic track properties, filled for *all* tracks, regardless of type
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2;
  std::vector<int>*   m_trk_charge;
  std::vector<int>*   m_trk_nstub;

=======
>>>>>>> my_dev
  // sim track properties, filled for *all* sim tracks
  std::vector<float>* m_simtrk_pt;
  std::vector<float>* m_simtrk_eta;
  std::vector<float>* m_simtrk_phi;
  std::vector<float>* m_simtrk_z0;
  std::vector<int>*   m_simtrk_id;
  std::vector<int>*   m_simtrk_type;

  // *sim track* properties, for sim tracks that are matched to a L1 track using simtrackID
  std::vector<float>* m_matchID_simtrk_pt;
  std::vector<float>* m_matchID_simtrk_eta;
  std::vector<float>* m_matchID_simtrk_phi;
  std::vector<float>* m_matchID_simtrk_z0;
  std::vector<int>*   m_matchID_simtrk_id;
  std::vector<int>*   m_matchID_simtrk_type;

  // *L1 track* properties, for sim tracks that are matched to an L1 track using simtrackID
  std::vector<float>* m_matchID_trk_pt;
  std::vector<float>* m_matchID_trk_eta;
  std::vector<float>* m_matchID_trk_phi;
  std::vector<float>* m_matchID_trk_z0;
  std::vector<float>* m_matchID_trk_chi2; 
<<<<<<< HEAD
  std::vector<int>*   m_matchID_trk_charge;
=======
>>>>>>> my_dev
  std::vector<int>*   m_matchID_trk_nstub;
  std::vector<int>*   m_matchID_trk_nmatch;

  // *sim track* properties, for sim tracks that are matched to an L1 track using dR<0.1
  std::vector<float>* m_matchDR_simtrk_pt;
  std::vector<float>* m_matchDR_simtrk_eta;
  std::vector<float>* m_matchDR_simtrk_phi;
  std::vector<float>* m_matchDR_simtrk_z0;
  std::vector<int>*   m_matchDR_simtrk_id;
  std::vector<int>*   m_matchDR_simtrk_type;

  // *L1 track* properties, for sim tracks that are matched to an L1 track using dR<0.1
  std::vector<float>* m_matchDR_trk_pt;
  std::vector<float>* m_matchDR_trk_eta;
  std::vector<float>* m_matchDR_trk_phi;
  std::vector<float>* m_matchDR_trk_z0;
  std::vector<float>* m_matchDR_trk_chi2; 
<<<<<<< HEAD
  std::vector<int>*   m_matchDR_trk_charge;
=======
>>>>>>> my_dev
  std::vector<int>*   m_matchDR_trk_nstub;
  std::vector<int>*   m_matchDR_trk_nmatch;


};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackNtupleMaker::L1TrackNtupleMaker(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{
<<<<<<< HEAD
=======

  MyProcess = iConfig.getParameter< int >("MyProcess");

>>>>>>> my_dev
}

/////////////
// DESTRUCTOR
L1TrackNtupleMaker::~L1TrackNtupleMaker()
{
}  

//////////
// END JOB
void L1TrackNtupleMaker::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackNtupleMaker::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackNtupleMaker::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackNtupleMaker::beginJob" << endl;


  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;


  // initilize
<<<<<<< HEAD
  m_trk_pt     = new std::vector<float>;
  m_trk_eta    = new std::vector<float>;
  m_trk_phi    = new std::vector<float>;
  m_trk_z0     = new std::vector<float>;
  m_trk_chi2   = new std::vector<float>;
  m_trk_charge = new std::vector<int>;
  m_trk_nstub  = new std::vector<int>;

=======
>>>>>>> my_dev
  m_simtrk_pt   = new std::vector<float>;
  m_simtrk_eta  = new std::vector<float>;
  m_simtrk_phi  = new std::vector<float>;
  m_simtrk_z0   = new std::vector<float>;
  m_simtrk_id   = new std::vector<int>;
  m_simtrk_type = new std::vector<int>;

  m_matchID_simtrk_pt   = new std::vector<float>;
  m_matchID_simtrk_eta  = new std::vector<float>;
  m_matchID_simtrk_phi  = new std::vector<float>;
  m_matchID_simtrk_z0   = new std::vector<float>;
  m_matchID_simtrk_id   = new std::vector<int>;
  m_matchID_simtrk_type = new std::vector<int>;

  m_matchID_trk_pt     = new std::vector<float>;
  m_matchID_trk_eta    = new std::vector<float>;
  m_matchID_trk_phi    = new std::vector<float>;
  m_matchID_trk_z0     = new std::vector<float>;
  m_matchID_trk_chi2   = new std::vector<float>;
<<<<<<< HEAD
  m_matchID_trk_charge = new std::vector<int>;
=======
>>>>>>> my_dev
  m_matchID_trk_nstub  = new std::vector<int>;
  m_matchID_trk_nmatch = new std::vector<int>;

  m_matchDR_simtrk_pt   = new std::vector<float>;
  m_matchDR_simtrk_eta  = new std::vector<float>;
  m_matchDR_simtrk_phi  = new std::vector<float>;
  m_matchDR_simtrk_z0   = new std::vector<float>;
  m_matchDR_simtrk_id   = new std::vector<int>;
  m_matchDR_simtrk_type = new std::vector<int>;

  m_matchDR_trk_pt     = new std::vector<float>;
  m_matchDR_trk_eta    = new std::vector<float>;
  m_matchDR_trk_phi    = new std::vector<float>;
  m_matchDR_trk_z0     = new std::vector<float>;
  m_matchDR_trk_chi2   = new std::vector<float>;
<<<<<<< HEAD
  m_matchDR_trk_charge = new std::vector<int>;
=======
>>>>>>> my_dev
  m_matchDR_trk_nstub  = new std::vector<int>;
  m_matchDR_trk_nmatch = new std::vector<int>;


  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

<<<<<<< HEAD
  eventTree->Branch("trk_pt",    &m_trk_pt);
  eventTree->Branch("trk_eta",   &m_trk_eta);
  eventTree->Branch("trk_phi",   &m_trk_phi);
  eventTree->Branch("trk_z0",    &m_trk_z0);
  eventTree->Branch("trk_chi2",  &m_trk_chi2);
  eventTree->Branch("trk_charge",&m_trk_charge);
  eventTree->Branch("trk_nstub", &m_trk_nstub);

=======
>>>>>>> my_dev
  eventTree->Branch("simtrk_pt",  &m_simtrk_pt);
  eventTree->Branch("simtrk_eta", &m_simtrk_eta);
  eventTree->Branch("simtrk_phi", &m_simtrk_phi);
  eventTree->Branch("simtrk_z0",  &m_simtrk_z0);
  eventTree->Branch("simtrk_id",  &m_simtrk_id);
  eventTree->Branch("simtrk_type",&m_simtrk_type);

  eventTree->Branch("matchID_simtrk_pt",  &m_matchID_simtrk_pt);
  eventTree->Branch("matchID_simtrk_eta", &m_matchID_simtrk_eta);
  eventTree->Branch("matchID_simtrk_phi", &m_matchID_simtrk_phi);
  eventTree->Branch("matchID_simtrk_z0",  &m_matchID_simtrk_z0);
  eventTree->Branch("matchID_simtrk_id",  &m_matchID_simtrk_id);
  eventTree->Branch("matchID_simtrk_type",&m_matchID_simtrk_type);

  eventTree->Branch("matchID_trk_pt",    &m_matchID_trk_pt);
  eventTree->Branch("matchID_trk_eta",   &m_matchID_trk_eta);
  eventTree->Branch("matchID_trk_phi",   &m_matchID_trk_phi);
  eventTree->Branch("matchID_trk_z0",    &m_matchID_trk_z0);
  eventTree->Branch("matchID_trk_chi2",  &m_matchID_trk_chi2);
<<<<<<< HEAD
  eventTree->Branch("matchID_trk_charge",&m_matchID_trk_charge);
=======
>>>>>>> my_dev
  eventTree->Branch("matchID_trk_nstub", &m_matchID_trk_nstub);
  eventTree->Branch("matchID_trk_nmatch",&m_matchID_trk_nmatch);

  eventTree->Branch("matchDR_simtrk_pt",  &m_matchDR_simtrk_pt);
  eventTree->Branch("matchDR_simtrk_eta", &m_matchDR_simtrk_eta);
  eventTree->Branch("matchDR_simtrk_phi", &m_matchDR_simtrk_phi);
  eventTree->Branch("matchDR_simtrk_z0",  &m_matchDR_simtrk_z0);
  eventTree->Branch("matchDR_simtrk_id",  &m_matchDR_simtrk_id);
  eventTree->Branch("matchDR_simtrk_type",&m_matchDR_simtrk_type);

  eventTree->Branch("matchDR_trk_pt",    &m_matchDR_trk_pt);
  eventTree->Branch("matchDR_trk_eta",   &m_matchDR_trk_eta);
  eventTree->Branch("matchDR_trk_phi",   &m_matchDR_trk_phi);
  eventTree->Branch("matchDR_trk_z0",    &m_matchDR_trk_z0);
  eventTree->Branch("matchDR_trk_chi2",  &m_matchDR_trk_chi2);
<<<<<<< HEAD
  eventTree->Branch("matchDR_trk_charge",&m_matchDR_trk_charge);
  eventTree->Branch("matchDR_trk_nstub", &m_matchDR_trk_nstub);
  eventTree->Branch("matchDR_trk_nmatch",&m_matchDR_trk_nmatch);


=======
  eventTree->Branch("matchDR_trk_nstub", &m_matchDR_trk_nstub);
  eventTree->Branch("matchDR_trk_nmatch",&m_matchDR_trk_nmatch);

>>>>>>> my_dev
}

//////////
// ANALYZE
void L1TrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

<<<<<<< HEAD
  cerr << "L1TrackNtupleMaker:  Start in analyze()" << endl;


  // clear variables
  m_trk_pt->clear();
  m_trk_eta->clear();
  m_trk_phi->clear();
  m_trk_z0->clear();
  m_trk_chi2->clear();
  m_trk_charge->clear();
  m_trk_nstub->clear();

=======
  //cout << "L1TrackNtupleMaker:  Start in analyze()" << endl;

  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    return;
  }


  // clear variables
>>>>>>> my_dev
  m_simtrk_pt->clear();
  m_simtrk_eta->clear();
  m_simtrk_phi->clear();
  m_simtrk_z0->clear();
<<<<<<< HEAD
  m_simtrk_id->clear();
  m_simtrk_type->clear();
=======
>>>>>>> my_dev

  m_matchID_simtrk_pt->clear();
  m_matchID_simtrk_eta->clear();
  m_matchID_simtrk_phi->clear();
  m_matchID_simtrk_z0->clear();
  m_matchID_simtrk_id->clear();
  m_matchID_simtrk_type->clear();

  m_matchID_trk_pt->clear();
  m_matchID_trk_eta->clear();
  m_matchID_trk_phi->clear();
  m_matchID_trk_z0->clear();
  m_matchID_trk_chi2->clear();
<<<<<<< HEAD
  m_matchID_trk_charge->clear();
=======
>>>>>>> my_dev
  m_matchID_trk_nstub->clear();
  m_matchID_trk_nmatch->clear();
  
  m_matchDR_simtrk_pt->clear();
  m_matchDR_simtrk_eta->clear();
  m_matchDR_simtrk_phi->clear();
  m_matchDR_simtrk_z0->clear();
  m_matchDR_simtrk_id->clear();
  m_matchDR_simtrk_type->clear();

  m_matchDR_trk_pt->clear();
  m_matchDR_trk_eta->clear();
  m_matchDR_trk_phi->clear();
  m_matchDR_trk_z0->clear();
  m_matchDR_trk_chi2->clear();
<<<<<<< HEAD
  m_matchDR_trk_charge->clear();
  m_matchDR_trk_nstub->clear();
  m_matchDR_trk_nstub->clear();
=======
  m_matchDR_trk_nstub->clear();
  m_matchDR_trk_nmatch->clear();
>>>>>>> my_dev



  //-----------------------------------------------------------------------------------------------
  // retrieve various containers
  //-----------------------------------------------------------------------------------------------

  // get sim tracks & sim vertices
<<<<<<< HEAD
  edm::Handle<edm::SimTrackContainer>   simTrackHandle;
  edm::Handle<edm::SimVertexContainer>  simVtxHandle;
  iEvent.getByLabel( "g4SimHits", simTrackHandle );
  iEvent.getByLabel( "g4SimHits", simVtxHandle );

  // get the L1 tracks
  edm::Handle<L1TkTrack_PixelDigi_Collection> L1TrackHandle;
  iEvent.getByLabel("L1Tracks", "Level1TkTracks", L1TrackHandle);



  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  L1TkTrack_PixelDigi_Collection::const_iterator iterL1Track;
  for (iterL1Track = L1TrackHandle->begin(); iterL1Track != L1TrackHandle->end(); ++iterL1Track) {
        
    float tmp_trk_pt  = iterL1Track->getMomentum().perp();
    float tmp_trk_eta = iterL1Track->getMomentum().eta();
    float tmp_trk_phi = iterL1Track->getMomentum().phi();
    float tmp_trk_z0  = iterL1Track->getVertex().z();
    float tmp_trk_chi2   = iterL1Track->getChi2();
    float tmp_trk_charge = 0;
    if (iterL1Track->getRInv() > 0) tmp_trk_charge = 1.0;
    else if (iterL1Track->getRInv() < 0) tmp_trk_charge = -1.0;

    if (tmp_trk_pt < 2.0) continue;

    m_trk_pt ->push_back(tmp_trk_pt);
    m_trk_eta->push_back(tmp_trk_eta);
    m_trk_phi->push_back(tmp_trk_phi);
    m_trk_z0 ->push_back(tmp_trk_z0);
    m_trk_chi2  ->push_back(tmp_trk_chi2);
    m_trk_charge->push_back(tmp_trk_charge);


    // ----------------------------------------------------------------------------------------------
    // get pointers to stubs associated to the L1 track
    std::vector< edm::Ptr< L1TkStub_PixelDigi_ > > theStubs = iterL1Track->getStubPtrs();

    int tmp_trk_nstub = (int) theStubs.size();
    m_trk_nstub->push_back(tmp_trk_nstub);

    /*
    // loop over the stubs
    for (unsigned int i=0; i<(unsigned int)theStubs.size(); i++) {
      bool genuine = theStubs.at(i)->isGenuine();
      if (genuine) {
	bool combine = theStubs.at(i)->isCombinatoric();
	int type = theStubs.at(i)->findType();
	unsigned int simid = theStubs.at(i)->findSimTrackId();
      }
    }
    */

  } //end loop L1tracks

=======
  edm::Handle<edm::SimTrackContainer>  simTrackHandle;
  edm::Handle<edm::SimVertexContainer> simVtxHandle;
  iEvent.getByLabel( "g4SimHits", simTrackHandle );
  iEvent.getByLabel( "g4SimHits", simVtxHandle );

  // sim hits for electrons
  edm::Handle<edm::PSimHitContainer> simHitslowHandle;
  edm::Handle<edm::PSimHitContainer> simHitslowHandleEC;
  iEvent.getByLabel( "g4SimHits", "TrackerHitsPixelBarrelLowTof", simHitslowHandle);
  iEvent.getByLabel( "g4SimHits", "TrackerHitsPixelEndcapLowTof", simHitslowHandleEC);

  // L1 tracks & stubs
  edm::Handle<L1TkStub_PixelDigi_Collection> PixelDigiL1TkStubHandle;
  edm::Handle<L1TkTrack_PixelDigi_Collection> L1TrackHandle;
  iEvent.getByLabel("L1TkStubsFromPixelDigis", "StubsPass", PixelDigiL1TkStubHandle);
  iEvent.getByLabel("L1Tracks", "Level1TkTracks", L1TrackHandle);

  // generator information
  edm::Handle<edm::HepMCProduct> genpHandle;
  iEvent.getByLabel("generator", genpHandle);


  // ----------------------------------------------------------------------------------------------
  // truth z position for electrons/pions where an interaction happens before the first measurement point
  // ----------------------------------------------------------------------------------------------

  const HepMC::GenEvent *myGenEvent = genpHandle->GetEvent();

  float my_zpos = -9999.0;
  for (HepMC::GenEvent::vertex_const_iterator vt = myGenEvent->vertices_begin(); vt != myGenEvent->vertices_end(); ++vt) {
    my_zpos = (float) (*vt)->position().z()/10.0; //must convert from mm to cm here!
    break;
  }
  if (my_zpos < -1000.0) cout << "WARNING, didn't find a true event vertex!" << endl;


  // ----------------------------------------------------------------------------------------------
  // loop over sim hits for electrons & pions
  // ----------------------------------------------------------------------------------------------

  vector<unsigned int> Electron_SimTrackIds;
  vector<unsigned int> Pion_SimTrackIds;
  PSimHitContainer::const_iterator iterSimHits;
  PSimHitContainer::const_iterator iterSimHitsEC;
  
  bool stopHere = false;
  unsigned int iprevious = -1;
  bool pion_stopHere = false;
  unsigned int pion_iprevious = -1;


  if (MyProcess==11 || MyProcess==211) {

    // -----------------------------------------------------------------------------------------------------
    // loop over sim hits in the barrel
    for (iterSimHits = simHitslowHandle->begin(); iterSimHits != simHitslowHandle->end(); ++iterSimHits) {
      unsigned short processType = iterSimHits->processType();
      int type = abs(iterSimHits->particleType());

      //electrons
      if (type == 11 && (processType == 2 || processType == 13)) {
	if (!stopHere) {
	  unsigned int trackId = iterSimHits->trackId();
	  if (trackId != iprevious) Electron_SimTrackIds.push_back(trackId);
	  iprevious = trackId;
	}
      }
      else stopHere = true;
      
      //pions
      if (type == 211 && (processType == 2 || processType == 13)) {
	if (!pion_stopHere) {
	  unsigned int trackId = iterSimHits->trackId();
	  if (trackId != pion_iprevious) Pion_SimTrackIds.push_back(trackId);
	  pion_iprevious = trackId;
	}
      }
      else pion_stopHere = true;
      
    } //end loop over simhits

    
    stopHere = false;
    iprevious = -1;
    pion_stopHere = false;
    pion_iprevious = -1;
  
  
    // -----------------------------------------------------------------------------------------------------
    // loop over sim hits in the endcap
    for (iterSimHitsEC = simHitslowHandleEC->begin(); iterSimHitsEC != simHitslowHandleEC->end(); ++iterSimHitsEC) {
      unsigned short processType = iterSimHitsEC->processType();
      int type = abs(iterSimHitsEC->particleType());
      
      // electrons
      if (type == 11 && (processType == 2 || processType == 13)) {
	if (!stopHere) {
	  unsigned int trackId = iterSimHitsEC->trackId();
	  if (trackId != iprevious) Electron_SimTrackIds.push_back(trackId);
	  iprevious = trackId;
	}
      }
      else stopHere = true;
      
      // pions
      if (type == 211 && (processType == 2 || processType == 13)) {
	if (!pion_stopHere) {
	  unsigned int trackId = iterSimHitsEC->trackId();
	  if (trackId != pion_iprevious) Pion_SimTrackIds.push_back(trackId);
	  pion_iprevious = trackId;
	}
      }
      else pion_stopHere = true;
      
    } //end loop over simhits
    
  }//end MyProcess==11/211
  

  if (MyProcess==11 && Electron_SimTrackIds.size() == 0) return;
  if (MyProcess==211 && Pion_SimTrackIds.size() == 0) return;
  
>>>>>>> my_dev


  // ----------------------------------------------------------------------------------------------
  // loop over sim tracks
  // ----------------------------------------------------------------------------------------------

  SimTrackContainer::const_iterator iterSimTracks;
  SimVertexContainer::const_iterator iterSimVtx;

  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {  
 
    float tmp_simtrk_pt  = iterSimTracks->momentum().pt();
    float tmp_simtrk_eta = iterSimTracks->momentum().eta();
    float tmp_simtrk_phi = iterSimTracks->momentum().phi();
    unsigned int tmp_simtrk_id = iterSimTracks->trackId();
    int tmp_simtrk_type = iterSimTracks->type();

<<<<<<< HEAD
    if (tmp_simtrk_pt < 2.0) continue;

    float matchID_trk_pt  = -9999;
    float matchID_trk_eta = -9999;
    float matchID_trk_phi = -9999;
    float matchID_trk_z0  = -9999;
    float matchID_trk_chi2   = -9999;
    float matchID_trk_charge = -9999;
    int   matchID_trk_nstub  = -9999;
    
    float matchDR_trk_pt  = -9999;
    float matchDR_trk_eta = -9999;
    float matchDR_trk_phi = -9999;
    float matchDR_trk_z0  = -9999;
    float matchDR_trk_chi2   = -9999;
    float matchDR_trk_charge = -9999;
    int   matchDR_trk_nstub  = -9999;

    int nMatchID = 0;
    int nMatchDR = 0;
=======
    if (MyProcess==13 && abs(tmp_simtrk_type)!=13) continue;
    if (MyProcess==11 && abs(tmp_simtrk_type)!=11) continue;
    if (MyProcess==211 && abs(tmp_simtrk_type)!=211) continue;
    if (MyProcess==6 && abs(tmp_simtrk_type)!=211) continue;

    if (tmp_simtrk_pt < 0.2) continue;
    if (fabs(tmp_simtrk_eta) > 2.5) continue;
    if (MyProcess==6 && tmp_simtrk_pt < 2.0) continue;


    // for electrons & pions, consider primary tracks
    if (MyProcess==11 && Electron_SimTrackIds.at(0) != tmp_simtrk_id) continue;
    if (MyProcess==211 && Pion_SimTrackIds.at(0) != tmp_simtrk_id) continue;


    float matchID_trk_pt    = -9999;
    float matchID_trk_eta   = -9999;
    float matchID_trk_phi   = -9999;
    float matchID_trk_z0    = -9999;
    float matchID_trk_chi2  = -9999;
    int   matchID_trk_nstub = -9999;
    
    float matchDR_trk_pt    = -9999;
    float matchDR_trk_eta   = -9999;
    float matchDR_trk_phi   = -9999;
    float matchDR_trk_z0    = -9999;
    float matchDR_trk_chi2  = -9999;
    int   matchDR_trk_nstub = -9999;

    int nMatchID = 0;
    int nMatchDR = 0;
    float dRmin = 5;
>>>>>>> my_dev


    // ----------------------------------------------------------------------------------------------
    // get sim vertex
    int sim_vtxid = iterSimTracks->vertIndex();
    float tmp_simtrk_z0 = -999.9;
<<<<<<< HEAD
=======

>>>>>>> my_dev
    if (sim_vtxid > -1) {
      const SimVertex& theSimVertex = (*simVtxHandle)[sim_vtxid];
      math::XYZTLorentzVectorD trkVtxPos = theSimVertex.position();
      tmp_simtrk_z0 = trkVtxPos.z();
    }
<<<<<<< HEAD
    else {
      cout << "MEBUG: warning, cannot acces sim vertex !?" << endl;
    }
=======
    else cout << "WARNING: cannot acces sim vertex !?" << endl;

    // if the simtrack isn't the primary, the recorded z position will be that of the first interaction 
    // more accurate in this scenario is to instead use the z position of the generated vertex
    if (MyProcess==11 && Electron_SimTrackIds.at(0) != 1) tmp_simtrk_z0 = my_zpos;
    if (MyProcess==211 && Pion_SimTrackIds.at(0) != 1) tmp_simtrk_z0 = my_zpos; 

    if (fabs(tmp_simtrk_z0) > 30.0) continue;
>>>>>>> my_dev

    m_simtrk_pt  ->push_back(tmp_simtrk_pt);
    m_simtrk_eta ->push_back(tmp_simtrk_eta);
    m_simtrk_phi ->push_back(tmp_simtrk_phi);
    m_simtrk_z0  ->push_back(tmp_simtrk_z0);
    m_simtrk_id  ->push_back(tmp_simtrk_id);
    m_simtrk_type->push_back(tmp_simtrk_type);
<<<<<<< HEAD
    

    // ----------------------------------------------------------------------------------------------
    // loop over L1 tracks for match
    L1TkTrack_PixelDigi_Collection::const_iterator iterL1Track;
    for (iterL1Track = L1TrackHandle->begin(); iterL1Track != L1TrackHandle->end(); ++iterL1Track) {
      
      float tmp_trk_pt  = iterL1Track->getMomentum().perp();
      float tmp_trk_eta = iterL1Track->getMomentum().eta();
      float tmp_trk_phi = iterL1Track->getMomentum().phi();
      float tmp_trk_z0  = iterL1Track->getVertex().z();
      float tmp_trk_chi2   = iterL1Track->getChi2();
      unsigned int tmp_trk_simtrackid = iterL1Track->getSimTrackId();
 
      float tmp_trk_charge = 0;
      if (iterL1Track->getRInv() > 0) tmp_trk_charge = 1.0;
      else if (iterL1Track->getRInv() < 0) tmp_trk_charge = -1.0;
      
      if (tmp_trk_pt < 2.0) continue;
      
=======

    
    // ----------------------------------------------------------------------------------------------
    // loop over L1 tracks for match
    // ----------------------------------------------------------------------------------------------

    L1TkTrack_PixelDigi_Collection::const_iterator iterL1Track;
    for (iterL1Track = L1TrackHandle->begin(); iterL1Track != L1TrackHandle->end(); ++iterL1Track) {
      
      float tmp_trk_pt   = iterL1Track->getMomentum().perp();
      float tmp_trk_eta  = iterL1Track->getMomentum().eta();
      float tmp_trk_phi  = iterL1Track->getMomentum().phi();
      float tmp_trk_z0   = iterL1Track->getVertex().z();
      float tmp_trk_chi2 = iterL1Track->getChi2();
      unsigned int tmp_trk_simtrackid = iterL1Track->getSimTrackId();
            

      // ----------------------------------------------------------------------------------------------
>>>>>>> my_dev
      // get pointers to stubs associated to the L1 track
      std::vector< edm::Ptr< L1TkStub_PixelDigi_ > > theStubs = iterL1Track->getStubPtrs();
      int tmp_trk_nstub = (int) theStubs.size();

<<<<<<< HEAD
      
      // ----------------------------------------------------------------------------------------------
      // matching based on dR < 0.1
      float deltaEta = fabs(tmp_trk_eta-tmp_simtrk_eta);
      float deltaPhi = fabs(tmp_trk_phi-tmp_simtrk_phi);
      while (deltaPhi > 3.14159) deltaPhi = fabs(2*3.14159 - deltaPhi);
      float dR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi);
      
      if (dR < 0.1) { //match
	nMatchDR++;
	if (tmp_trk_pt > matchDR_trk_pt) {
=======

      // ----------------------------------------------------------------------------------------------
      // matching based on dR < 0.1
      float deltaEta = tmp_trk_eta-tmp_simtrk_eta;
      float deltaPhi = tmp_trk_phi-tmp_simtrk_phi;
      while (deltaPhi > 3.14159) deltaPhi = fabs(2*3.14159 - deltaPhi);
      float dR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
      
      if (dR < 0.1) { //match
	nMatchDR++;
	if (dR < dRmin) {
>>>>>>> my_dev
	  matchDR_trk_pt  = tmp_trk_pt;
	  matchDR_trk_eta = tmp_trk_eta;
	  matchDR_trk_phi = tmp_trk_phi;
	  matchDR_trk_z0  = tmp_trk_z0;
<<<<<<< HEAD
	  matchDR_trk_chi2   = tmp_trk_chi2;
	  matchDR_trk_charge = tmp_trk_charge;
	  matchDR_trk_nstub  = tmp_trk_nstub;
=======
	  matchDR_trk_chi2  = tmp_trk_chi2;
	  matchDR_trk_nstub = tmp_trk_nstub;
	  dRmin = dR;
>>>>>>> my_dev
	}
      }// end if match dR<0.1
      

      // ----------------------------------------------------------------------------------------------
      // matching based on sim track ID
<<<<<<< HEAD
      if ( (tmp_simtrk_id == tmp_trk_simtrackid)) {
	nMatchID++;
	if (tmp_trk_pt > matchID_trk_pt) {
=======
      bool matchedToPrimary = false;
      if (MyProcess==11) {
	matchedToPrimary = true;
	if (std::find(Electron_SimTrackIds.begin(), Electron_SimTrackIds.end(), tmp_trk_simtrackid) == Electron_SimTrackIds.end()) matchedToPrimary = false;
      }
      else if (MyProcess==211) {
	matchedToPrimary = true;
	if (std::find(Pion_SimTrackIds.begin(), Pion_SimTrackIds.end(), tmp_trk_simtrackid) == Pion_SimTrackIds.end()) matchedToPrimary = false;
      }
      else {
	if (tmp_simtrk_id == tmp_trk_simtrackid) matchedToPrimary = true;
      }

      if (matchedToPrimary) {
	nMatchID++;
	
       	if (nMatchID < 2) {
>>>>>>> my_dev
	  matchID_trk_pt  = tmp_trk_pt;
	  matchID_trk_eta = tmp_trk_eta;
	  matchID_trk_phi = tmp_trk_phi;
	  matchID_trk_z0  = tmp_trk_z0;
<<<<<<< HEAD
	  matchID_trk_chi2   = tmp_trk_chi2;
	  matchID_trk_charge = tmp_trk_charge;
	  matchID_trk_nstub  = tmp_trk_nstub;
	}
      }
      
      
=======
	  matchID_trk_chi2  = tmp_trk_chi2;
	  matchID_trk_nstub = tmp_trk_nstub;
	}
      }
      
>>>>>>> my_dev
    }// end loop L1 tracks


    if (matchID_trk_pt > -1) {
      m_matchID_simtrk_pt  ->push_back(tmp_simtrk_pt);
      m_matchID_simtrk_eta ->push_back(tmp_simtrk_eta);
      m_matchID_simtrk_phi ->push_back(tmp_simtrk_phi);
      m_matchID_simtrk_z0  ->push_back(tmp_simtrk_z0);
<<<<<<< HEAD
      m_matchID_simtrk_type->push_back((int)tmp_simtrk_type);
=======
      m_matchID_simtrk_type->push_back(tmp_simtrk_type);
>>>>>>> my_dev
      m_matchID_simtrk_id  ->push_back(tmp_simtrk_id);
      
      m_matchID_trk_pt ->push_back(matchID_trk_pt);
      m_matchID_trk_eta->push_back(matchID_trk_eta);
      m_matchID_trk_phi->push_back(matchID_trk_phi);
      m_matchID_trk_z0 ->push_back(matchID_trk_z0);
      m_matchID_trk_chi2  ->push_back(matchID_trk_chi2);
<<<<<<< HEAD
      m_matchID_trk_charge->push_back(matchID_trk_charge);
=======
>>>>>>> my_dev
      m_matchID_trk_nstub ->push_back(matchID_trk_nstub);
      m_matchID_trk_nmatch->push_back(nMatchID);
    }

    if (matchDR_trk_pt > -1) {
      m_matchDR_simtrk_pt  ->push_back(tmp_simtrk_pt);
      m_matchDR_simtrk_eta ->push_back(tmp_simtrk_eta);
      m_matchDR_simtrk_phi ->push_back(tmp_simtrk_phi);
      m_matchDR_simtrk_z0  ->push_back(tmp_simtrk_z0);
      m_matchDR_simtrk_id  ->push_back(tmp_simtrk_id);
<<<<<<< HEAD
      m_matchDR_simtrk_type->push_back((int)tmp_simtrk_type);
=======
      m_matchDR_simtrk_type->push_back(tmp_simtrk_type);
>>>>>>> my_dev
      
      m_matchDR_trk_pt ->push_back(matchDR_trk_pt);
      m_matchDR_trk_eta->push_back(matchDR_trk_eta);
      m_matchDR_trk_phi->push_back(matchDR_trk_phi);
      m_matchDR_trk_z0 ->push_back(matchDR_trk_z0);
      m_matchDR_trk_chi2  ->push_back(matchDR_trk_chi2);
<<<<<<< HEAD
      m_matchDR_trk_charge->push_back(matchDR_trk_charge);
=======
>>>>>>> my_dev
      m_matchDR_trk_nstub ->push_back(matchDR_trk_nstub);
      m_matchDR_trk_nmatch->push_back(nMatchDR);
    }
    
<<<<<<< HEAD
    

  } //end loop simtracks
  
=======

  } //end loop simtracks
  


>>>>>>> my_dev
  eventTree->Fill();


} // end of analyze()


///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackNtupleMaker);
