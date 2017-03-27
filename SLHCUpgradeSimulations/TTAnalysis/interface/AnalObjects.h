/// ////////////////////////////////////////
/// Objects to be stored in the ntuble
/// ////////////////////////////////////////

#ifndef ANALOBJECTS_H
#define ANALOBJECTS_H

#include "TObject.h"
#include <vector>

namespace TTStudy {
  class Event : public TObject {

  public:
    Event(); 
    ~Event(){}
    
    int event;
    int run;
    int nPileUp;
    float genParticleEt;
    float genParticleEta;
    float genParticlePhi;
    float beamSpotX0;
    float beamSpotY0;
    float beamSpotZ0;
    int   nL1PV;
    float zL1PV;
    float sumL1PV;
    int   nOffPV;
    float zOffPV;
    float sumOffPV;

    ClassDef(Event, 1)
  };

  //  class Tracklet;
  class Tracklet : public TObject {
    
  public:
    Tracklet();
    ~Tracklet(){}

    int layer1;
    float pt1; 
    float phi1;
    float eta1;
    float x1;
    float y1;
    float z1;
    float r1;
    int trackIndex1;
    int particleId1;
    float truePt1;
    float deltaPhi1;
    float zIntercept1;
    
    int layer2;
    float pt2; 
    float phi2;
    float eta2;
    float x2;
    float y2;
    float z2;
    float r2;
    int trackIndex2;
    int particleId2;
    float truePt2;  
    float deltaPhi2;
    float zIntercept2;
    
    float phiMiss;
    float zMiss;
    float twoPointPt;
    float twoPointZIntercept;
    
    ClassDef(Tracklet, 1)
  };
  class Electron : public TObject {
    
  public:
    Electron(); 
    ~Electron(){}
    
    float e;
    float et;
    float phi;
    float eta;
    float x;
    float y;
    float z;
    float r;
    float vx;
    float vy;
    float vz;
    int   simTkIndx;
    float bStrength;
    
    std::vector<Tracklet> matchedTracklets;

    ClassDef(Electron, 1)
  };
  
  class SimTrack : public TObject {
    
  public:
    SimTrack();
    ~SimTrack(){}
    
    float pt;
    float eta;
    float phi;
    float vx;
    float vy;
    float vz;
    int vtxIndx;
    int type;
    
    ClassDef(SimTrack, 1)
  };

  class Track : public TObject {
    
  public:
    Track();
    ~Track() {}
    
    float pt;
    float eta;
    float phi;
    float curvature;
    float chiSquare;
    float chiSquareRed;
    float vertexX;
    float vertexY;
    float vertexZ;

    float pt_p4;
    float eta_p4;
    float phi_p4;
    float curvature_p4;
    float vertexX_p4;
    float vertexY_p4;
    float vertexZ_p4;
    float chiSquare_p4;
    float chiSquareRed_p4;
    
    float ptFromStub;
    int   nStub;
    int   nStub_PS;
    int   nStub_SS;
    int   pdgId;
    int   vertexId;
    bool  matchedSimTrack;

    float d0;
    float z0;    
    float d0Err;
    float z0Err;

    float d0PV;
    float z0PV;
    float d0ErrPV;
    float z0ErrPV;

    ClassDef(Track, 1)
  };

  class GenParticle : public TObject {
    
  public:
    GenParticle();
    ~GenParticle() {}
    
    float eta;
    float phi;
    float p;
    float px;
    float py;
    float pz;
    float pt;
    float energy;
    float vx;
    float vy;
    float vz;
    int status;
    int pdgId;
    int charge;
    int motherIndex;
    std::vector<int> daughterIndices;

    ClassDef(GenParticle, 1)
  };
  class L1Jet : public TObject {
    
  public:
    L1Jet(); 
    ~L1Jet(){}
    
    float pt;
    float phi;
    float eta;
    float et;
    float zvtx;
    float zvtx_tk;
    int   nTk;
    ClassDef(L1Jet, 1)
  };
  class L1Muon : public TObject {
    
  public:
    L1Muon();
    ~L1Muon() {}
    
    float pt;
    float eta;
    float phi;
      int isMip; 
    float isolation;
    float pt_tk;
    float eta_tk;
    float phi_tk;
    float curvature;
    float chiSquare;
    float chiSquareRed;
    float vertexX;
    float vertexY;
    float vertexZ;

    float d0;
    float z0;    
    float d0Err;
    float z0Err;

    ClassDef(L1Muon, 1)
  };
}
#endif
