/// ////////////////////////////////////////
/// Objects to be stored in the ntuble
/// ////////////////////////////////////////

#ifndef ANALOBJECTS_H
#define ANAlOBJECTS_H

#include "TObject.h"
#include <vector>

namespace TTIStudy {
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
    float ptFromStub;
    float eta;
    float phi;
    float curvature;
    float chiSquare;
    float vertexX;
    float vertexY;
    float vertexZ;
    int   nStub;
    int   pdgId;
    int   vertexId;
    bool  matchedSimTrack;
    
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
}
#endif
