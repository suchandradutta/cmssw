#include "SLHCUpgradeSimulations/ElectronAnalysis/interface/AnalObjects.h"
#include <iostream>

TTIStudy::Event::Event() :
  event(-1),
  run(-1),
  nPileUp(-1) 
{}

// Electron
TTIStudy::Electron::Electron() :
  e(-999.9),
  et(-999.9),
  phi(-999.9),
  eta(-999.9),
  x(-999.9),
  y(-999.9),
  z(-999.9),
  r(-999.9),
  vx(-999.9),
  vy(-999.9),
  vz(-999.9),
  simTkIndx(-1)
{
  matchedTracklets.clear();
  matchedTrackletCounts.clear();
  matchedIsoTrackletCounts.clear();
  matchedTrackCounts.clear();
}
// Information about SimTracks
TTIStudy::SimTrack::SimTrack() :
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  vx(-999.9),
  vy(-999.9),
  vz(-999.9),
  vtxIndx(-1),
  type(-1)
{}
// Information about Reconstruceted Tracks
TTIStudy::Track::Track() :
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  curvature(-999.9),
  chiSquare(-999.9),
  vertexEta(-999.9),
  vertexPhi(-999.9),
  vertexZ(-999.9),
  nStub(-1),
  pdgId(-1),
  vertexId(-1),
  matchedSimTrack(false)
{} 
// Information about Tracklets 
TTIStudy::Tracklet::Tracklet() :
  layer1(-1),
  pt1(-999.9), 
  phi1(-999.9),
  eta1(-999.9),
  x1(-999.9),
  y1(-999.9),
  z1(-999.9),
  r1(-999.9),
  trackIndex1(-5),
  particleId1(-5),
  truePt1(-999.9),
  deltaPhi1(-999.9),
  zIntercept1(-999.9),
  layer2(-1),
  pt2(-999.9), 
  phi2(-999.9),
  eta2(-999.9),
  x2(-999.9),
  y2(-999.9),
  z2(-999.9),
  r2(-999.9),
  trackIndex2(-5),
  particleId2(-5),
  truePt2(-999.9),  
  deltaPhi2(-999.9),
  zIntercept2(-999.9),
  phiMiss(-999.9),
  zMiss(-999.9),
  twoPointPt(-999.9),
  twoPointZIntercept(-999.9)
{}
