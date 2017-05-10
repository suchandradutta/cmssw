#ifndef L1TkTrigger_L1JetParticle_h
#define L1TkTrigger_L1JetParticle_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1TkJetParticle
// 

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"



namespace l1extra {
  
  class L1TkJetParticle : public reco::LeafCandidate
    {
      
    public:
      
      typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
      typedef std::vector< L1TTTrackType > L1TTTrackCollection;
      
      L1TkJetParticle();
      
      L1TkJetParticle( const LorentzVector& p4,
		       const edm::Ref< L1JetParticleCollection >& jetRef,
		       const std::vector< edm::Ptr< L1TTTrackType > >& trkPtrs,
		       float jetvtx = -999. );
      
      virtual ~L1TkJetParticle() {}
      
      // ---------- const member functions ---------------------
      
      const edm::Ref< L1JetParticleCollection >& getJetRef() const { return jetRef_ ; }
      
      const std::vector< edm::Ptr< L1TTTrackType > >& getTrkPtrs() const { return trkPtrs_; }
      
      float getJetVtx() const  { return JetVtx_ ; }

      // ---------- member functions ---------------------------
      
      void setJetVtx(float JetVtx) { JetVtx_ = JetVtx ; }
            
      int bx() const;
      
    private:
      edm::Ref< L1JetParticleCollection > jetRef_ ;
      std::vector< edm::Ptr< L1TTTrackType > >  trkPtrs_;
      
      float JetVtx_ ;
      
    };
}

#endif

