#ifndef _SimTracker_SiPhase2Digitizer_SSDigitizerAlgorithm_h
#define _SimTracker_SiPhase2Digitizer_SSDigitizerAlgorithm_h

#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"


// forward declarations
class TrackerTopology;

class SSDigitizerAlgorithm :public Phase2TrackerDigitizerAlgorithm {
 public:
  SSDigitizerAlgorithm(const edm::ParameterSet& conf);
  ~SSDigitizerAlgorithm() override;

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es) override;
  
  //run the algorithm to digitize a single det
  void accumulateSimHits(const std::vector<PSimHit>::const_iterator inputBegin,
                         const std::vector<PSimHit>::const_iterator inputEnd,
                         const size_t inputBeginGlobalIndex,
			 const unsigned int tofBin,
                         const Phase2TrackerGeomDetUnit* pixdet,
                         const GlobalVector& bfield) override;
  bool select_hit(const PSimHit& hit, double tCorr) override;

private:
  enum { SquareWindow, SampledMode, HitDetectMode, SampledOrHitDetectMode, HIPFindingMode};
  double nFactorial(int n);
  double aScalingConstant( int N , int i);
  double cbc3PulsePolarExpansion(double x);
  double signalShape(double x);
  double getSignalScale(double xval);
  void storeSignalShape();
  bool select_hit_sampledMode(const PSimHit& hit, double tCorr);
  bool select_hit_hitDetectMode(const PSimHit& hit, double tCorr);

  int hitDetectionMode_;
  std::vector<double> pulseShapeVec_;
  std::vector<double> pulseShapeParameters_;
  static constexpr float bx_time{25};

};
#endif
