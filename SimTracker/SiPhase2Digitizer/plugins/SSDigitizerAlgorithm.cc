#include <iostream>
#include <cmath>

#include "SimTracker/SiPhase2Digitizer/plugins/SSDigitizerAlgorithm.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

using namespace edm;

void SSDigitizerAlgorithm::init(const edm::EventSetup& es) { es.get<TrackerDigiGeometryRecord>().get(geom_); }
SSDigitizerAlgorithm::SSDigitizerAlgorithm(const edm::ParameterSet& conf)
    : Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
                                      conf.getParameter<ParameterSet>("SSDigitizerAlgorithm")),
      hitDetectionMode_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<int>("HitDetectionMode")),
      pulseShapeParameters_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<std::vector<double> >("PulseShapeParameters")),
      deadTime_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CBCDeadTime")) {
  
  pixelFlag_ = false;
  LogInfo("SSDigitizerAlgorithm ") << "SSDigitizerAlgorithm constructed "
                                   << "Configuration parameters:"
                                   << "Threshold/Gain = "
                                   << "threshold in electron Endcap = " << theThresholdInE_Endcap_
                                   << "threshold in electron Barrel = " << theThresholdInE_Barrel_ << " "
                                   << theElectronPerADC_ << " " << theAdcFullScale_ << " The delta cut-off is set to "
                                   << tMax_ << " pix-inefficiency " << addPixelInefficiency_;
}
SSDigitizerAlgorithm::~SSDigitizerAlgorithm() { LogDebug("SSDigitizerAlgorithm") << "SSDigitizerAlgorithm deleted"; }
void SSDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
                                             std::vector<PSimHit>::const_iterator inputEnd,
                                             const size_t inputBeginGlobalIndex,
                                             const uint32_t tofBin,
                                             const Phase2TrackerGeomDetUnit* pixdet,
                                             const GlobalVector& bfield) {
  // produce SignalPoint's for all SimHit's in detector
  // Loop over hits
  uint32_t detId = pixdet->geographicalId().rawId();
  size_t simHitGlobalIndex = inputBeginGlobalIndex;  // This needs to be stored to create the digi-sim link later

  // find the relevant hits
  std::vector<PSimHit> matchedSimHits;
  std::copy_if(inputBegin, inputEnd, std::back_inserter(matchedSimHits), [detId](auto const& hit) -> bool {
    return hit.detUnitId() == detId;
  });
  // loop over a much reduced set of SimHits
  for (auto const& hit : matchedSimHits) {
    LogDebug("SSDigitizerAlgorithm") << hit.particleType() << " " << hit.pabs() << " " << hit.energyLoss() << " "
                                     << hit.tof() << " " << hit.trackId() << " " << hit.processType() << " "
                                     << hit.detUnitId() << hit.entryPoint() << " " << hit.exitPoint();

    std::vector<DigitizerUtility::EnergyDepositUnit> ionization_points;
    std::vector<DigitizerUtility::SignalPoint> collection_points;

    double signalScale = 1.0;
    // fill collection_points for this SimHit, indpendent of topology
    if (select_hit(hit, (pixdet->surface().toGlobal(hit.localPosition()).mag()/30.), signalScale)) {
      primary_ionization(hit, ionization_points);  // fills ionization_points

      // transforms ionization_points -> collection_points
      drift(hit, pixdet, bfield, ionization_points, collection_points);

      // compute induced signal on readout elements and add to _signal
      // hit needed only for SimHit<-->Digi link
      induce_signal(hit, simHitGlobalIndex, tofBin, pixdet, collection_points);
    }
    ++simHitGlobalIndex;
  }
}
//
// -- Select the Hit for Digitization
//
bool SSDigitizerAlgorithm::select_hit(const PSimHit& hit, double tCorr, double& sigScale) {
  bool result = false; 
  if (hitDetectionMode_ == SSDigitizerAlgorithm::SampledMode) result = select_hit_sampledMode(hit, tCorr, sigScale);
  else if (hitDetectionMode_ == SSDigitizerAlgorithm::LatchedMode) result = select_hit_latchedMode(hit, tCorr, sigScale);
  else {
    double toa = hit.tof() - tCorr;
    if (toa > theTofLowerCut_ && toa < theTofUpperCut_) result = true;
  }
  return result;
}
//
// -- Select Hits in Sampled Mode
//
bool SSDigitizerAlgorithm::select_hit_sampledMode(const PSimHit& hit, double tCorr, double& sigScale)  {
  double toa = hit.tof() - tCorr;
  double sampling_time = bx_time;
  
  DetId det_id = DetId(hit.detUnitId());
  float theThresholdInE = (det_id.subdetId() == StripSubdetector::TOB) 
           ? theThresholdInE_Barrel_ : theThresholdInE_Endcap_;
  
  sigScale = getSignalScale(sampling_time-toa); 
  if (sigScale*hit.energyLoss()/GeVperElectron_ > theThresholdInE) return true;

  return false;
}
//
// -- Select Hits in Hit Detection Mode
//
bool SSDigitizerAlgorithm::select_hit_latchedMode(const PSimHit& hit, double tCorr, double& sigScale) {
  float toa = hit.tof() - tCorr;
  toa -= hit.eventId().bunchCrossing() * bx_time;

  float sampling_time = (-1)*(hit.eventId().bunchCrossing()+1) * bx_time;
  
  DetId det_id = DetId(hit.detUnitId());
  float theThresholdInE = (det_id.subdetId() == StripSubdetector::TOB) 
           ? theThresholdInE_Barrel_
           : theThresholdInE_Endcap_;

  bool lastPulse = true;
  bool aboveThr = false;
  for (float i = deadTime_; i <= bx_time; i++) {
    sigScale = getSignalScale(sampling_time -toa + i); 
    
    if (sigScale*hit.energyLoss()/GeVperElectron_ > theThresholdInE) aboveThr = true;
    else aboveThr = false;
    if (!lastPulse && aboveThr) {
      return true;
    }
    lastPulse = aboveThr;
  }
  return  false;
}

double SSDigitizerAlgorithm::nFactorial(int n) {
  return TMath::Gamma(n+1);
}
double SSDigitizerAlgorithm::aScalingConstant(int N , int i) {
  return std::pow(-1,(double)i) * nFactorial(N) * nFactorial(N+2)/ (nFactorial(N-i) * nFactorial(N+2-i) * nFactorial(i));
}
double SSDigitizerAlgorithm::cbc3PulsePolarExpansion(double x) {
  if (pulseShapeParameters_.size() < 6) return -1;
  double xOffset = pulseShapeParameters_[0];
  double tau     = pulseShapeParameters_[1];
  double r       = pulseShapeParameters_[2];
  double theta   = pulseShapeParameters_[3];
  int nTerms     = static_cast<int>(pulseShapeParameters_[4]);
  
  double fN = 0;
  double xx = x - xOffset;
  if (xx < 0)
    return 0;

  for (int i = 0; i < nTerms; i++) {
    double angularTerm = 0;
    double temporalTerm = 0;
    double rTerm = std::pow(r,i)/(std::pow(tau,2.*i)*nFactorial(i+2));
    for (int j = 0; j <= i; j++) {
      double aij = nFactorial(i)*nFactorial(i+2)/(nFactorial(i-j)*nFactorial(i+2-j)*nFactorial(j));
      angularTerm += std::pow(std::cos(theta),(double)(i-j))*std::pow(std::sin(theta),(double)j);
      temporalTerm += std::pow(-1., (double)j)*aij*std::pow(xx,(double)(i-j))*std::pow(tau,(double)j);
    }
    double fi = rTerm*angularTerm*temporalTerm;
      
    fN += fi;
  }
  return fN;
}
double SSDigitizerAlgorithm::signalShape(double x) {
  double xOffset   = pulseShapeParameters_[0];
  double tau       = pulseShapeParameters_[1];
  //  double r = pulseShapeParameters_[2];
  //  double theta = pulseShapeParameters_[3];
  //  double nTerms = pulseShapeParameters_[4];
  double maxCharge = pulseShapeParameters_[5];
                                                                                                  
  double xx = x - xOffset;
  return maxCharge*(TMath::Exp(-xx/tau)*std::pow(xx/tau,2.)*cbc3PulsePolarExpansion(x));
}
void SSDigitizerAlgorithm::storeSignalShape() {
  for (int i = 0; i < 1000; i++) {
    float val =  i*0.1;
    pulseShapeVec_.push_back(signalShape(val));    
  } 
}
double SSDigitizerAlgorithm::getSignalScale(double xval) {
  double res = 0.0;
  int len = pulseShapeVec_.size();

  if (xval*10 >= len) return res;
  if (xval < 0.0) return res;

  unsigned int lower = std::floor(xval) * 10;
  unsigned int upper = std::ceil(xval) * 10;
  for (size_t i = lower+1; i < upper*10; i++) {
    float val = i*0.1;
    if (val > xval) {
      res = pulseShapeVec_[i-1];
      break;
    }
  } 
  return res;
}
bool SSDigitizerAlgorithm::isAboveThreshold(const DigitizerUtility::SimHitInfo*const hisInfo, float charge, float thr) {
  if (charge > thr) return true;
  else return false;
}
