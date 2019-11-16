#include <iostream>
#include <cmath>

#include "SimTracker/SiPhase2Digitizer/plugins/SSDigitizerAlgorithm.h"
#include "SimTracker/Common/interface/SiG4UniversalFluctuation.h"
#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineSimService.h"


// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"


using namespace edm;

void SSDigitizerAlgorithm::init(const edm::EventSetup& es) {
  es.get<TrackerDigiGeometryRecord>().get(geom_);
}
SSDigitizerAlgorithm::SSDigitizerAlgorithm(const edm::ParameterSet& conf) :
  Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
				  conf.getParameter<ParameterSet>("SSDigitizerAlgorithm")),
  hitDetectionMode_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<int>("HitDetectionMode")),
  pulseShapeParameters_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<std::vector<double> >("PulseShapeParameters")),
  deadTime_(conf.getParameter<ParameterSet>("SSDigitizerAlgorithm").getParameter<double>("CBCDeadTime"))
{
  pixelFlag = false;
  LogInfo("SSDigitizerAlgorithm ") << "SSDigitizerAlgorithm constructed "
				   << "Configuration parameters:"
				   << "Threshold/Gain = "
				   << "threshold in electron Endcap = "
				   << theThresholdInE_Endcap
				   << "threshold in electron Barrel = "
				   << theThresholdInE_Barrel
				   <<" " << theElectronPerADC << " " << theAdcFullScale
				   << " The delta cut-off is set to " << tMax
				   << " pix-inefficiency "<<AddPixelInefficiency;


  storeSignalShape();
}
SSDigitizerAlgorithm::~SSDigitizerAlgorithm() {
  LogDebug("SSDigitizerAlgorithm") << "SSDigitizerAlgorithm deleted";
}
void SSDigitizerAlgorithm::accumulateSimHits(std::vector<PSimHit>::const_iterator inputBegin,
					     std::vector<PSimHit>::const_iterator inputEnd,
					     const size_t inputBeginGlobalIndex,
                                             const unsigned int tofBin,
					     const Phase2TrackerGeomDetUnit* pixdet,
					     const GlobalVector& bfield) {
  // produce SignalPoint's for all SimHit's in detector
  // Loop over hits
  uint32_t detId = pixdet->geographicalId().rawId();
  size_t simHitGlobalIndex = inputBeginGlobalIndex; // This needs to be stored to create the digi-sim link later
  for (auto it = inputBegin; it != inputEnd; ++it, ++simHitGlobalIndex) {
    // skip hits not in this detector.
    if ((*it).detUnitId() != detId)
      continue;

    LogDebug("SSDigitizerAlgorithm")
      << (*it).particleType() << " " << (*it).pabs() << " "
      << (*it).energyLoss() << " " << (*it).tof() << " "
      << (*it).trackId() << " " << (*it).processType() << " "
      << (*it).detUnitId()
      << (*it).entryPoint() << " " << (*it).exitPoint() ;
      
    std::vector<DigitizerUtility::EnergyDepositUnit> ionization_points;
    std::vector<DigitizerUtility::SignalPoint> collection_points;
    // fill collection_points for this SimHit, indpendent of topology
    if (select_hit(*it, (pixdet->surface().toGlobal((*it).localPosition()).mag()/30.))) { 
      primary_ionization(*it, ionization_points); // fills _ionization_points
      drift(*it, pixdet, bfield, ionization_points, collection_points);  // transforms _ionization_points to collection_points

      // compute induced signal on readout elements and add to _signal

      induce_signal(*it, simHitGlobalIndex, tofBin, pixdet, collection_points); // *ihit needed only for SimHit<-->Digi link
    }
  }
}
//
// -- Select the Hit for Digitization
//
bool SSDigitizerAlgorithm::select_hit(const PSimHit& hit, double tCorr) {
  bool result = false; 
  if (hitDetectionMode_ == SSDigitizerAlgorithm::SampledMode) result = select_hit_sampledMode(hit, tCorr);
  else if (hitDetectionMode_ == SSDigitizerAlgorithm::LachedMode) result = select_hit_lachedMode(hit, tCorr);
  else {
    if ((hit.tof() - tCorr) >= theTofLowerCut && (hit.tof() - tCorr) <= theTofUpperCut) result = true;
  }
  return result;
}
//
// -- Select Hits in Sampled Mode
//
bool SSDigitizerAlgorithm::select_hit_sampledMode(const PSimHit& hit, double tCorr)  {
  double toa = hit.tof() - tCorr;

  double sampling_time = bx_time;

  float theThresholdInE = 0.;
  DetId det_id = DetId(hit.detUnitId());
  if (det_id.subdetId() == StripSubdetector::TOB) theThresholdInE = theThresholdInE_Barrel;
  else theThresholdInE = theThresholdInE_Endcap;
  
  double scale = getSignalScale(sampling_time-toa); 
  if (scale*hit.energyLoss()/GeVperElectron  > theThresholdInE) return true;
  else return false;;
}
//
// -- Select Hits in Hit Detection Mode
//
bool SSDigitizerAlgorithm::select_hit_lachedMode(const PSimHit& hit, double tCorr) {
  float toa = hit.tof() - tCorr;
  toa -= hit.eventId().bunchCrossing() * bx_time;

  float sampling_time = (-1)*(hit.eventId().bunchCrossing()+1) * bx_time ;

  float theThresholdInE = 0.;
  DetId det_id = DetId(hit.detUnitId());
  if (det_id.subdetId() == StripSubdetector::TOB) theThresholdInE = theThresholdInE_Barrel;
  else theThresholdInE = theThresholdInE_Endcap;
  bool lastPulse = true;
  bool aboveThr = false;
  for (float i = deadTime_; i <= bx_time; i++) {
    double scale = getSignalScale(sampling_time -toa + i ); 

    if (scale*hit.energyLoss()/GeVperElectron > theThresholdInE) aboveThr = true;
    else aboveThr = false;
    if (!lastPulse && aboveThr) return true;
    lastPulse = aboveThr;
  }
  return  false;
}

double SSDigitizerAlgorithm::nFactorial(int n) {
  return TMath::Gamma(n+1);
}
double SSDigitizerAlgorithm::aScalingConstant( int N , int i) {

  return std::pow(-1,(double)i) * nFactorial(N) * nFactorial(N+2)/ (nFactorial(N-i) * nFactorial(N+2-i) * nFactorial(i));
}
double SSDigitizerAlgorithm::cbc3PulsePolarExpansion(double x)
{
  if (pulseShapeParameters_.size() < 6) return -1;
  double xOffset = pulseShapeParameters_[0];
  double tau = pulseShapeParameters_[1];
  double r = pulseShapeParameters_[2];
  double theta = pulseShapeParameters_[3];
  int nTerms = (int)pulseShapeParameters_[4];
  
  double fN=0;
  double xx = x - xOffset;
  if( xx  < 0 )
    return 0;

  for( int i = 0 ; i < nTerms ; i++)
    {
      double angularTerm=0;
      double temporalTerm=0;
      double rTerm = std::pow(r,i)/(std::pow(tau,2.*i)*nFactorial(i+2));
      for( int j = 0 ; j <= i ; j++ )
	{
	  double aij = nFactorial(i)*nFactorial(i+2)/(nFactorial(i-j)*nFactorial(i+2-j)*nFactorial(j));
	  angularTerm += std::pow(std::cos(theta),(double)(i-j))*std::pow(std::sin(theta),(double)j);
	  temporalTerm += std::pow(-1., (double)j)*aij*std::pow(xx,(double)(i-j))*std::pow(tau,(double)j);
	}
      double fi = rTerm*angularTerm*temporalTerm;
      
      fN += fi;
    }
  return fN;
}
double SSDigitizerAlgorithm::signalShape(double x)
{
  double xOffset = pulseShapeParameters_[0];
  double tau = pulseShapeParameters_[1];
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
  double res =0.0;
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
