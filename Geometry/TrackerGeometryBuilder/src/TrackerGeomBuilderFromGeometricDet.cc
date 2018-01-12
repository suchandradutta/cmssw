#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PlaneBuilderForGluedDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripTopologyBuilder.h"
#include "CondFormats/GeometryObjects/interface/PTrackerParameters.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cfloat>
#include <cassert>
using std::vector;
using std::string;

namespace {
  void verifyDUinTG(TrackerGeometry const & tg) {
    int off=0; int end=0;
    for ( int i=1; i!=7; i++) {
      auto det = GeomDetEnumerators::tkDetEnum[i];
      off = tg.offsetDU(det);
      end = tg.endsetDU(det); assert(end>=off); // allow empty subdetectors. Needed for upgrade
      for (int j=off; j!=end; ++j) {
	assert(tg.detUnits()[j]->geographicalId().subdetId()==i);
	assert(GeomDetEnumerators::subDetGeom[tg.detUnits()[j]->subDetector()]==det);
	assert(tg.detUnits()[j]->index()==j);
      }
    }
  }
}

//TrackerGeometry* TrackerGeomBuilderFromGeometricDet::build( const GeometricDet* gd, const PTrackerParameters& ptp, const TrackerTopology* tTopo ) {
TrackerGeometry* TrackerGeomBuilderFromGeometricDet::build( const GeometricDet* gd) {
  //int BIG_PIX_PER_ROC_X = ptp.vpars[2];
  //int BIG_PIX_PER_ROC_Y = ptp.vpars[3];

  int BIG_PIX_PER_ROC_X = 0;
  int BIG_PIX_PER_ROC_Y = 0; // 0 as defined in the vpars vector in trackerParameters.xml, which is all what PTrackerParameters are for!!
    
  thePixelDetTypeMap.clear();
  theStripDetTypeMap.clear();
   
  TrackerGeometry* tracker = new TrackerGeometry(gd);
  std::vector<const GeometricDet*> comp;
  gd->deepComponents(comp);
 
  //if(tTopo)  theTopo = tTopo;

  //define a vector which associate to the detid subdetector index -1 (from 0 to 5) the GeometridDet enumerator to be able to know which type of subdetector it is
  

  /* ALL THIS IS TO SORT THE GEOMETRICDET BY SUBDETECTOR !!
  std::vector<GeometricDet::GDEnumType> gdsubdetmap(6,GeometricDet::unknown); // hardcoded "6" should not be a surprise... 
  GeometricDet::ConstGeometricDetContainer subdetgd = gd->components();
  
  LogDebug("SubDetectorGeometricDetType") << "GeometriDet enumerator values of the subdetectors";
  for(unsigned int i=0;i<subdetgd.size();++i) {
    assert(subdetgd[i]->geographicalId().subdetId()>0 && subdetgd[i]->geographicalId().subdetId()<7);
    gdsubdetmap[subdetgd[i]->geographicalId().subdetId()-1]= subdetgd[i]->type();
    LogTrace("SubDetectorGeometricDetType") << "subdet " << i 
					    << " type " << subdetgd[i]->type()
					    << " detid " <<  subdetgd[i]->geographicalId().rawId()
					    << " subdetid " <<  subdetgd[i]->geographicalId().subdetId();
					    }
*/
  
  std::vector<const GeometricDet*> dets[3];
  std::vector<const GeometricDet*> & armNegSideDets = dets[0]; armNegSideDets.reserve(comp.size());
  std::vector<const GeometricDet*> & DUTDets = dets[1]; DUTDets.reserve(comp.size());
  std::vector<const GeometricDet*> & armPosSideDets  = dets[2];  armPosSideDets.reserve(comp.size());
  



  
  for(auto & i : comp) {
    dets[i->geographicalID().subdetId()-1].emplace_back(i);
  }
  
    /*
  std::vector<const GeometricDet*> armNegSideDets;
  std::vector<const GeometricDet*> DUTDets;
  std::vector<const GeometricDet*> armPosSideDets;
  for (auto & i : comp) {
    if (i->geographicalID().rawId() < 750) armNegSideDets.push_back(i);
    else if (i->geographicalID().rawId() > 750 && i->geographicalID().rawId() < 800) DUTDets.push_back(i);
    else if (i->geographicalID().rawId() > 800) armPosSideDets.push_back(i);
    }*/
  






  //loop on all the six elements of dets and firstly check if they are from pixel-like detector and call buildPixel, then loop again and check if they are strip and call buildSilicon. "unknown" can be filled either way but the vector of GeometricDet must be empty !!
  // this order is VERY IMPORTANT!!!!! For the moment I (AndreaV) understand that some pieces of code rely on pixel-like being before strip-like 
  
  // now building the Pixel-like subdetectors
  /*
  for(unsigned int i=0;i<6;++i) {
    if(gdsubdetmap[i] == GeometricDet::PixelBarrel) 
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::PixelBarrel,
		 false,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
    if(gdsubdetmap[i] == GeometricDet::PixelPhase1Barrel) 
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P1PXB,
		 false,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y);
    // Phase2 case
        if(gdsubdetmap[i] == GeometricDet::PixelPhase2Barrel) 
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P2PXB,
		 true,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y);
	//
    if(gdsubdetmap[i] == GeometricDet::PixelEndCap)
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::PixelEndcap,
		 false,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
    if(gdsubdetmap[i] == GeometricDet::PixelPhase1EndCap)
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P1PXEC,
		 false,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
    if(gdsubdetmap[i] == GeometricDet::PixelPhase2EndCap)
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P2PXEC,
		 true,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
    if(gdsubdetmap[i] == GeometricDet::OTPhase2Barrel) 
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P2OTB,
		 true,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
    if(gdsubdetmap[i] == GeometricDet::OTPhase2EndCap) 
      buildPixel(dets[i],tracker,GeomDetEnumerators::SubDetector::P2OTEC,
		 true,
		 BIG_PIX_PER_ROC_X,
		 BIG_PIX_PER_ROC_Y); 
  }
  //now building Strips
  for(unsigned int i=0;i<6;++i) {
    if(gdsubdetmap[i] == GeometricDet::TIB)   buildSilicon(dets[i],tracker,GeomDetEnumerators::SubDetector::TIB, "barrel");
    if(gdsubdetmap[i] == GeometricDet::TID)   buildSilicon(dets[i],tracker,GeomDetEnumerators::SubDetector::TID, "endcap");
    if(gdsubdetmap[i] == GeometricDet::TOB)   buildSilicon(dets[i],tracker,GeomDetEnumerators::SubDetector::TOB, "barrel");
    if(gdsubdetmap[i] == GeometricDet::TEC)   buildSilicon(dets[i],tracker,GeomDetEnumerators::SubDetector::TEC, "endcap");
  }  
  // and finally the "empty" subdetectors (maybe it is not needed)
  for(unsigned int i=0;i<6;++i) {
    if(gdsubdetmap[i] == GeometricDet::unknown) {
      if(!dets[i].empty()) throw cms::Exception("NotEmptyUnknownSubDet") << "Subdetector " << i+1 << " is unknown but it is not empty: " << dets[i].size();
      buildSilicon(dets[i],tracker,GeomDetEnumerators::tkDetEnum[i+1], "barrel"); // "barrel" is used but it is irrelevant
    }
  }
  */


  // Telescope arm, (-Z) side
  buildPixel(armNegSideDets,tracker,GeomDetEnumerators::SubDetector::P1PXEC,
	     false,
	     BIG_PIX_PER_ROC_X,
	     BIG_PIX_PER_ROC_Y);
  // Telescope DUT containers (contains only 1 DUT here) 
  buildPixel(DUTDets,tracker,GeomDetEnumerators::SubDetector::P2OTB,
	     true,
	     BIG_PIX_PER_ROC_X,
	     BIG_PIX_PER_ROC_Y); 
  // Telescope arm, (+Z) side
  buildPixel(armPosSideDets,tracker,GeomDetEnumerators::SubDetector::P1PXEC,
	     false,
	     BIG_PIX_PER_ROC_X,
	     BIG_PIX_PER_ROC_Y);
  




  buildGeomDet(tracker);//"GeomDet"

  //verifyDUinTG(*tracker);

  return tracker;
}

void TrackerGeomBuilderFromGeometricDet::buildPixel(std::vector<const GeometricDet*>  const & gdv, 
						    TrackerGeometry* tracker,
						    GeomDetType::SubDetector det,
						    bool upgradeGeometry,
						    int BIG_PIX_PER_ROC_X, // in x direction, rows. BIG_PIX_PER_ROC_X = 0 for SLHC
						    int BIG_PIX_PER_ROC_Y) // in y direction, cols. BIG_PIX_PER_ROC_Y = 0 for SLHC
{
  LogDebug("BuildingGeomDetUnits") << " Pixel type. Size of vector: " << gdv.size() 
				   << " GeomDetType subdetector: " << det 
				   << " logical subdetector: " << GeomDetEnumerators::subDetGeom[det]
				   << " big pix per ROC x: " << BIG_PIX_PER_ROC_X << " y: " << BIG_PIX_PER_ROC_Y
				   << " is upgrade: " << upgradeGeometry;
  
  tracker->setOffsetDU(GeomDetEnumerators::subDetGeom[det]);

  for(auto i : gdv){

    std::string const & detName = i->name().fullname();
    if (thePixelDetTypeMap.find(detName) == thePixelDetTypeMap.end()) {
      std::unique_ptr<const Bounds> bounds(i->bounds());
      
      PixelTopology* t = 
	  PixelTopologyBuilder().build(&*bounds,
				       upgradeGeometry,
				       i->pixROCRows(),
				       i->pixROCCols(),
				       BIG_PIX_PER_ROC_X,
				       BIG_PIX_PER_ROC_Y,
				       i->pixROCx(), i->pixROCy());
      
      thePixelDetTypeMap[detName] = new PixelGeomDetType(t,detName,det);
      tracker->addType(thePixelDetTypeMap[detName]);
    }

    PlaneBuilderFromGeometricDet::ResultType plane = buildPlaneWithMaterial(i);
    GeomDetUnit* temp =  new PixelGeomDetUnit(&(*plane),thePixelDetTypeMap[detName],i->geographicalID());

    tracker->addDetUnit(temp);
    tracker->addDetUnitId(i->geographicalID());
  }
  tracker->setEndsetDU(GeomDetEnumerators::subDetGeom[det]);
}

void TrackerGeomBuilderFromGeometricDet::buildSilicon(std::vector<const GeometricDet*>  const & gdv, 
						      TrackerGeometry* tracker,
						      GeomDetType::SubDetector det,
						      const std::string& part)
{ 
  LogDebug("BuildingGeomDetUnits") << " Strip type. Size of vector: " << gdv.size() 
				   << " GeomDetType subdetector: " << det 
				   << " logical subdetector: " << GeomDetEnumerators::subDetGeom[det]
				   << " part " << part;
  
  tracker->setOffsetDU(GeomDetEnumerators::subDetGeom[det]);

  for(auto i : gdv){

    std::string const & detName = i->name().fullname();
    if (theStripDetTypeMap.find(detName) == theStripDetTypeMap.end()) {
       std::unique_ptr<const Bounds> bounds(i->bounds());
       StripTopology* t =
	   StripTopologyBuilder().build(&*bounds,
				       i->siliconAPVNum(),
				       part);
      theStripDetTypeMap[detName] = new  StripGeomDetType( t,detName,det,
						   i->stereo());
      tracker->addType(theStripDetTypeMap[detName]);
    }
     
    double scale  = (partnerDetId(i->geographicalID())) ? 0.5 : 1.0 ;	// theTopo

    PlaneBuilderFromGeometricDet::ResultType plane = buildPlaneWithMaterial(i,scale);  
    GeomDetUnit* temp = new StripGeomDetUnit(&(*plane), theStripDetTypeMap[detName],i->geographicalID());
    
    tracker->addDetUnit(temp);
    tracker->addDetUnitId(i->geographicalID());
  }  
  tracker->setEndsetDU(GeomDetEnumerators::subDetGeom[det]);

}


void TrackerGeomBuilderFromGeometricDet::buildGeomDet(TrackerGeometry* tracker){

  PlaneBuilderForGluedDet gluedplaneBuilder;
  auto const & gdu = tracker->detUnits();
  auto const & gduId = tracker->detUnitIds();

  for(u_int32_t i=0;i<gdu.size();i++){

    tracker->addDet(gdu[i]);
    tracker->addDetId(gduId[i]);
    string gduTypeName = gdu[i]->type().name();

    //this step is time consuming >> TO FIX with a MAP?
    if( (gduTypeName.find("Ster")!=std::string::npos || 
         gduTypeName.find("Lower")!=std::string::npos) && 
        (glued(gduId[i])!=0 || stack(gduId[i])!=0 )) {  // theTopo
    

      std::cout << "TrackerGeomBuilderFromGeometricDet::buildGeomDet  gduId[i].rawId() = " << gduId[i].rawId() << std::endl;
      std::cout << "TrackerGeomBuilderFromGeometricDet::buildGeomDet  partnerDetId(gduId[i]).rawId() = " << partnerDetId(gduId[i]).rawId() << std::endl;


      int partner_pos=-1;
      for(u_int32_t jj=0;jj<gduId.size();jj++){
  	  if(partnerDetId(gduId[i]) == gduId[jj]) {   // theTopo
  	    partner_pos=jj;
  	    break;
  	  }
      }
      if(partner_pos==-1){
	  throw cms::Exception("Configuration") <<"Module Type is Stereo or Lower but no partner detector found \n"
					        <<"There is a problem on Tracker geometry configuration\n";
      }

      const GeomDetUnit* dus = gdu[i];
      const GeomDetUnit* dum = gdu[partner_pos];
      std::vector<const GeomDetUnit *> composed(2);
      composed[0]=dum;
      composed[1]=dus;
      DetId composedDetId;
      if(gduTypeName.find("Ster")!=std::string::npos){

        PlaneBuilderForGluedDet::ResultType plane = gluedplaneBuilder.plane(composed);
        composedDetId = glued(gduId[i]);  // theTopo
        GluedGeomDet* gluedDet = new GluedGeomDet(&(*plane),dum,dus,composedDetId);
        tracker->addDet((GeomDet*) gluedDet);
        tracker->addDetId(composedDetId);

      } else if (gduTypeName.find("Lower")!=std::string::npos){

        //FIXME::ERICA: the plane builder is built in the middle...
        PlaneBuilderForGluedDet::ResultType plane = gluedplaneBuilder.plane(composed);
        composedDetId = stack(gduId[i]);   // theTopo
        StackGeomDet* stackDet = new StackGeomDet(&(*plane),dum,dus,composedDetId);
	std::cout << "TrackerGeomBuilderFromGeometricDet::buildGeomDet  composedDetId.rawId() = " << composedDetId.rawId() << std::endl;
        tracker->addDet((GeomDet*) stackDet);
        tracker->addDetId(composedDetId);

      } 

    }

  }
}



uint32_t TrackerGeomBuilderFromGeometricDet::glued(const DetId &id) const {

    uint32_t subdet=id.subdetId();
    /*
    if ( subdet == PixelSubdetector::PixelBarrel )
      return 0;
    if ( subdet == PixelSubdetector::PixelEndcap )
      return 0;
    if ( subdet == StripSubdetector::TIB )
      return tibGlued(id);
    if ( subdet == StripSubdetector::TID )
      return tidGlued(id);
    if ( subdet == StripSubdetector::TOB )
      return tobGlued(id);
    if ( subdet == StripSubdetector::TEC )
      return tecGlued(id);

    throw cms::Exception("Invalid DetId") << "Unsupported DetId in TrackerTopology::glued";
    */

    // Telescope arm, (-Z) side
    if (subdet == 1) { return 0;  }    

    // Telescope DUT containers (contains only 1 DUT here)    
    else if (subdet == 2) {   
      if ( (id.rawId() % 2) == 1) { return DetId( id.rawId() - 1 ); }  // inner sensor = 1
      else { return DetId( id.rawId() - 2 ); } // outer sensor = 2
    }    

    // Telescope arm, (+Z) side
    else if (subdet == 3) { return 0; }  
  
    else { std::cout << "Unexpected subdet = " << subdet << std::endl; return 0; }



    return 0;
}


uint32_t TrackerGeomBuilderFromGeometricDet::stack(const DetId &id) const {

    uint32_t subdet=id.subdetId();
    /*
    if ( subdet == PixelSubdetector::PixelBarrel )
      return 0;
      if ( subdet == PixelSubdetector::PixelEndcap )
      return 0;
      if ( subdet == StripSubdetector::TIB )
      return tibStack(id);
      if ( subdet == StripSubdetector::TID )
      return tidStack(id);
      if ( subdet == StripSubdetector::TOB )
      return tobStack(id);
      if ( subdet == StripSubdetector::TEC )
      return tecStack(id);

      throw cms::Exception("Invalid DetId") << "Unsupported DetId in TrackerTopology::stack";
    */


    // Telescope arm, (-Z) side
    if (subdet == 1) { return 0;  }    

    // Telescope DUT containers (contains only 1 DUT here)    
    else if (subdet == 2) {   
      if ( (id.rawId() % 2) == 1) { return DetId( id.rawId() - 1 ); }  // inner sensor = 1
      else { return DetId( id.rawId() - 2 ); } // outer sensor = 2
    }    

    // Telescope arm, (+Z) side
    else if (subdet == 3) { return 0; }  
  
    else { std::cout << "Unexpected subdet = " << subdet << std::endl; return 0; }


}


DetId TrackerGeomBuilderFromGeometricDet::partnerDetId(const DetId &id) const {

  uint32_t subdet=id.subdetId();
  /*
    if ( subdet == PixelSubdetector::PixelBarrel )
    return 0;
    if ( subdet == PixelSubdetector::PixelEndcap )
    return 0;
    if ( subdet == StripSubdetector::TIB )
    return tibPartnerDetId(id);
    if ( subdet == StripSubdetector::TID )
    return tidPartnerDetId(id);
    if ( subdet == StripSubdetector::TOB )
    return tobPartnerDetId(id);
    if ( subdet == StripSubdetector::TEC )
    return tecPartnerDetId(id);

    throw cms::Exception("Invalid DetId") << "Unsupported DetId in TrackerTopology::partnerDetId";
  */

  // Telescope arm, (-Z) side
  if (subdet == 1) { return 0;  }    

  // Telescope DUT containers (contains only 1 DUT here)    
  else if (subdet == 2) {   
    if ( (id.rawId() % 2) == 1) { return DetId( id.rawId() + 1 ); }  // inner sensor, hence return outer
    else { return DetId( id.rawId() - 1 ); } // outer sensor, hence return inner
  }    

  // Telescope arm, (+Z) side
  else if (subdet == 3) { return 0; }  
  
  else { std::cout << "Unexpected subdet = " << subdet << std::endl; return 0; }

  return 0;
}




PlaneBuilderFromGeometricDet::ResultType
TrackerGeomBuilderFromGeometricDet::buildPlaneWithMaterial(const GeometricDet* gd,
							   double scale) const
{
  PlaneBuilderFromGeometricDet planeBuilder;
  PlaneBuilderFromGeometricDet::ResultType plane = planeBuilder.plane(gd);  
  //
  // set medium properties (if defined)
  //
  plane->setMediumProperties(MediumProperties(gd->radLength()*scale,gd->xi()*scale));

  return plane;
}
