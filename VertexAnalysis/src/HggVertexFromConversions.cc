#include "../interface/HggVertexFromConversions.h"
#include "../interface/PhotonInfo.h"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <assert.h>
#include "TVectorD.h"


using namespace std;



HggVertexFromConversions::HggVertexFromConversions(float s1, 
						   float s2, 
						   float s3, 
						   float s4, 
						   float s5, 
						   float s6 ):
  sigmaPix_(s1),
  sigmaTib_(s2),
  sigmaTob_(s3),
  sigmaFwd1_(s4),
  sigmaFwd2_(s5),
  sigmaFwd3_(s6)  
{}


double HggVertexFromConversions::vtxdZ(const PhotonInfo & pho)
{


  // attribute the error depending on the tracker region
  double dz=-99999;

  if ( pho.iDet() ==1 ) { // barrel
    if ( pho.conversionVertex().Perp() <=15 ) {
      dz=sigmaPix_;
    } else if ( pho.conversionVertex().Perp() > 15 && pho.conversionVertex().Perp() <=60 ) {
      dz=sigmaTib_;
    } else {
      dz=sigmaTob_;
    }

  } else { // endcap


    if ( fabs(pho.conversionVertex().Z() ) <=50 ) {
      dz=sigmaFwd1_;
    } else if ( fabs(pho.conversionVertex().Z() ) > 50 && fabs(pho.conversionVertex().Z()) <= 100 ) {
      dz=sigmaFwd2_;
    } else {
      dz=sigmaFwd3_;
    }
  }

  return dz;

}


double HggVertexFromConversions::vtxZ(const PhotonInfo & pho)
{

  // get the z from conversions
  double deltaX1 =  pho.caloPosition().X()- pho.conversionVertex().X();
  double deltaY1 =  pho.caloPosition().Y()- pho.conversionVertex().Y();
  double deltaZ1 =  pho.caloPosition().Z()- pho.conversionVertex().Z();
  double R1 = sqrt(deltaX1*deltaX1+deltaY1*deltaY1);
  double tantheta = R1/deltaZ1;
  
  double deltaX2 = pho.conversionVertex().X()-pho.beamSpot().X();
  double deltaY2 = pho.conversionVertex().Y()-pho.beamSpot().Y();
  double R2 = sqrt(deltaX2*deltaX2+deltaY2*deltaY2);
  double deltaZ2 = R2/tantheta;
  double higgsZ = pho.conversionVertex().Z()-deltaZ1-deltaZ2;
  return higgsZ;

}


