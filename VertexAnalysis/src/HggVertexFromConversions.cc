#include "../interface/HggVertexFromConversions.h"
#include "../interface/PhotonInfo.h"
#include "../interface/VertexAlgoParameters.h"

#include "stdio.h"
#include "math.h"

#include <algorithm>
#include <numeric>
#include <iterator>

#include <assert.h>
#include "TVectorD.h"


using namespace std;



HggVertexFromConversions::HggVertexFromConversions(VertexAlgoParameters &ap):
  sigma1Pix_ (ap.sigma1Pix ),
  sigma1Tib_ (ap.sigma1Tib ),
  sigma1Tob_ (ap.sigma1Tob ),
  sigma1PixFwd_(ap.sigma1PixFwd),
  sigma1Tid_(ap.sigma1Tid),
  sigma1Tec_(ap.sigma1Tec),
  
  sigma2Pix_ (ap.sigma2Pix ),
  sigma2Tib_ (ap.sigma2Tib ),
  sigma2Tob_ (ap.sigma2Tob ),
  sigma2PixFwd_(ap.sigma2PixFwd),
  sigma2Tid_(ap.sigma2Tid),
  sigma2Tec_(ap.sigma2Tec)
{
}


double HggVertexFromConversions::vtxdZ(const PhotonInfo & pho)
{


  // attribute the error depending on the tracker region
  double dz=-99999;

  if ( pho.iDet() ==1 ) { // barrel
    if ( pho.conversionVertex().Perp() <=15 ) {
      dz=sigma2Pix_;
    } else if ( pho.conversionVertex().Perp() > 15 && pho.conversionVertex().Perp() <=60 ) {
      dz=sigma2Tib_;
    } else {
      dz=sigma2Tob_;
    }

  } else { // endcap


    if ( fabs(pho.conversionVertex().Z() ) <=50 ) {
      dz=sigma2PixFwd_;
    } else if ( fabs(pho.conversionVertex().Z() ) > 50 && fabs(pho.conversionVertex().Z()) <= 100 ) {
      dz=sigma2Tid_;
    } else {
      dz=sigma2Tec_;
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
  double higgsZ =  pho.caloPosition().Z()-deltaZ1-deltaZ2;
  return higgsZ;

}


