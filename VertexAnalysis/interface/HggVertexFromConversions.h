#ifndef hgg_VertexFromConversions_h
#define hgg_VertexFromConversions_h

#include <vector>
#include <map>
#include <string>

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMatrixDSym.h"

class PhotonInfo;

class HggVertexFromConversions
{
public:

  HggVertexFromConversions(float s1, 
			   float s2, 
			   float s3, 
			   float s4, 
			   float s5, 
			   float s6 );
  
 private:
  double vtxZ(const PhotonInfo & pho);
  double vtxdZ(const PhotonInfo & pho);
  float sigmaPix_;
  float sigmaTib_;
  float sigmaTob_;
  float sigmaFwd1_;
  float sigmaFwd2_;
  float sigmaFwd3_;

};



#endif

