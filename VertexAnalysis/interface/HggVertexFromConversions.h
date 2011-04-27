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
class VertexAlgoParameters;

class HggVertexFromConversions
{
public:

  HggVertexFromConversions(VertexAlgoParameters & ap);
  
  double vtxZ(const PhotonInfo & pho);
  double vtxdZ(const PhotonInfo & pho);

 private:
  float & sigmaPix_;
  float & sigmaTib_;
  float & sigmaTob_;
  float & sigmaFwd1_;
  float & sigmaFwd2_;
  float & sigmaFwd3_;

};



#endif

