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
  float & sigma1Pix_;
  float & sigma1Tib_;
  float & sigma1Tob_;
  float & sigma1PixFwd_;
  float & sigma1Tid_;
  float & sigma1Tec_;

  float & sigma2Pix_;
  float & sigma2Tib_;
  float & sigma2Tob_;
  float & sigma2PixFwd_;
  float & sigma2Tid_;
  float & sigma2Tec_;

};



#endif

