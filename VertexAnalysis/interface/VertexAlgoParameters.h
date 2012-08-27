#ifndef hgg_VertexAlgoParameters_h
#define hgg_VertexAlgoParameters_h

/** parameters for the vertex selection algorithms */
class VertexAlgoParameters {
public:
	bool fixTkIndex;

	bool rescaleTkPtByError;
	float trackCountThr;
	bool highPurityOnly;
	float maxD0Signif, maxDzSignif;
	/// float maxD0, maxDz;
	/// float minPt;
        bool removeTracksInCone;
        float coneSize;

	std::string vtxProbFormula; 
	
	int useAllConversions;
	float sigma1Pix;
	float sigma1Tib;
	float sigma1Tob;
	float sigma1PixFwd;
	float sigma1Tid;
	float sigma1Tec;

  	float sigma2Pix;
	float sigma2Tib;
	float sigma2Tob;
	float sigma2PixFwd;
	float sigma2Tid;
	float sigma2Tec;

  	float singlelegsigma1Pix;
	float singlelegsigma1Tib;
	float singlelegsigma1Tob;
	float singlelegsigma1PixFwd;
	float singlelegsigma1Tid;
	float singlelegsigma1Tec;
  
	float singlelegsigma2Pix;
	float singlelegsigma2Tib;
	float singlelegsigma2Tob;
	float singlelegsigma2PixFwd;
	float singlelegsigma2Tid;
	float singlelegsigma2Tec;
};

#endif
