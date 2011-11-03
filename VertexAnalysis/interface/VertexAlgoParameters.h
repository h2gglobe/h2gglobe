#ifndef hgg_VertexAlgoParameters_h
#define hgg_VertexAlgoParameters_h

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

	bool useAllConversions;
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
};

#endif
