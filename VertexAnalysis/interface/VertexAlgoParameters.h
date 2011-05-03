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

	float sigmaPix;
	float sigmaTib;
	float sigmaTob;
	float sigmaFwd1;
	float sigmaFwd2;
	float sigmaFwd3;
};

#endif
