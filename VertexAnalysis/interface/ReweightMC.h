#ifndef hgg_ReweightMC_h
#define hgg_ReweightMC_h

#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include <vector>
#include <string>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
/**
 *
 * \class ReweightMC
 *
 * usage example:
 * <code>
 * </code>
 */
class ReweightMC
{
public:

	enum constants {nVtxDefault=30};
	typedef enum target {Gen=0, Reco=1} target_t;
	static const TString UID;

	ReweightMC(target_t target = Gen, int nVtx=nVtxDefault);
	ReweightMC(TH1* hMCGen, float GenPoissonMean=18.0);
	ReweightMC(TH1* hMCGen, TH1* hGenTarget);
	ReweightMC(TH2* hMCGenReco, TH1* hMCGen, float RecoPoissonMean=18.0);
	ReweightMC(TH2* hMCGenReco, TH1* hMCGen, TH1* hRecoTarget);
	~ReweightMC();

	int getnVtx(){ return nVtx_;};

	int fillMC(int nGen, int nReco, double weight=1.0);
	int fillMC(int nGen, double weight=1.0);
	void makeResponse();
	TH1 *getResponse(){ return hMCGenOrGenReco_;};
	TH1 *getMCGen(){ return hMCGen_;};

	int	fillTarget(int n, double weight=1.0);
	TH1 *getTarget(){return hTarget_;};

	double getWeight(int nGenInMC, int nRecoInMC=0);
	TH1 *getWeights(){ return hWeights_;};

private:
	void makeWeights();

	target_t target_;
	int nVtx_;
	TH1 *hMCGenOrGenReco_;
	union histo_t{
		TH2 *MCGenReco;
		TH1 *MCGen;
	} h_;
	TH1 *hMCGen_;
	TH1 *hTarget_;
	TH1 *hWeights_;
	vector<TH1* > toDelete_;
};


#endif

// Local Variables:
// mode: c++
// c-basic-offset: 8
// End:
