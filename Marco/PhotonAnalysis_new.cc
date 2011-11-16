#include "../interface/PhotonAnalysis.h"


#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#define PADEBUG 0

using namespace std;



// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::PhotonAnalysis()  : 
	runStatAnalysis(false), doTriggerSelection(false),
	name_("PhotonAnalysis"),
	vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams),
	tmvaPerVtxMethod("BDTG"),
	tmvaPerVtxWeights(""),
	tmvaPerEvtMethod("evtBTG"),
	tmvaPerEvtWeights(""),
	energyCorrectionMethod("DaunceyAndKenzie"), energyCorrected(0), energyCorrectedError(0)
{
	addConversionToMva=true;
	mvaVertexSelection=false;
	useDefaultVertex=false;
	forcedRho = -1.;

	keepPP = true;
	keepPF = true;
	keepFF = true; 
}

// ----------------------------------------------------------------------------------------------------
PhotonAnalysis::~PhotonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void readEnergyScaleOffsets(const std::string &fname, EnergySmearer::energySmearingParameters::eScaleVector &escaleOffsets, 
			    EnergySmearer::energySmearingParameters::phoCatVector &photonCategories, bool data=true
	)
{
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuMap(const char * fname, TDirectory * dir)
{
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::loadPuWeights(int typid, TDirectory * puFile)
{
    
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Init(LoopAll& l) 
{
	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis START"<<endl;

	if(PADEBUG) 
		cout << "InitRealPhotonAnalysis END"<<endl;

	// FIXME book of additional variables
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::Analysis(LoopAll& l, Int_t jentry) {

  //std::cout << "Matteo" << std::endl;
}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::StatAnalysis(LoopAll& l, Int_t jentry) 
{
}




// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
  //vtxAna_.setBranchAdresses(t,"vtx_std_");
  //vtxAna_.getBranches(t,"vtx_std_",s);
}


// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::PreselectPhotons(LoopAll& l, int jentry) 
{
 	
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::FillReductionVariables(LoopAll& l, int jentry) 
{
	if(PADEBUG) 
		cout<<"myFillReduceVar START"<<endl;
  	
	PreselectPhotons(l,jentry);
		
	if(PADEBUG) 
		cout<<"myFillReduceVar END"<<endl;

}

// ----------------------------------------------------------------------------------------------------
bool PhotonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) 
{

	

	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SkimEvents(LoopAll& l, int jentry)
{
	return true;
}

// ----------------------------------------------------------------------------------------------------

bool PhotonAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void PhotonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
}

void PhotonAnalysis::GetRegressionCorrections(LoopAll &l){

}
