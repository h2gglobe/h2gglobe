#include "../interface/EmptyAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
EmptyAnalysis::EmptyAnalysis()  
{
	name_ = "EmptyAnalysis";
}

// ----------------------------------------------------------------------------------------------------
EmptyAnalysis::~EmptyAnalysis() 
{
}

// ----------------------------------------------------------------------------------------------------
void EmptyAnalysis::Term(LoopAll& l) 
{
}

// ----------------------------------------------------------------------------------------------------
void EmptyAnalysis::Init(LoopAll& l) 
{
	PhotonAnalysis::Init(l);
}

// ----------------------------------------------------------------------------------------------------
void EmptyAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{
}

// ----------------------------------------------------------------------------------------------------
bool EmptyAnalysis::SelectEvents(LoopAll&, int)
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void EmptyAnalysis::FillReductionVariables(LoopAll& l, int jentry)
{
}
   
// ----------------------------------------------------------------------------------------------------
bool EmptyAnalysis::SelectEventsReduction(LoopAll&, int)
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
bool EmptyAnalysis::SkimEvents(LoopAll&, int)
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void EmptyAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{
	vtxAna_.branches(outputTree,"vtx_std_");	
}
