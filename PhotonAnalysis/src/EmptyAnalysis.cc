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
bool EmptyAnalysis::Analysis(LoopAll& l, Int_t jentry) 
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


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
