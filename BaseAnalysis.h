#ifndef __BASEANALYSIS__
#define __BASEANALYSIS__

#include "LoopAll.h"

// ------------------------------------------------------------------------------------
/**
 * \class BaseAnalysis
 *
 * Defines the interface to be implemented by analysed classes loaded 
 * by LoopAll.
 * 
 * The analysis model is made by the two steps:
 * - reduction step: events are skimmed and additional variables are added 
 *   to the analysis tree.
 * - analysis step: the reduced events are processed and further selected.
 *   Final plots and tuples are produced in this step.
 */
class BaseAnalysis 
{
public:
	// ! C-TOR
	BaseAnalysis();
	
	// ! D-TOR
	virtual ~BaseAnalysis();
	
	// ! the class name
	virtual const std::string & name() const = 0; 
	
	// !  
	operator const std::string & () const { return this->name(); };
	
	// ! Initialization method. Called before the events loop.
	virtual void Init(LoopAll&) = 0;
	// ! Finilization method. Called at the end of the event loop.
	virtual void Term(LoopAll&) = 0;
	// ! Method called at the start of each new file, call to reset any random variables etc
	virtual void ResetAnalysis() =0;
	
	// ! Method to use in order to book additional veriables in the output tree 
	virtual void ReducedOutputTree(LoopAll &, TTree *) = 0;
	// ! Method used to manually book additional input branches to be read
	virtual void GetBranches(TTree *, std::set<TBranch *>& ) = 0;
	
	// ! Quick event skimming. No inputs are read at this stage yet. So the user has to manually call the GetEntry methods.
	virtual bool SkimEvents(LoopAll&, int) = 0;
	// ! Fill additional variables used for the reduction step.
	virtual void FillReductionVariables(LoopAll&, int) = 0;
	// ! Select events for the reduction step.
	virtual bool SelectEventsReduction(LoopAll&, int) = 0;
	// ! Select events in the analysis step.
	virtual bool SelectEvents(LoopAll&, int) = 0;
	
	// ! Implemement your final analysis here
	virtual bool Analysis(LoopAll&, Int_t) = 0;
	
};

// ! Used to search analyzers by name 
bool operator == (BaseAnalysis * a, const std::string & b);

#endif
