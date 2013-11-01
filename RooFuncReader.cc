#include "RooFuncReader.h"

// -----------------------------------------------------------------------------------------------------------
RooFuncReader::RooFuncReader(RooWorkspace * ws, const std::string name, const std::string & trainingvars)
{
	ws_ = ws;
	trainingvars_ = new RooArgList(*ws_->set(trainingvars.c_str()));
	func_ = static_cast<RooAbsReal*>(ws_->function(name.c_str()))->functor(*trainingvars_);	
	varsbuf_.resize(trainingvars_->getSize(),0.);
}

// -----------------------------------------------------------------------------------------------------------
void RooFuncReader::bookVariable(const std::string &name, float * ptr)
{
	assert( trainingvars_->at(varsptrs_.size())->GetName() == name.c_str() );
	varsptrs_.push_back(ptr);
}
