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
<<<<<<< HEAD
	assert( trainingvars_->at(varsptrs_.size())->GetName() == name.c_str() );
=======
	std::cerr << "RooFuncReader::bookVariable " << trainingvars_->getSize() << " " << varsptrs_.size() << std::endl;
	if( std::string(trainingvars_->at(varsptrs_.size())->GetName()) != name ) {
		std::cerr << "RooFuncReader. Error booking variable: expecting " << trainingvars_->at(varsptrs_.size())->GetName() 
			  << " got " << name << std::endl;
		assert( 0 );
	}
>>>>>>> h2gglobe/master
	varsptrs_.push_back(ptr);
}
