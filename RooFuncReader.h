#ifndef _RooFuncReader_h_
#define _RooFuncReader_h_

#include "RooWorkspace.h"
#include "RooFunctor.h"

class RooFuncReader
{
public:
	RooFuncReader(RooWorkspace * ws, const std::string name, const std::string & trainingvars);

	void bookVariable(const std::string &name, float * ptr);
	
	double eval() {
		for(size_t iv=0; iv<varsbuf_.size(); ++iv) {
			varsbuf_[iv] = *(varsptrs_[iv]);
		}
		return func_->eval(&varsbuf_[0]);
	};
	
private:
	RooWorkspace * ws_;
	RooFunctor * func_;
	
	RooArgList * trainingvars_;

	std::vector<double> varsbuf_;
	std::vector<float*> varsptrs_;
	
};

#endif
