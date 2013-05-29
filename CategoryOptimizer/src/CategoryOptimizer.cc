#include "../interface/CategoryOptimizer.h"
#include "TMath.h"

// ---------------------------------------------------------------------------------------------
GenericFigureOfMerit::GenericFigureOfMerit(std::vector<AbsModelBuilder *> & sig, std::vector<AbsModelBuilder *> & bkg, 
					   AbsFomProvider * fom, 
					   int ndim, int nbound, const double * cutoffs,  int northocuts, 
					   bool addConstraint, 
					   bool telescopicBoundaries, 
					   const std::vector<HistoConverter *> & transformations) : 
	sigModels_(sig), bkgModels_(bkg), fom_(fom), ndim_(ndim), nbound_(nbound), northocuts_(northocuts),
	cutoffs_(cutoffs,cutoffs+ndim),
	addConstraint_(addConstraint), telescopicBoundaries_(telescopicBoundaries),
	transformations_(transformations)
{
	allModels_.reserve(sigModels_.size()+bkgModels_.size());
	std::copy(sigModels_.begin(),sigModels_.end(),std::back_inserter(allModels_));
	std::copy(bkgModels_.begin(),bkgModels_.end(),std::back_inserter(allModels_));
}


// ---------------------------------------------------------------------------------------------
double GenericFigureOfMerit::operator() (double *x, double *p) const
{
	// sort out the input parameters
	std::vector<double *> deltas(ndim_);
	std::vector<double> firstb(ndim_);
	std::vector<double> extfirstb(ndim_);
	std::vector<double> lastb(ndim_);
	std::vector<double> mindim(ndim_);
	const double * cutoffs = (p != 0 ? p : &cutoffs_[0]);

	for(int idim=0; idim<ndim_; ++idim) {
		deltas[idim] = &x[idim*nbound_+1];
		firstb[idim] = x[idim*nbound_];
		if( ! std::isfinite(firstb[idim]) ) { return 1000.; }
		extfirstb[idim] = x[idim*nbound_];
		lastb[idim] = x[idim*nbound_];
		mindim[idim] = ( !transformations_.empty() && transformations_[idim]!=0 ? 0. : sigModels_[0]->getMin(idim) );
	}
	
	std::vector<double> orthocuts(x+ndim_*(nbound_+addConstraint_),x+ndim_*(nbound_+addConstraint_)+northocuts_);
	
	// initialize the model builders
	std::ostream_iterator<double> output(cout, ",");
	CategoryOptimizer::doTransform(transformations_,&extfirstb[0]);
	for(std::vector<AbsModelBuilder *>::const_iterator imod=allModels_.begin(); imod!=allModels_.end(); ++imod ) {
		if( northocuts_ > 0 ) {
			(*imod)->setOrthoCuts(&orthocuts[0]);
		}
		(*imod)->beginIntegration(&firstb[0]);
	}
	
	double ret = 0.;
	// loop over the categories
	std::vector<double> currb = lastb;
	std::vector<double> distances;
	for(int ii=0; ii<nbound_; ++ii) {
		int jj=ii+1;
		std::vector<double> newb = currb;
		
		// work out the new boundary starting from the last one
		bool found = false;
		for(; jj<nbound_; ++jj) {
			for(int idim=0; idim<ndim_; ++idim) {
				if( telescopicBoundaries_ ) {
					newb[idim] -= deltas[idim][jj-1];
				} else {
					newb[idim]  = deltas[idim][jj-1];
				}
				if( ! std::isfinite(newb[idim]) ) { return 1000.; }
				// ignore boundaries outside of the ranges
				if( newb[idim] < mindim[idim] ) { newb[idim] = mindim[idim]; }
				// make sure that at least in one dimension the boundaries are far enough
				if( lastb[idim] > newb[idim] + cutoffs[idim] ) { found = true; }
			}
			if( found ) { break; }
			else { 
				double distance = 0.;
				for(int idim=0; idim<ndim_; ++idim) {
					/// double idist = std::max(lastb[idim]-newb[idim],0.)/cutoffs[idim]-1.;
					double idist = std::max(lastb[idim]-newb[idim],0.)/cutoffs[idim];
					distance += idist*idist;
			 	}
				distance = sqrt(distance);
				ret += 6.52896*(1.-TMath::Erf(100.*(distance-1.+0.01)))/(distance*distance*distance+0.3*0.3*0.3);
				//// ret += 100./distance; /// FIXME should be configurable and continuos at the boundary
			}
		}
		
		ii = jj - 1;
		if( found ) {
			lastb = newb;
		}
		currb = newb;
		
		// define a new set of boundaries
		CategoryOptimizer::doTransform(transformations_,&newb[0]);
		if( found ) {
			for(std::vector<AbsModelBuilder *>::const_iterator imod=allModels_.begin(); imod!=allModels_.end(); ++imod ) {
				(*imod)->addBoundary(&newb[0]);
			}
		} 
	}
	
	for(std::vector<AbsModelBuilder *>::const_iterator imod=allModels_.begin(); imod!=allModels_.end(); ++imod ) {
		(*imod)->endIntegration();
	}
	
	// retrieve the signal and background modes
	std::vector<AbsModel *> sigModels, bkgModels;
	for(std::vector<AbsModelBuilder *>::const_iterator isig=sigModels_.begin(); isig!=sigModels_.end(); ++isig ) {
		sigModels.push_back( (*isig)->getModel() );
	}
	for(std::vector<AbsModelBuilder *>::const_iterator ibkg=bkgModels_.begin(); ibkg!=bkgModels_.end(); ++ibkg ) {
		bkgModels.push_back( (*ibkg)->getModel() );
	}

	// compute the FOM
	ret += (*fom_)(sigModels,bkgModels);

	// additional constraint
	if( addConstraint_ ) {
		for(int idim=0; idim<ndim_; ++idim ) {
			ret +=  x[ndim_*nbound_+idim]*pow( (lastb[idim]-mindim[idim]), 2. );
		}
	}

	return ret;
}

// ---------------------------------------------------------------------------------------------
double CategoryOptimizer::optimizeNCat(int ncat, const double * cutoffs, bool dryrun, bool debug)
{
	int nbound = ncat+1;
	
	std::vector<double> tmpcutoffs(cutoffs,cutoffs+ndim_);
	bool build = transformations_.empty() && ! transformModels_.empty();
	if( build ) { 
		std::cout << "Buildinig variable transformations" << std::endl; 
	}
	for(int idim=0; idim<ndim_; ++idim) {
		if( build ) {
			TH1 * transformPdf = transformModels_[0]->getPdf(idim);
			for(size_t imod=1; imod<transformModels_.size(); ++imod) {
				TH1 * itransformPdf = transformModels_[imod]->getPdf(idim);
				transformPdf->Add(itransformPdf);
				delete itransformPdf;
			}
			HistoConverter * conv = cdfInv(transformPdf,
						       transformModels_[0]->getMin(idim),
						       transformModels_[0]->getMax(idim));
			setTransformation(idim,conv);
			delete transformPdf;
		}
		/// FIXME need to transform the cutoffs
		/// if( ! transformations_.empty() && transformations_[idim] != 0 ) {
		/// 	HistoConverter * conv = transformations_[idim];
		/// 	double median = conv->eval(0.5);
		/// 	tmpcutoffs[idim] = fabs( conv->eval(0.5+0.5*cutoffs[idim]) - 
		/// 				 conv->eval(0.5-0.5*cutoffs[idim]) );
		/// }
	}

	
	// Book the FOM
	GenericFigureOfMerit theFom(sigModels_, bkgModels_, fom_, ndim_, nbound, 
				    &tmpcutoffs[0], orthocuts_.size(),
				    addConstraint_, telescopicBoundaries_, transformations_);
	minimizer_->SetFunction(theFom);
	
	// Define category boundaries. 
	// Last boundary fixed to the maximum range in each dimension
	for(int idim=0; idim<ndim_; ++idim) {
		double min = sigModels_[0]->getMin(idim);
		double max = sigModels_[0]->getMax(idim);
		if( ! transformations_.empty() && transformations_[idim]!=0 ) {
			min = 0.; max = 1.;
		}
		double range = max - min;
		if( ! floatFirst_ ) { 
			minimizer_->SetFixedVariable(idim*nbound, Form("firstBoundDim%d",idim), 
						     max);
		} else {
			minimizer_->SetLimitedVariable(idim*nbound, Form("firstBoundDim%d",idim), 
						       max, tmpcutoffs[idim]*0.5, 
						       min+ncat*tmpcutoffs[idim], max );
		}
		for(int icat=0; icat<ncat; ++icat) {
			if( telescopicBoundaries_ ) {
				minimizer_->SetLimitedVariable(idim*nbound + icat +1,
							       Form( "deltaBoundDim%dBin%d", idim, icat ), 
							       range/(double)ncat,   // FIXME smarter initialization
							       tmpcutoffs[idim]*0.5, tmpcutoffs[idim], range );
			} else {
				minimizer_->SetLimitedVariable(idim*nbound + icat +1,
							       Form( "absBoundDim%dBin%d", idim, icat ), 
							       max - (double)(icat+1)*range/(double)ncat,
							       tmpcutoffs[idim]*0.5, 
							       min, max );
			}
		}
	}
	
	// penalty to constrain the lower boundaires
	if( addConstraint_ ) {
		for(int idim=0; idim<ndim_; ++idim) {
			if( floatingConstraint_ ) { 
				minimizer_->SetLimitedVariable(ndim_*nbound + idim, Form( "lambdaDim%d", idim ), 
							       minConstraint_, minConstraint_*0.1,
							       minConstraint_, 1e+3*minConstraint_);

			} else { 
				minimizer_->SetFixedVariable(ndim_*nbound + idim, Form( "lambdaDim%d", idim ), 
							     minConstraint_);
			}
		}
	}
	
	for(size_t iortho=0; iortho<orthocuts_.size(); ++iortho) {
		std::pair<std::string, std::vector<double> >& orthocut = orthocuts_[iortho];
		if( orthocut.second.size() == 1 ) { 
			minimizer_->SetFixedVariable( ndim_*(nbound+addConstraint_) + iortho, orthocut.first,
						      orthocut.second[0] );
		} else if( orthocut.second.size() == 2 ) { 
			minimizer_->SetVariable( ndim_*(nbound+addConstraint_) + iortho, orthocut.first,
						 orthocut.second[0], orthocut.second[1] );
		} else { 
			assert( orthocut.second.size() == 4 );
			minimizer_->SetLimitedVariable( ndim_*(nbound+addConstraint_) + iortho, orthocut.first,
							orthocut.second[0], orthocut.second[1], orthocut.second[2], 
							orthocut.second[3]
				);
		}
	}

	// Call to the minimization
	minimizer_->PrintResults();
	if( ! dryrun ) {
		minimizer_->Minimize();
		minimizer_->PrintResults();
	}

	///// // Refit last boundaries
	///// if( refitLast_ ) {
	///// 	std::vector<double> step0(minimizer_->X(),minimizer_->X()+minimizer_->NDim());
	///// 	int nrefit = 1;
	///// 	
	///// 	for(int idim=0; idim<ndim_; ++idim) {
	///// 		double range = sigModels_[0]->getMax(idim) - sigModels_[0]->getMin(idim);
	///// 		minimizer_->SetFixedVariable(idim*nbound, Form("fixedBoundDim%d",idim), 
	///// 					     step0[idim*nbound]);
	///// 		for(int icat=0; icat<ncat-nrefit; ++icat) {
	///// 			minimizer_->SetFixedVariable(idim*nbound + icat +1, 
	///// 						     Form( (telescopicBoundaries_?"deltaBoundDim%dBin%d":"absBoundDim%dBin%d"),
	///// 							   idim, icat ),
	///// 						     step0[idim*nbound + icat +1]);
	///// 		}
	///// 		
	///// 		for(int icat=ncat-nrefit; icat<ncat; ++icat) {
	///// 			if( telescopicBoundaries_ ) {
	///// 				minimizer_->SetLimitedVariable(idim*nbound + icat +1,
	///// 							       Form( "deltaBoundDim%dBin%d", idim, icat ), 
	///// 							       cutoffs[idim]*(ncat-icat+1.)*10.,   // FIXME smarter initialization
	///// 							       cutoffs[idim]*0.5, cutoffs[idim], range );
	///// 			} else {
	///// 				minimizer_->SetLimitedVariable(idim*nbound + icat +1,
	///// 							       Form( "absBoundDim%dBin%d", idim, icat ), 
	///// 							       sigModels_[0]->getMax(idim) - range/(double)ncat, 
	///// 							       cutoffs[idim]*0.5, 
	///// 							       sigModels_[0]->getMin(idim)+ncat*cutoffs[idim], 
	///// 							       sigModels_[0]->getMax(idim) );
	///// 			}
	///// 
	///// 		}
	///// 	}
	///// 	// penalty to constrain the lower boundaires
	///// 	if( addConstraint_ ) {
	///// 		for(int idim=0; idim<ndim_; ++idim) {
	///// 			if( floatingConstraint_ ) { 
	///// 				minimizer_->SetLimitedVariable(ndim_*nbound + idim, Form( "lambdaDim%d", idim ), 
	///// 							       minConstraint_, minConstraint_*0.1,
	///// 							       minConstraint_, 1e+3*minConstraint_);
	///// 				
	///// 			} else { 
	///// 				minimizer_->SetFixedVariable(ndim_*nbound + idim, Form( "lambdaDim%d", idim ), 
	///// 							     minConstraint_);
	///// 			}
	///// 		}
	///// 	}
	///// 	
	///// 	minimizer_->PrintResults();
	///// 	minimizer_->Minimize();
	///// 	minimizer_->PrintResults();
	///// }
	
	// store results
	std::vector<double> bestFit(minimizer_->X(),minimizer_->X()+minimizer_->NDim());
	double best = minimizer_->MinValue();
	minima_[ncat] = std::make_pair(best,bestFit);
	
	if( debug ) {
		theFom.debug();
		theFom.DoEval(&bestFit[0]);
		theFom.debug(false);
	}
	
	return best;
}

// ---------------------------------------------------------------------------------------------
void CategoryOptimizer::reduce(int ninput, const double * boundaries, const double * cutoffs, int ntarget, double threshold)
{
	assert( ndim_ == 1 && orthocuts_.size() == 0 ); // multi-dim case not yet implemented

	GenericFigureOfMerit startingFom(sigModels_,bkgModels_,fom_,ndim_,ninput,cutoffs,0,false,false,
					 std::vector<HistoConverter *>(0));
	double f0 = startingFom.DoEval(boundaries);
	double fn = f0;	
	std::vector<double> bn(boundaries,boundaries+ninput);
	std::ostream_iterator< double > output( std::cout, "," );
	while( (bn.size() > ntarget) && (fabs( (f0-fn)/fn ) < threshold) ) {
		if( fn < minima_[bn.size()-1].first ) {
			minima_[bn.size()-1] = std::make_pair(fn,bn);
		}
		if( bn.size() < 3 ) { break; }
		std::vector<double> bnm1(bn.size()-1);
		float fnm1;
		for(size_t itest=0; itest<bn.size(); ++itest) {
			std::vector<double> btest(bn.size()-1);
			if( itest > 0 ) {
				std::copy(bn.begin(),bn.begin()+itest,btest.begin());
			}
			std::copy(bn.begin()+itest+1,bn.end(),btest.begin()+itest);
			GenericFigureOfMerit nm1Fom(sigModels_,bkgModels_,fom_,ndim_,btest.size(),cutoffs,0,false,false,
						    std::vector<HistoConverter *>(0));
			double ftest = nm1Fom(&btest[0],0);
			if( itest == 0 || ftest < fnm1 ) {
				fnm1 = ftest;
				bnm1.swap(btest);
			}
		}
		bn.swap(bnm1);
		fn = fnm1;
	}
}


// ---------------------------------------------------------------------------------------------
double CategoryOptimizer::getBoundaries(int ncat, double * boundaries)
{
	int nbound = ncat+1;
	std::cout << "get Boundaries " << ncat << std::endl;
	std::pair<double,std::vector<double> > & res = minima_[ncat];
	std::copy(res.second.begin(),res.second.begin()+ndim_*nbound,boundaries);
	if( telescopicBoundaries_ ) {
		for(int idim=0; idim<ndim_; ++idim) {
			for(int icat=1; icat<=ncat; ++icat) {
				boundaries[idim*ncat+icat] = boundaries[idim*ncat+icat-1] - boundaries[idim*ncat+icat];
			}
		}
	}
	for(int idim=0; idim<ndim_; ++idim) {
		for(int icat=0; icat<=ncat; ++icat) {
			if( !transformations_.empty() && transformations_[idim]!=0 ) {
				boundaries[idim*ncat+icat] = transformations_[idim]->eval(boundaries[idim*ncat+icat]);
			}
		}
	}
	
	return res.first;

}
