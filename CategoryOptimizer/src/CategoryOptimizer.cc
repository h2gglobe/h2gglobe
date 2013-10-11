#include "../interface/CategoryOptimizer.h"
#include "TMath.h"
#include "TMinuitMinimizer.h"

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
double GenericFigureOfMerit::DoEval(const double * x) const 
{
	std::vector<double> xv(x,x+ndim_*(nbound_+(addConstraint_?1:0))+northocuts_);
	std::vector<double> pv(cutoffs_);
	//// std::cout << ndim_<< " " << nbound_ << " " << northocuts_ << " " << ndim_*(nbound_+(addConstraint_?1:0))+northocuts_ << std::endl;
	//// std::copy( xv.begin(), xv.end(), std::ostream_iterator<double>(std::cout, ",") );
	//// std::cout << std::endl;
	
	return this->operator()(&xv[0],&pv[0]); 
}

double penalty(double distance)
{
	return 6.52896*(1.-TMath::Erf(100.*(distance-1.+0.01)))/(distance*distance*distance+0.3*0.3*0.3);
}

// ---------------------------------------------------------------------------------------------
double GenericFigureOfMerit::operator() (double *x, double *p) const
{
	static double cached_ = 0.;
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
	
 	///// std::copy( extfirstb.begin(), extfirstb.end(), std::ostream_iterator<double>(std::cout, ",") );
	///// std::cout << std::endl;
	/// std::vector<double> orthocuts(&x[ndim_*(nbound_+(int)addConstraint_)],&[x+ndim_*(nbound_+(int)addConstraint_)+northocuts_]);
	double * orthocuts = &x[ndim_*(nbound_+(addConstraint_?1:0))];
	if( transformations_.size() > ndim_ ) { 
		extfirstb.resize(extfirstb.size()+northocuts_);
		std::copy(&orthocuts[0],&orthocuts[northocuts_],&extfirstb[ndim_]);
		orthocuts = &extfirstb[ndim_];
	}
	// initialize the model builders
	std::ostream_iterator<double> output(std::cout, ",");
	CategoryOptimizer::doTransform(transformations_,&extfirstb[0]);
	
	///// std::cout << ndim_ << " " << nbound_ << " " << ndim_*(nbound_+(addConstraint_?1:0)) << std::endl;
	///// std::copy( &orthocuts[0], &orthocuts[northocuts_], std::ostream_iterator<double>(std::cout, ",") );
	///// std::cout << std::endl;
 	///// std::copy( extfirstb.begin(), extfirstb.end(), std::ostream_iterator<double>(std::cout, ",") );
	///// std::cout << std::endl;
	for(std::vector<AbsModelBuilder *>::const_iterator imod=allModels_.begin(); imod!=allModels_.end(); ++imod ) {
		if( northocuts_ > 0 ) {
			(*imod)->setOrthoCuts(&orthocuts[0]);
		}
		(*imod)->beginIntegration(&extfirstb[0]);
	}
	
	double ret = 0.;
	// loop over the categories
	std::vector<double> currb = lastb;
	std::vector<double> distances;
	int nb = 0;
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
			for(int idim=0; idim<ndim_; ++idim) {
				if( lastb[idim] < newb[idim] ) { found=false; }
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
				ret += penalty(distance);
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
			/// std::copy( newb.begin(), newb.end(), std::ostream_iterator<double>(std::cout, ",") );
			/// std::cout << std::endl;
			for(std::vector<AbsModelBuilder *>::const_iterator imod=allModels_.begin(); imod!=allModels_.end(); ++imod ) {
				if( ! (*imod)->addBoundary(&newb[0]) ) { 
					std::cout << "adding penalty from model " 
						  << (*imod)->getModel()->getType() 
						  << " " << (*imod)->getModel()->getNcat() 
						  << " " << (*imod)->getModel()->getCategoryYield( 
							  (*imod)->getModel()->getNcat() -1 ) << std::endl;
					ret += penalty((*imod)->getPenalty());
				}
			}
			++nb;
		}
	}
	
	if( nb == 0 ) { ret+=penalty(0.); }
	if( ret > 0 ) { return cached_ + ret; };
	
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
	double fom =  (*fom_)(sigModels,bkgModels);
	if( fom < 0 ) { cached_ = fom; };
	ret += fom;

	// additional constraint
	if( addConstraint_ ) {
		for(int idim=0; idim<ndim_; ++idim ) {
			ret +=  x[ndim_*nbound_+idim]*pow( (lastb[idim]-mindim[idim]), 2. );
		}
	}

	return ret;
}

// ---------------------------------------------------------------------------------------------
double CategoryOptimizer::optimizeNCat(int ncat, const double * cutoffs, bool dryrun, bool debug, 
				       const double * initial_values)
{
	int nbound = ncat+1;
	
	std::vector<double> tmpcutoffs(cutoffs,cutoffs+ndim_);
	bool build = transformations_.empty() && ! transformModels_.empty();
	if( build ) { 
		std::cout << "Buildinig variable transformations" << std::endl; 
	}
	int ndim = ndim_;
	if( tranformOrtho_ ) { ndim += orthocuts_.size(); }
	for(int idim=0; idim<ndim; ++idim) {
		if( build ) {
			TH1 * transformPdf = transformModels_[0]->getPdf(idim);
			transformPdf->Print();
			for(size_t imod=1; imod<transformModels_.size(); ++imod) {
				TH1 * itransformPdf = transformModels_[imod]->getPdf(idim);
				transformPdf->Add(itransformPdf);
				delete itransformPdf;
			}
			transformPdf->Scale(1./transformPdf->Integral());
			/// transformPdf->Print("all");
			float max = transformPdf->GetMaximum();
			if( max > tmpcutoffs[idim] ) { 
				std::cout << "Moving cutofff for dimension " << idim << " " 
					  << tmpcutoffs[idim] << " -> " << max << std::endl;
				tmpcutoffs[idim] = max;
			}
			HistoConverter * conv = cdfInv(transformPdf,
						       transformModels_[0]->getMin(idim),
						       transformModels_[0]->getMax(idim));
			HistoConverter * convm1 = cdf(transformPdf,
						       transformModels_[0]->getMin(idim),
						       transformModels_[0]->getMax(idim));
			setTransformation(idim,conv,convm1);
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
	std::cout << "number of dimensions " << ndim << " " << ndim_ << " " << orthocuts_.size() << " " << tranformOrtho_ << " " 
		  << transformations_.size() << std::endl;
	
	
	// Book the FOM
	GenericFigureOfMerit theFom(sigModels_, bkgModels_, fom_, ndim_, nbound, 
				    &tmpcutoffs[0], orthocuts_.size(),
				    addConstraint_, telescopicBoundaries_, transformations_);
	minimizer_->SetFunction(theFom);
	std::vector<std::pair<int, std::pair<double,double> > > paramsToScan;
	const double * ival = initial_values;
	std::vector<double> bestFit;
	double best = 1.e+6;
	
	// Define category boundaries. 
	// Last boundary fixed to the maximum range in each dimension
	for(int idim=0; idim<ndim_; ++idim) {
		TString dimname = ( dimnames_[idim] != "" ? dimnames_[idim] :  Form("dim%d",idim) );
		double min = sigModels_[0]->getMin(idim);
		double max = sigModels_[0]->getMax(idim);
		if( ! transformations_.empty() && transformations_[idim]!=0 ) {
			min = 0.; max = 1.;
		}
		double first = max;
		if( ! floatFirst_ ) {
			max -= tmpcutoffs[idim];
		}
		double range = (max - min);
		if( ival ) {
			bestFit.push_back(inv_transformations_[idim]->eval(*ival)); ++ival;
		} else { 
			bestFit.push_back(first);
		}
		if( ! floatFirst_ ) { 
			minimizer_->SetFixedVariable(bestFit.size()-1, Form("%sBound%d",dimname.Data(),0),
						     bestFit.back());
		} else {
			minimizer_->SetLimitedVariable(bestFit.size()-1, Form("%sBound%d",dimname.Data(),0), 
						       bestFit.back(), range, // tmpcutoffs[idim]*speed_,
						       min, max );
			if( scan_ > 0 && scanBoundaries_ ) { 
				if( ival ) {
					paramsToScan.push_back(
						std::make_pair(bestFit.size()-1,
							       std::make_pair(inv_transformations_[idim]->eval(*ival),max)));
				} else {
					paramsToScan.push_back(std::make_pair(bestFit.size()-1,std::make_pair(min,max) ));
				}
			}
		}
		for(int ibound=1; ibound<nbound; ++ibound) {
			if( telescopicBoundaries_ ) {
				if( ival ) {
					bestFit.push_back(inv_transformations_[idim]->eval(*ival)); ++ival;
				} else { 
					bestFit.push_back(range/(double)ncat);
				}
				minimizer_->SetLimitedVariable(bestFit.size()-1,
							       Form( "%sDeltaBound%d",dimname.Data(), ibound ), 
							       bestFit.back(),
							       range, // tmpcutoffs[idim]*speed_,
							       tmpcutoffs[idim], range );
				if( scan_ > 0 && scanBoundaries_ ) { 
					paramsToScan.push_back(std::make_pair(bestFit.size()-1,std::make_pair(min,max)));
				}
				
			} else {
				if( ival ) {
					bestFit.push_back(inv_transformations_[idim]->eval(*ival)); ++ival;
				} else { 
					bestFit.push_back(first - (double)(ibound)*range/(double)ncat + tmpcutoffs[idim]);
				}
				minimizer_->SetLimitedVariable(bestFit.size()-1,
							       Form( "%sBound%d",dimname.Data(), ibound ), 
							       bestFit.back(),
							       range, // tmpcutoffs[idim]*speed_, 
							       min, max );
				if( scan_ > 0 && scanBoundaries_ ) { 
					if(ival) { 
						paramsToScan.push_back(
							std::make_pair(bestFit.size()-1,
								       std::make_pair(inv_transformations_[idim]->eval(*ival),
										      bestFit[bestFit.size()-2])));
					} else {
						paramsToScan.push_back(std::make_pair(bestFit.size()-1,std::make_pair(min,max)));
					}
				}
			}
		}
	}
	
	// penalty to constrain the lower boundaires
	if( addConstraint_ ) {
		for(int idim=0; idim<ndim_; ++idim) {
			TString dimname = ( dimnames_[idim] != "" ? dimnames_[idim] :  Form("dim%d",idim) );
			bestFit.push_back(minConstraint_);
			if( floatingConstraint_ ) { 
				minimizer_->SetLimitedVariable(bestFit.size()-1, Form( "%sLambda",dimname.Data() ), 
							       minConstraint_, minConstraint_*0.1,
							       minConstraint_, 1e+3*minConstraint_);
			} else { 
				minimizer_->SetFixedVariable(bestFit.size()-1, Form( "%sLambda",dimname.Data() ), 
							     minConstraint_);
			}
		}
	}
	
	for(size_t iortho=0; iortho<orthocuts_.size(); ++iortho) {
		std::pair<std::string, std::vector<double> >& orthocut = orthocuts_[iortho];
		double min = ( orthocut.second.size() > 2 ? orthocut.second[2] : 0 );
		double max = ( orthocut.second.size() > 3 ? orthocut.second[3] : 0 );
		double start = orthocut.second[0];
		if( transformations_.size() > ndim_ ) {
			min = 0.; max = 1.;
			start = inv_transformations_[ndim_+iortho]->eval(start);
		}
		double step = max-min;/// ( orthocut.second.size() > 1 ? orthocut.second[1] : 0. ) * speed_;
		//// std::cout << orthocut.first << " " << start << " " << step << " " << min << " " << max << std::endl;
		bestFit.push_back(start);

		if( orthocut.second.size() == 1 ) { 
			minimizer_->SetFixedVariable( bestFit.size()-1, orthocut.first,
						      start );
		} else if( orthocut.second.size() == 2 ) { 
			minimizer_->SetVariable( bestFit.size()-1, orthocut.first,
						 start, step );						
		} else { 
			assert( orthocut.second.size() == 4 );
			minimizer_->SetLimitedVariable( bestFit.size()-1, orthocut.first,
							start, step, min, max
				);
			if( scan_ > 0 ) { 
				paramsToScan.push_back(std::make_pair(bestFit.size()-1,
								      std::make_pair(min+0.2*(max-min),max-0.2*(max-min))));
			}
		}
	}
	
	/// std::vector<double> x(scan_), y(+1);
	double x[100], y[100];
	if( scan_>0 ) { 
		for( int irep=0; irep<repeat_; ++irep ) {
			std::cout << "Scanning parameters " << std::endl;
			minimizer_->PrintResults();
			unsigned int nstep = scan_;
			for(int ii=paramsToScan.size()-1; ii>=0; --ii) {
				int ipar = paramsToScan[ii].first;
				std::pair<double,double> rng = paramsToScan[ii].second;
				std::cout << ipar << " " << rng.first << " " << rng.second << std::endl;
				minimizer_->Scan(ipar,nstep,&x[0],&y[0],rng.first,rng.second);
				/// minimizer_->SetVariableValue(ipar, x[ std::min_element(y.begin(),y.end()) - y.begin()]);
				std::copy( &x[0], &x[nstep-1], std::ostream_iterator<double>(std::cout, ",") );
				std::cout << std::endl;
				std::copy( &y[0], &y[nstep-1], std::ostream_iterator<double>(std::cout, ",") );
				std::cout << std::endl;		       
				minimizer_->PrintResults();
			}
		}
	}
	std::cout << "here" << std::endl;

	// Call to the minimization
	std::cout << "Calling minimization (strategy: " << strategy_ << ")" << std::endl;
	std::copy( bestFit.begin(), bestFit.end(), std::ostream_iterator<double>(std::cout, ",") );
	std::cout << std::endl;
	minimizer_->SetStrategy(strategy_);
	minimizer_->PrintResults();
	minimizer_->SetStrategy(strategy_);
	if( ! dryrun ) {
		minimizer_->Minimize();
		const double * x = minimizer_->X();
		std::copy(x,x+minimizer_->NDim(),bestFit.begin());
		std::copy(x,x+minimizer_->NDim(), std::ostream_iterator<double>(std::cout, ",") );
		std::cout << std::endl;
		best = minimizer_->MinValue();
		minimizer_->PrintResults();
	} else {
		best = theFom.DoEval(&bestFit[0]);
	}

	std::copy( bestFit.begin(), bestFit.end(), std::ostream_iterator<double>(std::cout, ",") );
	std::cout << std::endl;
	///// // Refit last boundaries
	///// if( refitLast_ ) {
	///// 	assert( orthocuts_.size() == 0 );
	///// 	TMinuitMinimizer * minuit = dynamic_cast<TMinuitMinimizer*>(minimizer_);
	///// 	assert(minuit != 0);
	///// 
	///// 	std::vector<double> step0(minimizer_->X(),minimizer_->X()+minimizer_->NDim());
	///// 	int nrefit = 1;
	///// 	for(int idim=0; idim<ndim_; ++idim) {
	///// 		double min = sigModels_[0]->getMin(idim);
	///// 		double max = sigModels_[0]->getMax(idim);
	///// 		if( ! transformations_.empty() && transformations_[idim]!=0 ) {
	///// 			min = 0.; max = 1.;
	///// 		}
	///// 		double range = max - min;
	///// 		minuit->FixVariable(idim*nbound);
	///// 		minuit->FixVariable(idim*nbound + icat +1);
	///// 	
	///// 		for(int icat=ncat-nrefit; icat<ncat; ++icat) {
	///// 			if( telescopicBoundaries_ ) {
	///// 				minuit->SetLimitedVariable(idim*nbound + icat +1,
	///// 							       Form( "deltaBoundDim%dBin%d", idim, icat ), 
	///// 							       range/(double)ncat,   // FIXME smarter initialization
	///// 							       tmpcutoffs[idim]*0.5, tmpcutoffs[idim], range );
	///// 			} else {
	///// 				minuit->SetLimitedVariable(idim*nbound + icat +1,
	///// 							       Form( "absBoundDim%dBin%d", idim, icat ), 
	///// 							       step0[idim*nbound + nrefit] - tmpcutoffs[idim]*2.,
	///// 							       tmpcutoffs[idim]*0.5, 
	///// 							       min, max );
	///// 			}
	///// 
	///// 		}
	///// 	}
	///// 	
	///// 	minimizer_->PrintResults();
	///// 	minimizer_->Minimize();
	///// 	minimizer_->PrintResults();
	///// }
	
	// store results
	//// std::vector<double> bestFit(minimizer_->X(),minimizer_->X()+minimizer_->NDim());
	//// double best = minimizer_->MinValue();
	minima_[ncat] = std::make_pair(best,bestFit);
	
	if( debug ) {
		theFom.debug();
		theFom.DoEval(&bestFit[0]);
		theFom.debug(false);
	}
	
	return best;
}

// ---------------------------------------------------------------------------------------------
void CategoryOptimizer::addFloatingOrthoCut(const char * name, double val, double step, double min, double max)
{
	orthocuts_.push_back(std::make_pair(name,std::vector<double>(4)));
	orthocuts_.back().second[0] = val;
	orthocuts_.back().second[1] = step;
	orthocuts_.back().second[2] = min;
	orthocuts_.back().second[3] = max;
	
}

// ---------------------------------------------------------------------------------------------
void CategoryOptimizer::addFixedOrthoCut(const char * name, double val)
{
	orthocuts_.push_back(std::make_pair(name,std::vector<double>(1)));
	orthocuts_.back().second[0] = val;
}

// ---------------------------------------------------------------------------------------------
void CategoryOptimizer::reduce(int ninput, const double * boundaries, const double * cutoffs, int ntarget, double threshold)
{
	assert( ndim_ == 1 ); // multi-dim case not yet implemented

	std::vector<double> inv_boundaries(boundaries,boundaries+ninput);
	for(int idim=0; idim<ndim_; ++idim) {
		for(int ibound=0; ibound<ninput; ++ibound) {
			if( !inv_transformations_.empty() && inv_transformations_[idim]!=0 ) {
				inv_boundaries[idim*ninput+ibound] = 
					inv_transformations_[idim]->eval(inv_boundaries[idim*ninput+ibound]);
			}
		}
	}

	GenericFigureOfMerit startingFom(sigModels_,bkgModels_,fom_,ndim_,ninput,cutoffs,0,false,false,
					 transformations_);
	double f0 = startingFom.DoEval(&inv_boundaries[0]);
	double fn = f0;	
	std::vector<double> bn(inv_boundaries);
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
						    transformations_);
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
double CategoryOptimizer::getBoundaries(int ncat, double * boundaries, double * orthocuts)
{
	int nbound = ncat+1;
	std::cout << "get Boundaries " << ncat << std::endl;
	std::pair<double,std::vector<double> > & res = minima_[ncat];
	if( res.second.empty() ) { return 1.e+6; }
	std::copy(res.second.begin(),res.second.begin()+ndim_*nbound,boundaries);
	if( orthocuts_.size() > 0 ) {
		std::copy(res.second.begin()+ndim_*(nbound+addConstraint_),res.second.begin()+ndim_*(nbound+addConstraint_)+orthocuts_.size(),orthocuts);
		if( transformations_.size() > ndim_ ) {
			for(int iorto=0; iorto<orthocuts_.size(); ++iorto) {
				orthocuts[iorto] = transformations_[ndim_+iorto]->eval(orthocuts[iorto]);
			}
		}
	}
	if( telescopicBoundaries_ ) {
		for(int idim=0; idim<ndim_; ++idim) {
			for(int ibound=1; ibound<nbound; ++ibound) {
				boundaries[idim*nbound+ibound] = boundaries[idim*nbound+ibound-1] - boundaries[idim*nbound+ibound];
			}
		}
	}
	for(int idim=0; idim<ndim_; ++idim) {
		for(int ibound=0; ibound<nbound; ++ibound) {
			if( !transformations_.empty() && transformations_[idim]!=0 ) {
				std::cout << boundaries[idim*nbound+ibound] << " -> ";
				boundaries[idim*nbound+ibound] = transformations_[idim]->eval(boundaries[idim*nbound+ibound]);
				double invb = inv_transformations_[idim]->eval(boundaries[idim*nbound+ibound]);
				std::cout << "(" << invb << ")";
				if( invb >= 1. ) {
					boundaries[idim*nbound+ibound] = sigModels_[0]->getMax(idim);
				} else if( invb <= 0. ) {
					boundaries[idim*nbound+ibound] = sigModels_[0]->getMin(idim);
				}
				std::cout << boundaries[idim*nbound+ibound] << " " << std::endl;
			}
		}
	}
	
	return res.first;
}
