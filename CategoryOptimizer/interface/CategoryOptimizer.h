#ifndef _CategoryOptimizer_h_
#define _CategoryOptimizer_h_

#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include "RooAbsPdf.h"
#include <vector>
#include <map>

#include "FunctionHelpers.h"

// ------------------------------------------------------------------------------------------------
class AbsModel 
{
public:
	enum type_t { sig=0, bkg };
	
	type_t getType() { return type_; };
	
	virtual RooAbsPdf * getCategoryPdf(int icat) { return 0; };
	virtual RooRealVar * getX() { return 0; };
	virtual double getCategoryYield(int icat) { return categoryYields_[icat]; };
	
	virtual int getNcat() { return categoryYields_.size(); };
	
	virtual void clear() { categoryYields_.clear(); };

	virtual std::string name() const { return ""; };
	
protected:
	type_t type_;
	std::vector<double> categoryYields_;
};

// ------------------------------------------------------------------------------------------------
class AbsModelBuilder
{
public:
	virtual AbsModel * getModel() = 0;
	virtual void setOrthoCuts(double * cuts) { assert(0); };
	virtual void beginIntegration(double * theta) = 0;
	virtual void endIntegration() = 0;
	
	virtual bool addBoundary(double * theta) = 0;
	
	virtual double getMin(int idim) = 0;
	virtual double getMax(int idim) = 0;

	virtual TH1 * getPdf(int idim) {};
	virtual double getPenalty() { return 1e+10; };
};

// ------------------------------------------------------------------------------------------------
class AbsFomProvider
{
public:
	AbsFomProvider() : debug_(false) {};

	virtual double operator() ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const = 0;
	
	void debug(bool x) { debug_ = x; }

protected:
	bool debug_;
	
};


// ------------------------------------------------------------------------------------------------
class GenericFigureOfMerit : public ROOT::Math::IBaseFunctionMultiDim
{
public:
	GenericFigureOfMerit(std::vector<AbsModelBuilder *> & sig, std::vector<AbsModelBuilder *> & bkg, AbsFomProvider * fom, 
			     int ndim, int nbound, const double * cutoffs, int northocuts, 
			     bool addConstraint, bool telescopicBoundaries, const std::vector<HistoConverter *> & transformations);
	
	double operator() (double *x, double *p) const;

	virtual double DoEval(const double * x) const; 
	
	virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const { return new GenericFigureOfMerit(*this); }
	
	virtual unsigned int NDim() const { return ndim_*(nbound_+(addConstraint_?1:0))+northocuts_; };
	
	void debug(bool x=true) { fom_->debug(x); };
	
private:
	std::vector<AbsModelBuilder *> sigModels_, bkgModels_, allModels_;
	AbsFomProvider * fom_;
	
	int ndim_;
	int nbound_;
	int northocuts_;
	std::vector<double> cutoffs_;
	
	bool addConstraint_, telescopicBoundaries_;
	const std::vector<HistoConverter *> & transformations_;
};

// ------------------------------------------------------------------------------------------------
class CategoryOptimizer
{
public:
	CategoryOptimizer( ROOT::Math::Minimizer * minimizer, int ndim) : 
		minimizer_(minimizer), ndim_(ndim), strategy_(2), scan_(-1), scanBoundaries_(true), tranformOrtho_(false),
		addConstraint_(false), telescopicBoundaries_(true), floatFirst_(false), 
		refitLast_(false), speed_(0.5), transformations_(0), inv_transformations_(0), dimnames_(ndim_) {};
	
	void addSignal(AbsModelBuilder * sig, bool defineTransform=false) { 
		sigModels_.push_back(sig); 
		if( defineTransform ) { transformModels_.push_back(sig); }
	};
	void addBackground(AbsModelBuilder * bkg, bool defineTransform=false) { 
		bkgModels_.push_back(bkg); 
		if( defineTransform ) { transformModels_.push_back(bkg); }
	};
	void addConstraint(bool x=true, double minConstraint=10., bool floating=true) { 
		addConstraint_ = x; minConstraint_ = minConstraint; floatingConstraint_=floating;
	};
	void floatFirst(bool x=true) { floatFirst_ = x; };
	void refitLast(bool x=true) { refitLast_ = x; };
	void absoluteBoundaries(bool x=true) { telescopicBoundaries_ = !x; };
	void setFigureOfMerit(AbsFomProvider * fom) { fom_ = fom; };
	
	double optimizeNCat(int ncat, const double * cutoffs, bool dryrun=false, bool debug=false, const double * initial_values=0);
	double getBoundaries(int ncat, double * boundaries, double * orthocuts);

	void reduce(int ninput, const double * boundaries, const double * cutoffs, int ntarget=1, double threshold=1.);

	void addFloatingOrthoCut(const char * name, double val, double step, double mix=-1., double max=-1.);
	void addFixedOrthoCut(const char * name, double val);
	void setOrthoCut(int icut, double x) { orthocuts_[icut].second[0] = x; };
	int nOrthoCuts() { return orthocuts_.size(); };

	void setTransformation(int idim, HistoConverter * x, HistoConverter *y) { 
		transformations_.resize(ndim_+(tranformOrtho_?orthocuts_.size():0),0);
		transformations_[idim] = x; 
		inv_transformations_.resize(transformations_.size(),0);
		inv_transformations_[idim] = y; 
	};
	static void doTransform(const std::vector<HistoConverter *> & transformations, double* boundaries) {
		for(size_t ii=0; ii<transformations.size(); ++ii) {
			if( transformations[ii] ) { boundaries[ii] = (*transformations[ii])(&boundaries[ii],0); }
		}
	};

	void setStrategy(int x) { strategy_=x; };
	void setScan(int x, bool y, int repeat=1) { scan_=x; scanBoundaries_=y; repeat_=repeat;};
	void setTransformOrtho(bool x=true) { tranformOrtho_=x; };
	void setSpeed(double x) { speed_=x; };
	
	void setDimName(int idim, const char * name) { dimnames_[idim] = name; };
private:

	ROOT::Math::Minimizer * minimizer_;
	int ndim_, strategy_, scan_, repeat_;
	bool scanBoundaries_, tranformOrtho_;	
	
	std::vector<AbsModelBuilder *> sigModels_;
	std::vector<AbsModelBuilder *> bkgModels_;
	std::vector<AbsModelBuilder *> transformModels_;
	AbsFomProvider * fom_;
	
	std::map<int, std::pair<double,std::vector<double> > > minima_;
	
	bool addConstraint_, telescopicBoundaries_, floatingConstraint_, floatFirst_, refitLast_;
	double minConstraint_, speed_;
	std::vector<std::pair<std::string, std::vector<double> > > orthocuts_;
	
	std::vector<HistoConverter *> transformations_, inv_transformations_;
	std::vector<TString> dimnames_;
	
};

#endif
