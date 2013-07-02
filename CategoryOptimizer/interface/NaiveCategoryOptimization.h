#ifndef _NaiveCategoryOptimization_h_
#define _NaiveCategoryOptimization_h_

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/IFunction.h"

#include <algorithm>

#include "CategoryOptimizer/interface/CategoryOptimizer.h"
#include "CategoryOptimizer/interface/FunctionHelpers.h"

class CutAndCountModelBuiler;

// ------------------------------------------------------------------------------------------------
class CutAndCountModel : public AbsModel 
{
public:
	CutAndCountModel(AbsModel::type_t type=AbsModel::sig) { type_ = type; };
	
	void addCategory(double counts) { categoryYields_.push_back(counts); };
};

// ------------------------------------------------------------------------------------------------
class CutAndCountModelBuilder : public AbsModelBuilder
{
public:
	CutAndCountModelBuilder(AbsModel::type_t type, TH1 * pdf, double norm, double min,  double max) :
		norm_(norm), 
		converter_(new HistoToTF1(Form("%s_integrator",pdf->GetName()),integrate1D(pdf))),
		model_(type) {
		ranges_.push_back(std::make_pair(min,max));
	};
	
	CutAndCountModelBuilder(AbsModel::type_t type, TH2 * pdf, double norm, double xmin, double xmax, 
			       double ymin, double ymax)  : 
		norm_(norm),
		converter_(new SimpleHistoToTF2(Form("%s_integrator",pdf->GetName()),integrate2D(pdf))),
		model_(type) {
		ranges_.push_back(std::make_pair(xmin,xmax));
		ranges_.push_back(std::make_pair(ymin,ymax));
	};

	~CutAndCountModelBuilder() {
		delete converter_;
	};
	
	AbsModel * getModel() { return &model_; };
	void beginIntegration(double * boundaries) { 
		lastIntegral_ = (*converter_)(boundaries,0);
		model_.clear(); 
	};

	void endIntegration() {};
	
	bool addBoundary(double * boundaries) {
		bool ret = true;
		double integral = (*converter_)(boundaries,0);
		double norm = norm_*(integral - lastIntegral_);
		if( norm == 0 ) { ret = false; }
		model_.addCategory(norm);
		lastIntegral_ = integral;
		return ret;
	};
	
	double getMin(int idim) { return ranges_[idim].first;  };
	double getMax(int idim) { return ranges_[idim].second; };
	
	HistoConverter * getInputModel() { return converter_; };
	
private:
	CutAndCountModel model_;
	double norm_, lastIntegral_;
	std::vector<std::pair<double,double> > ranges_;
	
	HistoConverter * converter_;
};


// ------------------------------------------------------------------------------------------------
class NaiveCutAndCountFomProvider : public AbsFomProvider 
{
public:
	NaiveCutAndCountFomProvider(ROOT::Math::IBaseFunctionMultiDim * denom=0) : own_(false), denom_(denom) {};
	NaiveCutAndCountFomProvider(TF1 * denom) : own_(true), denom_(new TF1ToFunctor(denom)) {};
	
	~NaiveCutAndCountFomProvider() {
		if( own_ ) { delete denom_; }
	};
	
	double operator()  ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const;

private:
	bool own_;
	ROOT::Math::IBaseFunctionMultiDim * denom_;

};

// ------------------------------------------------------------------------------------------------
class PoissonCutAndCountFomProvider : public AbsFomProvider 
{
public:
	double operator()  ( std::vector<AbsModel *> sig, std::vector<AbsModel *> bkg) const;
};

#endif
