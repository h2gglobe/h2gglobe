#ifndef ROOCONTAINER
#define ROOCONTAINER

// ROOT includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TAxis.h"

// RooFit includes
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"

// RooStats includes
#include "RooWorkspace.h"

// standard includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>

class RooContainer {

  public:

   RooContainer(int ncat=1);
   ~RooContainer(){};
   void SetNCategories(int);
   void MakeSystematicStudy(std::vector<std::string>);
   void AddObservable(std::string,float,float);
   void AddRealVar(std::string,float,float,float);
   void AddGenericPdf(std::string,std::string,std::string,
		      std::vector<std::string> &, 
		      int form,double norm_guess=10
		     ,double norm_min=0., double norm_max=1.e6);
   void ComposePdf(std::string, std::string
			     ,std::vector<std::string> &,bool use_extended=true);
   void ConvolutePdf(std::string,std::string,std::string,std::string
			     ,double norm_guess=100);

   void CreateDataSet(std::string,std::string,int nbins,double x1=-990,double x2=-990); 
   void MakeSystematics(std::string,std::string,std::string);

   void FitToData(std::string,std::string 			
	         ,double x1,double x2,double x3,double x4);
   void FitToData(std::string,std::string,double x1,double x2);
   void FitToData(std::string,std::string); 

   void FitToSystematicSet(std::string,std::string,std::string
	                  ,double x1,double x2,double x3,double x4);
   void FitToSystematicSet(std::string,std::string,std::string,double x1,double x2);
   void FitToSystematicSet(std::string,std::string,std::string);

   void InputDataPoint(std::string,int,float,float w=1.);
   void InputSystematicSet(std::string,std::string,int
		,std::vector<float>,float w=1.);

   void WriteDataCard(std::string,std::string,std::string,std::string);
   void GenerateBinnedPdf(std::string,std::string,std::string,int,double x1=-999,double x2=-999);


   void Save();
   
  private:

   void addRealVar(std::string,float,float);
   void addRealVar(std::string,float,float,float);
   void addGenericPdf(std::string,std::string,std::string,
		      std::vector<std::string> &,
		      int form, 
		      double,double, double);
   void composePdf(std::string , std::string 
			     ,std::vector<std::string> &,bool);
   void convolutePdf(std::string,std::string,std::string,RooRealVar &,double norm_guess=100);

   void createDataSet(std::string,std::string,int,double x1,double x2);
   void makeSystematics(std::string,std::string,std::string);
   
   void fitToData(std::string,std::string,std::string
	         ,double,double,double,double);
   void fitToSystematicSet(std::string,std::string,std::string
	     ,double,double,double,double);
   void generateBinnedPdf(std::string,std::string,std::string,RooRealVar&,int,double,double);
   void writeRooDataHist(std::string, TH1F *);
   void writeRooPlot(RooPlot *);
   void removeDuplicateElements(std::vector<RooAbsPdf*> &);

   std::vector<std::string> systematics_;
   std::vector<std::string>::iterator it_sys;

   std::string getweightName(std::string);   
   std::string getnormName(std::string);   
   std::string getcatName(std::string,int);   
   std::string getsysName(std::string,std::string);   
   std::string getsysindexName(std::string,std::string
			 ,int,int);
   int ncat;
   int nsigmas;
   bool make_systematics;
   std::vector<RooAbsPdf*> pdf_saves_;

   std::map<std::string, RooRealVar> m_real_var_;
   std::map<std::string, RooAbsPdf*> m_gen_;
   std::map<std::string, RooExtendPdf> m_exp_;
   std::map<std::string, RooAddPdf> m_pdf_;
   std::map<std::string, TH1F> m_th1f_;

   std::map<std::string, float> m_var_min_;
   std::map<std::string, float> m_var_max_;

   std::map<std::string,RooDataSet> data_; 
   std::map<std::string,RooRealVar*> m_data_var_ptr_; 
   std::map<std::string,RooRealVar*> m_weight_var_ptr_; 

   std::map<std::string,std::vector<RooDataSet*> > data_up_; 
   std::map<std::string,std::vector<RooDataSet*> > data_dn_;

   std::map<std::string,std::vector<TH1F*> > m_th1f_up_; 
   std::map<std::string,std::vector<TH1F*> > m_th1f_dn_;
 
   std::map<std::string,int> bins_;
   std::map<std::string,float> inits_;
   std::map<RooPlot*,RooFitResult*> fit_res_;

   RooWorkspace ws;   
   

};


#endif
