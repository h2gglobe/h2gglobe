#ifndef ROOCONTAINER
#define ROOCONTAINER

// ROOT includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

// RooFit includes
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
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

    RooContainer(int ncat=1, int nsigmas=1);
    
    ~RooContainer(){};
    void SetNCategories(int);
    void Verobose(bool noisy=true);
    void AddGlobalSystematic(std::string,double,double);
    void SaveSystematicsData();
    void MakeSystematicStudy(std::vector<std::string>,std::vector<int>);
    void AddObservable(std::string,double,double);
    void AddConstant(std::string,double);
    void AddRealVar(std::string,double,double xmin=-10,double xmax=10);
    void AddFormulaVar(std::string,std::string,std::string);
    void AddGenericPdf(std::string,std::string,std::string,
		       std::vector<std::string> &, 
		       int form,double norm_guess=10
		       ,double norm_min=0., double norm_max=1.e6);
    void AddSpecificCategoryPdf(int *,std::string,std::string,std::string,
		       std::vector<std::string> &, 
		       int form,double norm_guess=10
		       ,double norm_min=0., double norm_max=1.e6);
    void ComposePdf(std::string, std::string
		    ,std::vector<std::string> &,bool use_extended=true);
    void ComposeSpecificCategoryPdf(int *,std::string, std::string
		    ,std::vector<std::string> &,bool use_extended=true);
   // void ConvolutePdf(std::string,std::string,std::string,std::string
   //		      ,double norm_guess=100);

   void SumBinnedDatasets(std::string,std::string,std::string,std::vector<double>, std::vector<double>, bool scale=true);
   void CreateDataSet(std::string,std::string,int nbins,double x1=-990,double x2=-990); 
   void MakeSystematics(std::string,std::string,int);

   void FitToData(std::string,std::string 			
	         ,double x1,double x2,double x3,double x4);
   void FitToData(std::string,std::string,double x1,double x2);
   void FitToData(std::string,std::string); 

   void FitToSystematicSet(std::string,std::string,std::string
	                  ,double x1,double x2,double x3,double x4);
   void FitToSystematicSet(std::string,std::string,std::string,double x1,double x2);
   void FitToSystematicSet(std::string,std::string,std::string);

   void InputDataPoint(std::string,int,double,double w=1.);
   void InputSystematicSet(std::string s_name, std::string sys_name, std::vector<int> cats
			   ,std::vector<double> x, std::vector<double> weights=std::vector<double>(0));

   void WriteDataCard(std::string,std::string,std::string,std::string);
   void WriteSpecificCategoryDataCards(std::string,std::string,std::string,std::string);
   void GenerateBinnedPdf(std::string,std::string,std::string,int,int,int,double x1=-999,double x2=-999);
   void CombineBinnedDatasets(std::string,std::string, double fraction=-1);
   std::vector<double> GetFitNormalisations(std::string,std::string,double,double);

   void Save();

   void AppendDataSet(std::string, RooDataSet*);
   void AppendTH1F(std::string, TH1F*);

   std::vector<std::string> GetDataSetNames();
   std::vector<std::string> GetTH1FNames();
   
   int ncat;
   int nsigmas;
   float sigmaRange;
   bool make_systematics;
   bool save_systematics_data;
   bool verbosity_;

  private:

   void addRealVar(std::string,double,double);
   void addRealVar(std::string,double,double,double);
   void addFormulaVar(std::string,std::string,std::string);
   void addGenericPdf(std::string,std::string,std::string,
		      std::vector<std::string> &,
		      int form, 
		      double,double, double);
   void composePdf(std::string , std::string 
			     ,std::vector<std::string> &,bool);
  // void convolutePdf(std::string,std::string,std::string,RooRealVar &,double norm_guess=100);

   void sumBinnedDatasets(std::string,std::string,std::string,double,double,bool);
   void createDataSet(std::string,std::string,int,double x1,double x2);
   void makeSystematics(std::string,std::string,int);
   
   void fitToData(std::string,std::string,std::string
	         ,double,double,double,double);
   void fitToSystematicSet(std::string,std::string,std::string
	     ,double,double,double,double);
   void generateBinnedPdf(int,std::string,std::string,std::string,std::string,RooRealVar*,RooDataSet&,int,int,int,double,double);
   void combineBinnedDatasets(std::string,std::string,double);
   void writeRooDataHist(std::string, TH1F *);
   void writeRooPlot(RooPlot *,double);
   void writeSpecificCategoryDataCard(int,std::string,std::string,std::string,std::string);
   void removeDuplicateElements(std::vector<RooAbsPdf*> &);
   void setAllParametersConstant();

   double getNormalisationFromFit(std::string,std::string,RooAbsPdf *,RooRealVar*,double,double,bool);

   void getArgSetParameters(RooArgSet*,std::vector<double> &);
   void setArgSetParameters(RooArgSet*,std::vector<double> &);

   std::map<std::string,int> systematics_;
   std::map<std::string,std::pair<double,double> > global_systematics_;
   std::map<std::string,int>::iterator it_sys;

   std::string getweightName(std::string);   
   std::string getnormName(std::string);   
   std::string getcatName(std::string,int);   
   std::string getsysName(std::string,std::string);   
   std::string getsysindexName(std::string,std::string
			 ,int,int);
   std::vector<RooAbsPdf*> pdf_saves_;

   std::map<std::string, RooRealVar> m_real_var_;
   std::map<std::string, RooFormulaVar> m_form_var_;
   std::map<std::string, RooAbsPdf*> m_gen_;
   std::map<std::string, RooExtendPdf> m_exp_;
   std::map<std::string, RooAddPdf> m_pdf_;
   std::map<std::string, std::vector<RooRealVar*> > m_comp_pdf_norm_;
   std::map<std::string, TH1F> m_th1f_;

   std::map<std::string, double> m_var_min_;
   std::map<std::string, double> m_var_max_;

   std::map<std::string,RooDataSet> data_; 
   std::map<std::string,std::string> data_obs_names_; 
   std::map<std::string,RooRealVar*> m_data_var_ptr_; 
   std::map<std::string,RooRealVar*> m_weight_var_ptr_; 

   std::map<std::string,std::vector<RooDataSet*> > data_up_; 
   std::map<std::string,std::vector<RooDataSet*> > data_dn_;

   std::map<std::string,std::vector<TH1F*> > m_th1f_up_; 
   std::map<std::string,std::vector<TH1F*> > m_th1f_dn_;
 
   std::map<std::string,int> bins_;
   std::map<std::string,double> inits_;
   std::map<RooPlot*,double> fit_res_;
   std::map<std::string,RooFitResult*> fit_results_;

   RooWorkspace ws;   
   

};


#endif
