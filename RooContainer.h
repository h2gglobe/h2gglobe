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
    void SaveSystematicsData(bool _do=true);
    void MakeSystematicPdfs(bool _do=true);
    void MakeSystematicStudy(std::vector<std::string>,std::vector<int>);
    void AddObservable(std::string,double,double);
    void AddConstant(std::string,double);
    void AddRealVar(std::string,double,double xmin=-10,double xmax=10);
    void AddFormulaVar(std::string,std::string,std::string);

   /** this is used (e.g. by StatAnalysis::Init(..)) to define
       the functional form of the background PDFs to be fitted to
       the data later on. 

       This calls addGenericPdf(..) (with a lower case a !) for
       each defined category.

   */
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

   /** calls createsDataSet(..) (with lower case c !) for all categories */
   void CreateDataSet(std::string name, std::string data_name, int nbins,double x1=-990,double x2=-990); 

   void MakeSystematics(std::string,std::string,int);

   //----------

   /** calls fitToData(..) (note the lower case first letter) for each category */
   void FitToData(std::string name_func, std::string name_var
	         ,double x1,double x2,double x3,double x4);

   void FitToData(std::string name_func, std::string name_var,double x1,double x2);

   /** Calls the most generic FitToData(..) method with some parameters
       set to default values.

       This method is for example called by StatAnalysis::Term(..) which in turn is called
       by LoopAll::TermReal(..) */
   void FitToData(std::string name_func,std::string name_var); 

   //----------

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
   bool fit_systematics;
   bool verbosity_;

  private:

   void addRealVar(std::string,double,double);
   void addRealVar(std::string,double,double,double);
   void addFormulaVar(std::string,std::string,std::string);

   /** this is called by AddGenericPdf(..) e.g. to define
       the functional form of the background parametrization
       (to be fitted to the data).

       @param form determines which functional form should be used
       (see also the comments in the implementation)
       e.g. 70.. corresponds to RooBernstein with form-70 parameters.

       @param obs_name is the name of the variable 
       @param var is the list of variable names to be used
         as parameters for the generated pdf

       @param formula only used for certain values of 'form'
   */
   void addGenericPdf(std::string name, std::string formula, std::string obs_name,
		      std::vector<std::string> &var,
		      int form, 
		      double norm_guess, double norm_min, double norm_max);
   void composePdf(std::string , std::string 
			     ,std::vector<std::string> &,bool);
  // void convolutePdf(std::string,std::string,std::string,RooRealVar &,double norm_guess=100);

   void sumBinnedDatasets(std::string,std::string,std::string,double,double,bool);

   /** creates a RooDataSet for all defined categories 

       @param data_name is the name and title of the RooDataSet to be created.
       @param name is the name of the variable 

       Note that there is another function CreateDataSet(..) which has a capital C but
       otherwise has the same list of parameters.
   */
   void createDataSet(std::string name, std::string data_name, int nbins, double x1, double x2);

   void makeSystematics(std::string,std::string,int);
   
   /** the most generic fitToData function, containing the implementation
       of the fit of the previously defined pdfs to the data */
   void fitToData(std::string name_func, std::string name_data, std::string name_var
	         ,double x1, double x2, double x3, double x4);

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

   /** the list of already known variables. Maps from the variable's
       name to the actual RooRealVar object */
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
