/****************************************************** 
   RooContainer.cc  
   Original Author - Nicholas Wardle - Imperial College
*******************************************************/
   
#include "RooContainer.h"
#include "RooMsgService.h"

// CombinedLimit includes
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace RooFit;

RooContainer::RooContainer(int n, int s):ncat(n),nsigmas(s),make_systematics(false),save_systematics_data(false),verbosity_(false),fit_systematics(false),save_roodatahists(true){
	
// Set up the arrays which may be needed
signalVector1 = new double[25];
backgroundVector1 = new double[25];

blind_data = true;
}

RooContainer::~RooContainer(){
delete [] signalVector1;
delete [] backgroundVector1;
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::SetNCategories(int n){
   ncat = n;
}
void RooContainer::Verbose(bool noisy){
	verbosity_=noisy;
	if (!noisy) RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}
void RooContainer::AddGlobalSystematic(std::string name,double val_sig, double val_bkg){
  global_systematics_[name] = std::pair<double,double>(val_sig,val_bkg);
}
void RooContainer::AddNormalisationSystematics(std::string name,std::vector<std::pair<double,double> > vals, int type){

  if (vals.size() != ncat) {
	std::cerr << "WARNING! -- RooContainer::AddNormalisationsSystematics -- expected input of errors"
		  << "to be equal in size to number of categories=" << ncat << std::endl;
	
  } else {
  	std::vector<std::pair<double,double> > sys;
  	if (type ==-1){ // Signal Systematics
		
		for (int cat=0;cat<ncat;cat++){	
			double value = 1.+(vals[cat].second/vals[cat].first);
			sys.push_back(std::pair<double,double> (value,1.0));
		}
    	} else if (type ==1){ // Background Systematics
		for (int cat=0;cat<ncat;cat++){	
			double value = 1.+(vals[cat].second/vals[cat].first);
			sys.push_back(std::pair<double,double> (1.0,value));
		}
    	} else if (type ==0){ // Both
		for (int cat=0;cat<ncat;cat++){	
			double value = 1.+(vals[cat].second/vals[cat].first);
			sys.push_back(std::pair<double,double> (value,value));
		}
    	} else std::cerr << "WARNING! -- RooContainer::AddNormalisationsSystematics -- type not understood " 
			<< type << std::endl;
 	normalisation_systematics_.insert(std::pair <std::string,std::vector <std::pair< double, double > > > (name,sys));
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::SaveSystematicsData(bool save){
   save_systematics_data = save;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::BlindData(bool blind){
   blind_data = blind;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematicPdfs(bool save){
   fit_systematics = save;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::SaveRooDataHists(bool save){
   save_roodatahists = save;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematicStudy(std::vector<std::string> sys_names,std::vector<int> sys_types){
   make_systematics = true;
   std::vector<std::string>::iterator it = sys_names.begin();
   std::vector<int>::iterator it_typ = sys_types.begin();
   for (;it!=sys_names.end();it++,it_typ++)
     systematics_.insert(std::pair<std::string,int>(*it,*it_typ));
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddConstant(std::string name,double init){
  RooRealVar temp(name.c_str(),name.c_str(),init,init,init);
  ws.import(temp,Silence());
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddObservable(std::string name,double xmin,double xmax){
    addRealVar(name,xmin,xmax);
    m_real_var_[name].setRange("FullObservableRange",xmin,xmax);
    addRealVar(getweightName(name),0.,1.0e6);
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddRealVar(std::string name,double init,double xmin, double xmax){
  for (int cat=0;cat<ncat;cat++){
    
    addRealVar(getcatName(name,cat),init,xmin,xmax);
    if (make_systematics && fit_systematics){
	for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
	  for (int sys=1;sys<=nsigmas;sys++)
	    addRealVar(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),init,xmin,xmax);
	  for (int sys=1;sys<=nsigmas;sys++)
	    addRealVar(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),init,xmin,xmax);
	}
    }
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddFormulaVar(std::string name,std::string formula,std::string var){
  for (int cat=0;cat<ncat;cat++){
    addFormulaVar(getcatName(name,cat),formula,getcatName(var,cat));
    if (make_systematics && fit_systematics){
	for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
	  for (int sys=1;sys<=nsigmas;sys++)
	    addFormulaVar(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),formula,getsysindexName(getcatName(var,cat),(it_sys->first),sys,-1));
	  for (int sys=1;sys<=nsigmas;sys++)
	    addFormulaVar(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),formula,getsysindexName(getcatName(var,cat),(it_sys->first),sys,1));
	}
    }
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddSpecificCategoryPdf(int *categories,std::string name,std::string formula,std::string obs_name
				,std::vector<std::string> & var, int form
				,double norm_guess, double norm_min, double norm_max){
   
  for (int cat=0;cat<ncat;cat++){
    if (!categories[cat]) continue;
    std::vector<std::string> cat_var;
    for (std::vector<std::string>::iterator it=var.begin()
	;it!=var.end()
	;++it){
      cat_var.push_back(getcatName(*it,cat));
    }  
    addGenericPdf(getcatName(name,cat),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);

    if (make_systematics && fit_systematics){
      
      for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
       for (int sys=1;sys<=nsigmas;sys++){
        std::vector<std::string> cat_var;
        for (std::vector<std::string>::iterator it=var.begin()
	  ;it!=var.end()
	  ;++it){
          cat_var.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,-1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
       } 
       for (int sys=1;sys<=nsigmas;sys++){
        std::vector<std::string> cat_var;
        for (std::vector<std::string>::iterator it=var.begin()
	  ;it!=var.end()
	  ;++it){
          cat_var.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
       } 
      }
   }
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddGenericPdf(std::string name,std::string formula,std::string obs_name
				,std::vector<std::string> & var, int form
				,double norm_guess, double norm_min, double norm_max){
  for (int cat=0;cat<ncat;cat++){
    
    std::vector<std::string> cat_var;
    for (std::vector<std::string>::iterator it=var.begin()
	;it!=var.end()
	;++it){
      cat_var.push_back(getcatName(*it,cat));
    }  
    addGenericPdf(getcatName(name,cat),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);

    if (make_systematics && fit_systematics){
      
      for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
       for (int sys=1;sys<=nsigmas;sys++){
        std::vector<std::string> cat_var;
        for (std::vector<std::string>::iterator it=var.begin()
	  ;it!=var.end()
	  ;++it){
          cat_var.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,-1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
       } 
       for (int sys=1;sys<=nsigmas;sys++){
        std::vector<std::string> cat_var;
        for (std::vector<std::string>::iterator it=var.begin()
	  ;it!=var.end()
	  ;++it){
          cat_var.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
       } 
      }
   }
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::ComposeSpecificCategoryPdf(int *categories,std::string name, std::string  composition
			     ,std::vector<std::string> & formula, bool use_extended){

  for (int cat=0;cat<ncat;cat++){
    if (!categories[cat]) continue;

    std::vector<std::string> cat_formula;
    for (std::vector<std::string>::iterator it=formula.begin()
	;it!=formula.end()
	;++it){
      cat_formula.push_back(getcatName(*it,cat));
    }  
    composePdf(getcatName(name,cat),composition,cat_formula,use_extended);
	
    if (make_systematics && fit_systematics){	
     for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
      for (int sys=1;sys<=nsigmas;sys++){
       std::vector<std::string> cat_formula;
       for (std::vector<std::string>::iterator it=formula.begin()
	   ;it!=formula.end()
	   ;++it){
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,-1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),composition,cat_formula,use_extended);
      }
      for (int sys=1;sys<=nsigmas;sys++){
       std::vector<std::string> cat_formula;
       for (std::vector<std::string>::iterator it=formula.begin()
	   ;it!=formula.end()
	   ;++it){
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),composition,cat_formula,use_extended);
      }
     }
   }
 }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::ComposePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula, bool use_extended){
  for (int cat=0;cat<ncat;cat++){
    std::vector<std::string> cat_formula;
    for (std::vector<std::string>::iterator it=formula.begin()
	;it!=formula.end()
	;++it){
      cat_formula.push_back(getcatName(*it,cat));
    }  
    composePdf(getcatName(name,cat),composition,cat_formula,use_extended);
	
    if (make_systematics && fit_systematics){	
     for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
      for (int sys=1;sys<=nsigmas;sys++){
       std::vector<std::string> cat_formula;
       for (std::vector<std::string>::iterator it=formula.begin()
	   ;it!=formula.end()
	   ;++it){
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,-1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),composition,cat_formula,use_extended);
      }
      for (int sys=1;sys<=nsigmas;sys++){
       std::vector<std::string> cat_formula;
       for (std::vector<std::string>::iterator it=formula.begin()
	   ;it!=formula.end()
	   ;++it){
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),(it_sys->first),sys,1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),composition,cat_formula,use_extended);
      }
     }
   }
 }
}

// ----------------------------------------------------------------------------------------------------
/*
void RooContainer::ConvolutePdf(std::string name, std::string f_pdf, std::string g_pdf, std::string observable, double norm_guess){
  std::map<std::string,RooRealVar>::iterator obs = m_real_var_.find(observable);
  if (obs == m_real_var_.end()){
    std::cerr << "WARNING -- RooContainer::ConvolutePdf -- No Observable found named "
	      << observable << std::endl;
    return;
  }

  for (int cat=0;cat<ncat;cat++){
    convolutePdf(getcatName(name,cat),getcatName(f_pdf,cat),getcatName(g_pdf,cat),(*obs).second,norm_guess);
	
    if (make_systematics){	
     for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
      for (int sys=1;sys<=nsigmas;sys++){
       convolutePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,-1),getsysindexName(getcatName(f_pdf,cat),it_sys->first,sys,-1)
		   ,getsysindexName(getcatName(g_pdf,cat),(it_sys->first),sys,-1),(*obs).second,norm_guess);
      }
      for (int sys=1;sys<=nsigmas;sys++){
       convolutePdf(getsysindexName(getcatName(name,cat),(it_sys->first),sys,1),getsysindexName(getcatName(f_pdf,cat),it_sys->first,sys,1)
		   ,getsysindexName(getcatName(g_pdf,cat),(it_sys->first),sys,1),(*obs).second,norm_guess);
      }
     }
   }
 }
}
*/
// ----------------------------------------------------------------------------------------------------
void RooContainer::CreateDataSet(std::string name,std::string data_name,int nbins, double x1, double x2){
  for (int cat=0;cat<ncat;cat++){
    std::string cat_name = getcatName(data_name,cat);
    createDataSet(name,cat_name,nbins,x1,x2);  
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::GenerateBinnedPdf(std::string hist_name,std::string pdf_name,std::string data_name,int effect_type,int nbins,int mode,double x1,double x2){

  for (int cat=0;cat<ncat;cat++){
   std::map<std::string,RooRealVar*>::iterator obs_var = m_data_var_ptr_.find(getcatName(data_name,cat));
   if (obs_var == m_data_var_ptr_.end()) {
    std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  No Dataset named "
	      << data_name << std::endl;
    return;
   }
    std::string obs_name = data_obs_names_[getcatName(data_name,cat)];
    generateBinnedPdf(cat,hist_name,getcatName(hist_name,cat),getcatName(pdf_name,cat),obs_name,obs_var->second,data_[getcatName(data_name,cat)],effect_type,nbins,mode,x1,x2);  
  }
  
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematics(std::string observable,std::string s_name, int effect){

  std::map<std::string,RooRealVar>::iterator obs = m_real_var_.find(observable);
  if (obs == m_real_var_.end()){
    std::cerr << "RooContainer::MakeSystematics -- WARNING!!! -- No Observable named "
	      << observable << std::endl;
    return;
  }

  for (int cat=0;cat<ncat;cat++){
      makeSystematics(observable,getcatName(s_name,cat),effect);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToData(std::string name_func, std::string name_var,double x1, double x2, double x3, double x4){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat),getcatName(name_var,cat)
	     ,x1,x2,x3,x4);
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToData(std::string name_func, std::string name_var,double x1, double x2){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat),getcatName(name_var,cat)
	     ,x1,x2,-999,-999);
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToData(std::string name_func, std::string name_var){
  for (int cat=0;cat<ncat;cat++){
    fitToData(getcatName(name_func,cat),getcatName(name_var,cat),getcatName(name_var,cat)
	     ,-999,-999,-999,-999);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name,double x1,double x2,double x3,double x4){
  for (int cat=0;cat<ncat;cat++) {
    fitToSystematicSet(getcatName(name_func,cat),getcatName(name_var,cat)
		      ,sys_name,x1,x2,x3,x4);
  }    
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name,double x1,double x2){
  for (int cat=0;cat<ncat;cat++) {
    fitToSystematicSet(getcatName(name_func,cat),getcatName(name_var,cat)
		      ,sys_name,x1,x2,-999,-999);
  }    
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::FitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name){
  for (int cat=0;cat<ncat;cat++) {
    fitToSystematicSet(getcatName(name_func,cat),getcatName(name_var,cat)
		      ,sys_name,-999,-999,-999,-999);
  }    
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::Save(){

  std::cout << "RooContainer::Save -- Saving To File "
            << std::endl;

  std::cout << "RooContainer::Save -- Saving Plots "
            << std::endl;

/*
  std::map<RooPlot*,double>::iterator it;
  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
      
       writeRooPlot((*it).first,(*it).second);
  }
*/
  std::vector<TCanvas*>::iterator it;
  for(it = fit_canvases_.begin();it != fit_canvases_.end();it++ ) (*it)->Write();
  
  std::cout << "RooContainer::Save -- Saving Pdfs "
            << std::endl;

  // check first and remove duplicates of the pointers. 
  removeDuplicateElements(pdf_saves_);

  std::vector<RooAbsPdf*>::iterator it_pdf;

  for(it_pdf  = pdf_saves_.begin()
     ;it_pdf != pdf_saves_.end()
     ;it_pdf++ ){
     
	  ws.import(**it_pdf,Silence());
  }
/*
  std::map<std::string,RooAbsPdf*>::iterator it_gen;

  for(it_gen  = m_gen_.begin()
     ;it_gen != m_gen_.end()
     ;it_gen++ ){
     
       ws.import(*(it_gen->second));
  }
  
*/
  for (std::map<std::string,RooDataSet>::iterator it_data = data_.begin()
      ;it_data!=data_.end();it_data++)	{


	ws.import(it_data->second);
  }

  std::cout << "RooContainer::Save -- Saving Histos "
            << std::endl;

  std::map<std::string,TH1F>::iterator it_d = m_th1f_.begin();
  std::map<std::string,TH1F>::iterator it_e = m_th1f_.end();

  for(;it_d != it_e; it_d++){
     writeRooDataHist((*it_d).first,&((*it_d).second));
  }
  ws.SetName("cms_hgg_workspace");
  
  // Make sure all parameters of the pdfs are set constant
  //setAllParametersConstant();

  ws.Write();
}
std::vector< std::pair<double,double> > RooContainer::GetFitNormalisationsAndErrors(std::string pdf_name, std::string data_name, double r1,double r2, bool external_fit){
 
  std::vector< std::pair< double,double> > normalisations;

  for (int cat=0;cat<ncat;cat++){
   std::map<std::string,RooRealVar*>::iterator obs_var = m_data_var_ptr_.find(getcatName(data_name,cat));
   if (obs_var == m_data_var_ptr_.end()) {
    std::cerr << "WARNING!!! -- RooContainer::GetFitNormalisations --  No Dataset named "
	      << data_name << std::endl;
   } else {
     bool multi_pdf;
     RooAbsPdf *pdf_ptr;

     std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(getcatName(pdf_name,cat));

     if (exp != m_exp_.end()) {
      pdf_ptr = &(exp->second);
      multi_pdf = false;
     } else {
      std::map<std::string,RooAddPdf>::iterator pdf  = m_pdf_.find(getcatName(pdf_name,cat));
     if (pdf != m_pdf_.end()) {
      multi_pdf = true;
      pdf_ptr = &(pdf->second);
     }
      else {
        std::cerr << "WARNING!!! -- RooContainer::GetFitNormalisations --  No Pdf named "
	         << pdf_name << std::endl;  
      } 
     }
     std::cout << getcatName(pdf_name,cat)<< std::endl;
     std::string obs_name = data_obs_names_[getcatName(data_name,cat)];
     // Just use a weird name (getcatName(obs_name)) for the histogram since we dont need it
     normalisations.push_back(getNormalisationAndErrorFromFit(getcatName(pdf_name,cat),getcatName(obs_name,cat),pdf_ptr,obs_var->second,r1,r2,multi_pdf,external_fit));
   }
  }

  return normalisations;

}
std::vector<double> RooContainer::GetFitNormalisations(std::string pdf_name, std::string data_name, double r1,double r2, bool external_fit){
 
  std::vector<double> normalisations;

  for (int cat=0;cat<ncat;cat++){
   std::map<std::string,RooRealVar*>::iterator obs_var = m_data_var_ptr_.find(getcatName(data_name,cat));
   if (obs_var == m_data_var_ptr_.end()) {
    std::cerr << "WARNING!!! -- RooContainer::GetFitNormalisations --  No Dataset named "
	      << data_name << std::endl;
   } else {
     bool multi_pdf;
     RooAbsPdf *pdf_ptr;

     std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(getcatName(pdf_name,cat));

     if (exp != m_exp_.end()) {
      pdf_ptr = &(exp->second);
      multi_pdf = false;
     } else {
      std::map<std::string,RooAddPdf>::iterator pdf  = m_pdf_.find(getcatName(pdf_name,cat));
     if (pdf != m_pdf_.end()) {
      multi_pdf = true;
      pdf_ptr = &(pdf->second);
     }
      else {
        std::cerr << "WARNING!!! -- RooContainer::GetFitNormalisations --  No Pdf named "
	         << pdf_name << std::endl;  
      } 
     }
     std::cout << getcatName(pdf_name,cat)<< std::endl;
     std::string obs_name = data_obs_names_[getcatName(data_name,cat)];
     // Just use a weird name (getcatName(obs_name)) for the histogram since we dont need it
     normalisations.push_back(getNormalisationFromFit(getcatName(pdf_name,cat),getcatName(obs_name,cat),pdf_ptr,obs_var->second,r1,r2,multi_pdf,external_fit));
   }
  }

  return normalisations;

}
// ----------------------------------------------------------------------------------------------------
void RooContainer::InputDataPoint(std::string var_name, int cat, double x, double w){

  if (cat < 0) return;
  if (cat>-1 && cat<ncat){
    std::string name = getcatName(var_name,cat);
    std::map<std::string, RooDataSet>::iterator it_var  = data_.find(name);
    if (it_var == data_.end()) 
      std::cerr << "WARNING -- RooContainer::InputDataPoint -- No DataSet named "<< name << std::endl;
    else{
      double min_x = m_var_min_[name];
      double max_x = m_var_max_[name];

      if (x > min_x && x < max_x){
        *(m_data_var_ptr_[name]) = x;
        ((*it_var).second).add(RooArgSet(*(m_data_var_ptr_[name])),w);

        m_th1f_[name].Fill(x,w);
      }
    }
  }

  else {
    std::cerr << "WARNING -- RooContainer::InputDataPoint -- No Category Number " << cat 
              << ", category must be from 0 to " << ncat-1
	      << std::endl;
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::InputBinnedDataPoint(std::string var_name, int cat, double x, double w){
 
  if (cat>-1 && cat<ncat){
    std::string name = getcatName(var_name,cat);
    std::map<std::string, RooDataSet>::iterator it_var  = data_.find(name);
    if (it_var == data_.end()) 
      std::cerr << "WARNING -- RooContainer::InputDataPointBinned -- No DataSet named "<< name << std::endl;
    else{
      double min_x = m_var_min_[name];
      double max_x = m_var_max_[name];

      if (x > min_x && x < max_x){
        m_th1f_[name].Fill(x,w);
      }
    }
  }

/*
  else {
    std::cerr << "WARNING -- RooContainer::InputDataPointBinned -- No Category Number " << cat 
              << ", category must be from 0 to " << ncat-1
	      << std::endl;
  }
*/
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::InputSystematicPoint(std::string s_name, std::string sys_name,int ishift, int cat,double val,double w){

  
  
	if (cat>-1 && cat<ncat) {

	  std::string cat_name = getcatName(s_name,cat);
	  std::string name = getsysindexName( cat_name, sys_name, abs(ishift), (ishift > 0 ? 1 : -1) );
	  
	  std::map<std::string,RooDataSet>::iterator it_var  = data_.find(cat_name);
	  
	  if (it_var == data_.end()) 
	    std::cerr << "WARNING -- RooContainer::InpusSystematicSet -- No DataSet named "<< cat_name << std::endl;
	  
	  else {
	    
	    // Safe to use this since creation of systematic set guaranteed to have come from an already existing dataset
	    RooRealVar *ptr = m_data_var_ptr_[cat_name];
	    double min_x = m_var_min_[cat_name];
	    double max_x = m_var_max_[cat_name];
	    
	    RooDataSet & data_set = data_[name];
	    TH1F & th1f_set = m_th1f_[name];
	    
	    if (val > min_x &&val < max_x ){

		// Only make datasets if the save_systematics_data is on
		if (save_systematics_data){
	    	   *ptr = val;
	    	   data_set.add(RooArgSet(*ptr),w);
		}
	    	th1f_set.Fill(val,w);
	    }
	  }
	  
	}
      
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::InputSystematicSet(std::string s_name, std::string sys_name, std::vector<int> cats
				      ,std::vector<double> x, std::vector<double> weights){

  if( weights.empty() ) { weights.resize(x.size(),1.); }
  
  assert( x.size() == cats.size() );
  assert( x.size() == weights.size() );

  if ((int)x.size() != 2*nsigmas){
    std::cerr << "WARNING -- RooContainer::InputSystematicSet -- Size of vector must be equal to "
	      << 2*nsigmas << " -- WARNING" << std::endl;
    assert(0);
    // return;	
  }

  else{
  
    for(size_t istep=0; istep < x.size(); ++istep ) {
    
      int cat = cats[istep];
      double w = weights[istep];
  
	if (cat>-1 && cat<ncat) {

	  /// int ishift = istep - nsigmas;
	  /// if (ishift==0) ishift++;
	  int ishift = (int) istep < nsigmas ? (int) istep - nsigmas : (int) istep - nsigmas + 1;

	  std::string cat_name = getcatName(s_name,cat);
	  std::string name = getsysindexName( cat_name, sys_name, abs(ishift), (ishift > 0 ? 1 : -1) );
	  
	  std::map<std::string,RooDataSet>::iterator it_var  = data_.find(cat_name);
	  
	  if (it_var == data_.end()) 
	    std::cerr << "WARNING -- RooContainer::InpusSystematicSet -- No DataSet named "<< cat_name << std::endl;
	  
	  else {
	    
	    // Safe to use this since creation of systematic set guaranteed to have come from an already existing dataset
	    RooRealVar *ptr = m_data_var_ptr_[cat_name];
	    double min_x = m_var_min_[cat_name];
	    double max_x = m_var_max_[cat_name];
	    
	    	RooDataSet & data_set = data_[name];
	    	TH1F & th1f_set = m_th1f_[name];
	    
	    	double val = x[istep];
	    if (val > min_x &&val < max_x ){

		// Only make datasets if the save_systematics_data is on
		if (save_systematics_data){
	    	   *ptr = val;
	    	   data_set.add(RooArgSet(*ptr),w);
		}
	    	th1f_set.Fill(val,w);
	    }
	  }
	  
	}
      
	//else {
	//  std::cerr << "WARNING -- RooContainer::InputDataPoint -- No Category Number " << cat 
	//	    << ", category must be from 0 to " << ncat-1
	//	    << std::endl;
	//}
    }
    
    //// else{
    ////   
    ////   if (cat>-1 && cat<ncat){
    ////     
    ////     std::string cat_name = getcatName(s_name,cat);
    ////     std::string name = getsysName(cat_name,sys_name);
    ////     
    ////     std::map<std::string,RooDataSet>::iterator it_var  = data_.find(cat_name);
    ////     
    ////     if (it_var == data_.end()) 
    //// 	std::cerr << "WARNING -- RooContainer::InpusSystematicSet -- No DataSet named "<< cat_name << std::endl;
    ////     
    ////     else {
    //// 	
    //// 	// Safe to use this since creation of systematic set guaranteed to have come from an already existing dataset
    //// 	RooRealVar *ptr = m_data_var_ptr_[cat_name];
    //// 	double min_x = m_var_min_[cat_name];
    //// 	double max_x = m_var_max_[cat_name];
    //// 	
    //// 	std::vector<RooDataSet*>::iterator data_set_up = data_up_[name].begin();
    //// 	std::vector<RooDataSet*>::iterator data_set_dn = data_dn_[name].begin();
    //// 	std::vector<TH1F*>::iterator th1f_set_up = m_th1f_up_[name].begin();
    //// 	std::vector<TH1F*>::iterator th1f_set_dn = m_th1f_dn_[name].begin();
    //// 	
    //// 	std::vector<double>::iterator val = x.begin();
    //// 	
    //// 	
    //// 	// Loop over the first nsigmas elements as the -1 -> dn sys
    //// 	for (;data_set_dn != data_dn_[name].end()
    //// 	       ;data_set_dn++,th1f_set_dn++,val++){
    //// 	  
    //// 	  if (*val > min_x && *val < max_x){
    //// 	    *ptr = *val;
    //// 	    (*data_set_dn)->add(RooArgSet(*ptr),w);
    //// 	    (*th1f_set_dn)->Fill(*val,w);
    //// 	  }
    //// 	}
    //// 	
    //// 	// Loop over the second nsigmas elements as the +1 -> up sys
    //// 	for (;data_set_up != data_up_[name].end()
    //// 	       ;data_set_up++,th1f_set_up++,val++){
    //// 	  
    //// 	  if (*val > min_x && *val < max_x){
    //// 	    *ptr = *val;
    //// 	    (*data_set_up)->add(RooArgSet(*ptr),w);
    //// 	    (*th1f_set_up)->Fill(*val,w);
    //// 	  }
    //// 	}
    ////     }
    ////   }
    ////   else {
    ////     std::cerr << "WARNING -- RooContainer::InputDataPoint -- No Category Number " << cat 
    //// 		<< ", category must be from 0 to " << ncat-1
    //// 		<< std::endl;
    ////   }
    //// }
  }
}
// -----------------------------------------------------------------------------------------
void RooContainer::SumBinnedDatasets(std::string new_name, std::string data_one,std::string data_two, std::vector<double> coefficients_one, std::vector<double> coefficients_two, bool scale){

   if (coefficients_one.size() != ncat || coefficients_two.size() != ncat ){
	std::cerr << "WARNING!! -- RooContainer::SumBinnedDataSets -- number of coefficients should be the same as number of categories " << std::endl;
   } else {
	
	for (int cat=0;cat<ncat;cat++){
	   sumBinnedDatasets(getcatName(new_name,cat),getcatName(data_one,cat),getcatName(data_two,cat),coefficients_one[cat],coefficients_two[cat],scale);
	}
	
   }
}
// -----------------------------------------------------------------------------------------
void RooContainer::SumBinnedDatasets(std::string new_name, std::string data_one,std::string data_two,double coefficients_one,double coefficients_two, bool scale){

	for (int cat=0;cat<ncat;cat++){
	   sumBinnedDatasets(getcatName(new_name,cat),getcatName(data_one,cat),getcatName(data_two,cat),coefficients_one,coefficients_two,scale);
	}
	
}
// -----------------------------------------------------------------------------------------
void RooContainer::MergeHistograms(std::string data_one,std::string data_two, bool systematics){

	// This will take 2 histograms and stick them one after another on the same x-axis, (as long as xmax_1=xmin_2)
	// WARNING- this will not work on Datasets because It replaces the histograms

	for (int cat=0;cat<ncat;cat++){

	  std::string catNameOne = getcatName(data_one,cat);
	  std::string catNameTwo = getcatName(data_two,cat);
	  std::map<std::string,RooDataSet>::iterator it_data = data_.find(catNameOne);

	  if (it_data!=data_.end()){

		std::cerr << "WARNING -- RooContainer::MergeHistograms -- Cannot Merge Original Datasets "
			  << catNameOne
			  << std::endl;
		return;
	  }
	
	  std::map<std::string,TH1F>::iterator itOne = m_th1f_.find(catNameOne);
	  std::map<std::string,TH1F>::iterator itTwo = m_th1f_.find(catNameTwo);

	  if (itOne!=m_th1f_.end() && itTwo!=m_th1f_.end())
		mergeHistograms(catNameOne,&(itOne->second),&(itTwo->second));
	  else {
		std::cerr << "WARNING -- RooContainer::MergeHistograms -- No Histograms found named "
			  << catNameOne << " or " << catNameTwo
			  << std::endl;
		return;
	  }

	  if (systematics){
		for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
		  for (int sys=1;sys<=nsigmas;sys++){
		     std::string sysDNameOne = getsysindexName(catNameOne,it_sys->first,sys,-1);
		     std::string sysUNameOne = getsysindexName(catNameOne,it_sys->first,sys,1);
		     std::string sysDNameTwo = getsysindexName(catNameTwo,it_sys->first,sys,-1);
		     std::string sysUNameTwo = getsysindexName(catNameTwo,it_sys->first,sys,1);
		     std::map<std::string,TH1F>::iterator itOneD = m_th1f_.find(sysDNameOne);
		     std::map<std::string,TH1F>::iterator itOneU = m_th1f_.find(sysUNameOne);
		     std::map<std::string,TH1F>::iterator itTwoD = m_th1f_.find(sysDNameTwo);
		     std::map<std::string,TH1F>::iterator itTwoU = m_th1f_.find(sysUNameTwo);
		     if (itOneD!=m_th1f_.end()) {
			mergeHistograms(sysDNameOne,&(itOneD->second),&(itTwoD->second));
			mergeHistograms(sysUNameOne,&(itOneU->second),&(itTwoU->second));
		     } else {
			  std::cerr << "WARNING -- RooContainer::MergeHistograms -- No (systematic) Histograms found named "
			  	    << sysDNameOne
			  	    << std::endl;

		     }	
		  }
	        }
	  }

	}
	
}
// -----------------------------------------------------------------------------------------
void RooContainer::mergeHistograms(std::string nameHist, TH1F* hist1, TH1F* hist2){
   
   //Get Bin Low edges of histogram 1
   int nbins1 = hist1->GetNbinsX();
   int nbins2 = hist2->GetNbinsX();
   int nbinsTot = nbins1+nbins2;

   double *arrBins1 = new double[nbins1];
   double *arrBins2 = new double[nbins2+1];
   double *arrBinsTot = new double[nbinsTot+1];

   for (int i=1;i<=nbins1;i++){
     arrBinsTot[i-1]=hist1->GetBinLowEdge(i);
   }
   for (int i=1;i<=nbins2+1;i++){	// Include upper edge in 2nd hist
     arrBinsTot[i+nbins1-1]=hist2->GetBinLowEdge(i);
   }

   const char *histoname = hist1->GetName();
   const char *histotitle = hist1->GetTitle();

   TH1F *newHist = new TH1F(Form("NUMPTYNAME%s",hist1->GetName()),histotitle,nbinsTot,arrBinsTot);
   newHist->SetName(hist1->GetName());
   for (int i=1;i<=nbins1;i++){
	newHist->SetBinContent(i,hist1->GetBinContent(i));
	newHist->SetBinError(i,hist1->GetBinError(i));
   } 
   for (int i=1;i<=nbins2;i++){
	newHist->SetBinContent(i+nbins1,hist2->GetBinContent(i));
	newHist->SetBinError(i+nbins1,hist2->GetBinError(i));
   } 

   std::cout << "RooContainer::MergeHistograms -- Replacing histogram - " 
	     << newHist->GetName()
	     << std::endl;
   // Now the dangerous part!
   *hist1 = *newHist;
   
   
}
// -----------------------------------------------------------------------------------------

void RooContainer::SumMultiBinnedDatasets(std::string new_name, std::vector<std::string > data, std::vector<double> normalisation,bool scale){
	// this version takes a vector of histograms and sums them to the total normalisations. 
	// if scale == true, then each histogram is summed together and then the total scaled to N
	// if scale == false then each is scaled to N/numberofhists then added
	// if N < 0 then total histogram is not scaled (if scale=false then N=1)
   if (normalisation.size() != ncat ){
	std::cerr << "WARNING!! -- RooContainer::SumMultiBinnedDataSets -- number of coefficients should be the same as number of categories " << std::endl;
   } else {

	for (int cat=0;cat<ncat;cat++){
	   std::vector<std::string> catNames;
	   for (std::vector<std::string>::iterator it = data.begin() ; it!=data.end();it++){
		catNames.push_back(getcatName(*it,cat));
	   }
	   sumMultiBinnedDatasets(getcatName(new_name,cat),catNames,normalisation[cat],scale);
	}
	
   }
}

void RooContainer::SumMultiBinnedDatasets(std::string new_name, std::vector<std::string > data, double normalisation,bool scale){
	// this version takes a vector of histograms and sums them to the total normalisations. 
	// if scale == true, then each histogram is summed together and then the total scaled to N
	// if scale == false then each is scaled to N/numberofhists then added
	// if N < 0 then total histogram is not scaled (if scale=false then N=1)

	for (int cat=0;cat<ncat;cat++){
	   std::vector<std::string> catNames;
	   for (std::vector<std::string>::iterator it = data.begin() ; it!=data.end();it++){
		catNames.push_back(getcatName(*it,cat));
	   }
	   sumMultiBinnedDatasets(getcatName(new_name,cat),catNames,normalisation,scale);
	}
	
}

void RooContainer::sumMultiBinnedDatasets(std::string new_name, std::vector<std::string > data, double normalisation,bool scale){
	
	int nHists = data.size();
   	std::map<std::string,TH1F>::iterator it_one = m_th1f_.find(data[0]); // Get The first histogram

	if (!scale && normalisation<0) normalisation=1.0;
	
   	if (it_one !=m_th1f_.end() ){

       	  TH1F *histOne = (TH1F*)((*it_one).second).Clone();
          histOne->SetName(Form("th1f_%s",new_name.c_str()));

	  if (!scale) histOne->Scale(1.0/histOne->Integral());

	  for (int i=1;i<nHists;i++){
	
   	   std::map<std::string,TH1F>::iterator it_two = m_th1f_.find(data[i]); // Get The first histogram
      	   if (scale) {  // histograms are weighted by their own integrals
		histOne->Add(&(it_two->second));
	   } else {	// histograms equally weighted
		histOne->Add(&(it_two->second),1.0/(it_two->second).Integral());
		
	   }
	  } 
	  if (normalisation > 0)  histOne->Scale(normalisation/histOne->Integral());	
   	  m_th1f_.insert(std::pair<std::string,TH1F>(new_name,*histOne));
	  std::cout << "RooContainer::SumBinnedDatasets -- Created New Histogram called " << new_name << std::endl;
  
  	} else {
	   std::cerr << "WARNING -- RooContainer::SumBinnedDatasets -- The following Histogram wasn't found " 
	             << data[0] << std::endl;
   	}
}
	

// -----------------------------------------------------------------------------------------
void RooContainer::sumBinnedDatasets(std::string new_name,std::string data_one,std::string data_two,double c1, double c2, bool scale){

   std::map<std::string,TH1F>::iterator it_one = m_th1f_.find(data_one);
   std::map<std::string,TH1F>::iterator it_two = m_th1f_.find(data_two);

   if (it_one !=m_th1f_.end() && it_two !=m_th1f_.end() ){
   
      TH1F *histOne = (TH1F*)((*it_one).second).Clone();
      histOne->SetName(Form("th1f_%s",new_name.c_str()));
      if (scale) {  // coefficients are just multiples each histogram
	histOne->Scale(c1);
	histOne->Add(&(it_two->second),c2);
      } else {
	histOne->Scale((c1/histOne->Integral()));
	histOne->Add(&(it_two->second),(c2/((*it_two).second).Integral()));
      }

   m_th1f_.insert(std::pair<std::string,TH1F>(new_name,*histOne));
   std::cout << "RooContainer::SumBinnedDatasets -- Created New Histogram called " << new_name << std::endl;
  
   } else {
	std::cerr << "WARNING -- RooContainer::SumBinnedDatasets -- One of the following Histograms wasn't found " 
	          << data_one << ", " << data_two << std::endl;
   }
   
}
// -----------------------------------------------------------------------------------------
void RooContainer::CombineBinnedDatasets(std::string data_one, std::string data_two, double fraction){
  for (int cat=0;cat<ncat;cat++) {
    combineBinnedDatasets(getcatName(data_one,cat),getcatName(data_two,cat),fraction);
  }    
}
// -----------------------------------------------------------------------------------------
void RooContainer::combineBinnedDatasets(std::string data_one, std::string data_two, double fraction){

   std::map<std::string,TH1F>::iterator it_one = m_th1f_.find(data_one);
   std::map<std::string,TH1F>::iterator it_two = m_th1f_.find(data_two);

     if (it_one !=m_th1f_.end() && it_two !=m_th1f_.end() ){

        if (fraction > -1){
          // currently the histograms are in the ratio 1:fraction_mc, must have  instead weighted to 1:fraction_data
	  // assume fraction is fraction_data/fraction_mc
	  double N1 = (it_one->second).Integral();
	  double N2 = (it_two->second).Integral();
	  double f_mc = N2/(N1+N2);
      double total = N1 + N2;

 	  double f_data = fraction*f_mc;

          // scale so that the two are in f_data proportion.
	  (it_one->second).Scale(total*(1-f_data)/N1);
	  (it_two->second).Scale(total*f_data/N2);
	  (it_one->second).Add(&(it_two->second));

	  std::cout << "RooContainer::CombinedBinnedDatasets -- Added second histogram to first with proportion corrected by " << fraction << std::endl;
	  // scale second one back to its original 

	  (it_two->second).Scale(N2/(it_two->second).Integral());
	  std::cout << "  seccond histogram is left unchanged " << std::endl;
	
	  
        } else {

	  (it_one->second).Add(&(it_two->second));
        }

     } else {

	std::cerr << "WARNING -- RooContainer::CombineBinnedDatasets -- Cannot find one of "
		  << data_one << " " << data_two << " "
		  << " Datasets will not be combined -- WARNING!" << std::endl;
	return;
     }
}
// -----------------------------------------------------------------------------------------
void RooContainer::WriteSpecificCategoryDataCards(std::string filename,std::string data_name
						 ,std::string sig_name, std::string bkg_name){

   for (int cat=0;cat<ncat;cat++){
	writeSpecificCategoryDataCard(cat,filename,data_name,sig_name,bkg_name);
   }
}
// -----------------------------------------------------------------------------------------
void RooContainer::writeSpecificCategoryDataCard(int cat,std::string filename,std::string data_name
						  ,std::string sig_name, std::string bkg_name){

   bool parameterisedBackground = false;
   bool multi_pdf = false;
   RooAbsPdf * pdf_ptr;
   RooRealVar *obs;

     std::map<std::string,TH1F>::iterator it_data = m_th1f_.find(getcatName(data_name,cat));

     if (!(it_data !=m_th1f_.end() )){

	std::cerr << "WARNING -- RooContainer::WriteDataCard -- Cannot find "
		  << data_name 
		  << " DataCard will NOT be Writen -- WARNING!" << std::endl;
	return;
     }

     std::map<std::string,TH1F>::iterator it_bkg  = m_th1f_.find(getcatName(bkg_name,cat));
     if (it_bkg==m_th1f_.end() ){
	parameterisedBackground = true;	
        std::map<std::string,RooAddPdf>::iterator pdf_bkg_a  = m_pdf_.find(getcatName(bkg_name,cat));
        std::map<std::string,RooExtendPdf>::iterator pdf_bkg_e  = m_exp_.find(getcatName(bkg_name,cat));
        if (pdf_bkg_e==m_exp_.end() && pdf_bkg_a==m_pdf_.end()){ 
	   std::cerr << "WARNING -- RooContainer::WriteDataCard -- Cannot find "
		     << bkg_name 
		     << " DataCard will NOT be Writen -- WARNING!" << std::endl;
	   return;
        } else if(pdf_bkg_e!=m_exp_.end() ){
	  multi_pdf = false;
   	  pdf_ptr = m_gen_[getcatName(bkg_name,cat)];
   	  obs = m_data_var_ptr_[getcatName(data_name,cat)] ;
	} else if(pdf_bkg_a!=m_pdf_.end() ){
	  multi_pdf = true;
   	  pdf_ptr = m_gen_[getcatName(bkg_name,cat)];
   	  obs = m_data_var_ptr_[getcatName(data_name,cat)] ;
	}
     }
   std::string signal_mass_name = (std::string) Form("%s_m$MASS",sig_name.c_str());
   ofstream file;
   if (parameterisedBackground){
     file.open(Form("%s_cms-hgg-datacard_parBKG_cat%d.txt",filename.c_str(),cat));
     file << "CMS-HGG DataCard for Binned Limit Setting with RooDataHist+Parameterised Background\n";
   }
   else{
     file.open(Form("%s_cms-hgg-datacard_cat%d.txt",filename.c_str(),cat));
     file << "CMS-HGG DataCard for Binned Limit Setting with RooDataHist\n";
   }

   file << "Run with: combine cms-hgg-datacard_cat"<< cat<<".txt -M Routine -D "<< data_name << " -m MASS --generateBinnedWorkaround -S 1"; 
   file << "\n---------------------------------------------\n";
   file << "imax *\n";
   file << "jmax *\n";
   file << "kmax *\n";
   file << "---------------------------------------------\n";
   if (parameterisedBackground){
     file << "shapes "<<data_name <<" * "<<filename<<" cms_hgg_workspace:roohist_"<< data_name<<"_$CHANNEL\n";
     file << "shapes "<<signal_mass_name<<" * "<<filename<<" cms_hgg_workspace:roohist_"<<signal_mass_name << "_$CHANNEL cms_hgg_workspace:roohist_"<< signal_mass_name <<"_$CHANNEL_$SYSTEMATIC01_sigma\n";
     file << "shapes "<<bkg_name<<" * "<<filename<<" cms_hgg_workspace:pdf_"<<bkg_name << "_$CHANNEL cms_hgg_workspace:roohist_"<< bkg_name <<"_$CHANNEL_$SYSTEMATIC01_sigma\n";
   }
   else 
     file << "shapes * * "<<filename<<" cms_hgg_workspace:roohist_$PROCESS_$CHANNEL cms_hgg_workspace:roohist_$PROCESS_$CHANNEL_$SYSTEMATIC01_sigma\n";
   file << "---------------------------------------------\n";

   // Write the data info	 
   file << "bin	           ";
    file << "cat"<<cat << " ";
   file << "\nobservation  ";
    file << " -1 ";
   file << "\n---------------------------------------------\n";
   // Write the model info	 
   file << "bin	       ";
    file << "cat"<<cat << "  cat" <<cat<<" ";
   file << "\nprocess  ";
    file << signal_mass_name<<" "<<bkg_name<< " ";
   file << "\nprocess  ";
    file << " -1   1  ";
   file << "\nrate     ";
   if (! parameterisedBackground)
   	 file << " -1  -1 ";
   else {
	 {
       //   std::string bcatName = getcatName(bkg_name,cat);
       //   std::string dcatName = getcatName(data_name,cat);
//	  std::string rngeName = getcatName("datacard_",cat);
//	  double r1 = m_var_min_[dcatName] ;
//	  double r2 = m_var_max_[dcatName] ;
//	  double norm = getNormalisationFromFit(bcatName,rngeName,m_gen_[bcatName],obs,r1,r2,multi_pdf);
	  file << " -1   1 " << " ";   
	}
   }
   file << "\n---------------------------------------------\n";

   // Now write the systematics lines:
   // first Global Systematics
   for (std::map<std::string,std::pair<double,double> >::iterator it_g_sys = global_systematics_.begin();it_g_sys!=global_systematics_.end(); it_g_sys++){ 
     file << it_g_sys->first << " lnN  ";
      file << (it_g_sys->second).first << " " << (it_g_sys->second).second << " ";	
     file <<std::endl;
   }

  
   // Write Category Normalisations 
   for (std::map<std::string,std::vector < std::pair< double,double> > >::iterator it_n_sys = normalisation_systematics_.begin();it_n_sys!=normalisation_systematics_.end(); it_n_sys++){
     file << it_n_sys->first << " lnN ";
     file << (((*it_n_sys).second)[cat].first) << " " << (((*it_n_sys).second)[cat].second) << " ";	
     file <<std::endl;

   }
    
   double sigmaUnitInv = 1./(sigmaRange/nsigmas);
   for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
    
     file << it_sys->first << " shape  ";
     if (it_sys->second == 0) // both background and signal effected
       file << sigmaUnitInv<< " " << sigmaUnitInv<< " ";
     else if (it_sys->second == -1) // signal effected
       file << sigmaUnitInv <<" 0 ";
     else if (it_sys->second == 1) // background effected
       file << " 0 "<<sigmaUnitInv;
     file << std::endl;
   }


   file.close();
}
// -----------------------------------------------------------------------------------------
void RooContainer::WriteDataCard(std::string filename,std::string data_name
				,std::string sig_name, std::string bkg_name){

   bool parameterisedBackground = false;
   bool multi_pdf = false;
   RooAbsPdf * pdf_ptr;
   RooRealVar *obs;

   for (int cat=0;cat<ncat;cat++){
     std::map<std::string,TH1F>::iterator it_data = m_th1f_.find(getcatName(data_name,cat));

     if (!(it_data !=m_th1f_.end() )){

	std::cerr << "WARNING -- RooContainer::WriteDataCard -- Cannot find "
		  << data_name 
		  << " DataCard will NOT be Writen -- WARNING!" << std::endl;
	return;
     }
   }

   for (int cat=0;cat<ncat;cat++){
     std::map<std::string,TH1F>::iterator it_bkg  = m_th1f_.find(getcatName(bkg_name,cat));
     if (it_bkg==m_th1f_.end() ){
	parameterisedBackground = true;	
        std::map<std::string,RooAddPdf>::iterator pdf_bkg_a  = m_pdf_.find(getcatName(bkg_name,cat));
        std::map<std::string,RooExtendPdf>::iterator pdf_bkg_e  = m_exp_.find(getcatName(bkg_name,cat));
        if (pdf_bkg_e==m_exp_.end() && pdf_bkg_a==m_pdf_.end()){ 
	   std::cerr << "WARNING -- RooContainer::WriteDataCard -- Cannot find "
		     << bkg_name 
		     << " DataCard will NOT be Writen -- WARNING!" << std::endl;
	   return;
        } else if(pdf_bkg_e!=m_exp_.end() ){
	  multi_pdf = false;
   	  pdf_ptr = m_gen_[getcatName(bkg_name,cat)];
   	  obs = m_data_var_ptr_[getcatName(data_name,cat)] ;
	} else if(pdf_bkg_a!=m_pdf_.end() ){
	  multi_pdf = true;
   	  pdf_ptr = m_gen_[getcatName(bkg_name,cat)];
   	  obs = m_data_var_ptr_[getcatName(data_name,cat)] ;
	}
     }
   }
   std::string signal_mass_name = (std::string) Form("%s_m$MASS",sig_name.c_str());
   ofstream file;
   if (parameterisedBackground){
     file.open(Form("%s_cms-hgg-datacard_parBKG.txt",filename.c_str()));
     file << "CMS-HGG DataCard for Binned Limit Setting with RooDataHist+Parameterised Background\n";
   }
   else{
     file.open(Form("%s_cms-hgg-datacard.txt",filename.c_str()));
     file << "CMS-HGG DataCard for Binned Limit Setting with RooDataHist\n";
   }

   file << "Run with: combine cms-hgg-datacard.txt -M Routine -D "<< data_name << " -m MASS --generateBinnedWorkaround -S 1"; 
   file << "\n---------------------------------------------\n";
   file << "imax *\n";
   file << "jmax *\n";
   file << "kmax *\n";
   file << "---------------------------------------------\n";
   if (parameterisedBackground){
     file << "shapes "<<data_name <<" * "<<filename<<" cms_hgg_workspace:roohist_"<< data_name<<"_$CHANNEL\n";
     file << "shapes "<<signal_mass_name<<" * "<<filename<<" cms_hgg_workspace:roohist_"<<signal_mass_name << "_$CHANNEL cms_hgg_workspace:roohist_"<< signal_mass_name <<"_$CHANNEL_$SYSTEMATIC01_sigma\n";
     file << "shapes "<<bkg_name<<" * "<<filename<<" cms_hgg_workspace:pdf_"<<bkg_name << "_$CHANNEL cms_hgg_workspace:roohist_"<< bkg_name <<"_$CHANNEL_$SYSTEMATIC01_sigma\n";
   }
   else 
     file << "shapes * * "<<filename<<" cms_hgg_workspace:roohist_$PROCESS_$CHANNEL cms_hgg_workspace:roohist_$PROCESS_$CHANNEL_$SYSTEMATIC01_sigma\n";
   file << "---------------------------------------------\n";

   // Write the data info	 
   file << "bin	           ";
   for (int cat=0;cat<ncat;cat++) file << "cat"<<cat << " ";
   file << "\nobservation  ";
   for (int cat=0;cat<ncat;cat++) file << " -1 ";
   file << "\n---------------------------------------------\n";
   // Write the model info	 
   file << "bin	       ";
   for (int cat=0;cat<ncat;cat++) file << "cat"<<cat << "  cat" <<cat<<" ";
   file << "\nprocess  ";
   for (int cat=0;cat<ncat;cat++) file << signal_mass_name<<" "<<bkg_name<< " ";
   file << "\nprocess  ";
   for (int cat=0;cat<ncat;cat++) file << " -1   1  ";
   file << "\nrate     ";
   if (! parameterisedBackground)
   	for (int cat=0;cat<ncat;cat++) file << " -1  -1 ";
   else {
	for (int cat=0;cat<ncat;cat++) {
       //   std::string bcatName = getcatName(bkg_name,cat);
       //   std::string dcatName = getcatName(data_name,cat);
//	  std::string rngeName = getcatName("datacard_",cat);
//	  double r1 = m_var_min_[dcatName] ;
//	  double r2 = m_var_max_[dcatName] ;
//	  double norm = getNormalisationFromFit(bcatName,rngeName,m_gen_[bcatName],obs,r1,r2,multi_pdf);
	  file << " -1   1 " << " ";   
	}
   }
   file << "\n---------------------------------------------\n";

   // Now write the systematics lines:
   // first Global Systematics
   for (std::map<std::string,std::pair<double,double> >::iterator it_g_sys = global_systematics_.begin();it_g_sys!=global_systematics_.end(); it_g_sys++){ 
     file << it_g_sys->first << " lnN  ";
     for (int cat=0;cat<ncat;cat++) file << (it_g_sys->second).first << " " << (it_g_sys->second).second << " ";	
     file <<std::endl;
   }

   // Write Category Normalisations 
   for (std::map<std::string,std::vector < std::pair< double,double> > >::iterator it_n_sys = normalisation_systematics_.begin();it_n_sys!=normalisation_systematics_.end(); it_n_sys++){
     file << it_n_sys->first << " lnN ";
     for (int cat=0;cat<ncat;cat++) file << ((*it_n_sys).second[cat].first) << " " << ((*it_n_sys).second[cat].second) << " ";	
     file <<std::endl;

   }


  
   // Finally Shape Systematics    
   double sigmaUnitInv = 1./(sigmaRange/nsigmas);
   for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
    
     file << it_sys->first << " shape  ";
     if (it_sys->second == 0) // both background and signal effected
      for (int cat=0;cat<ncat;cat++) file << sigmaUnitInv<< " " << sigmaUnitInv<< " ";
     else if (it_sys->second == -1) // signal effected
      for (int cat=0;cat<ncat;cat++) file << sigmaUnitInv <<" 0 ";
     else if (it_sys->second == 1) // background effected
      for (int cat=0;cat<ncat;cat++) file << " 0 "<<sigmaUnitInv;
     file << std::endl;
   }


   file.close();
}
// ----------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
std::string RooContainer::getcatName(std::string name, int i){
  char* char_string = Form("%s_cat%d",name.c_str(),i);
  std::string output(char_string);
  return output;
}
// ---------------------------------------------------------------------------------
std::string RooContainer::getweightName(std::string name){
  char* char_string = Form("weight_%s",name.c_str());
  std::string output(char_string);
  return output;
}
// ---------------------------------------------------------------------------------
std::string RooContainer::getnormName(std::string name){
  char* char_string = Form("norm_%s",name.c_str());
  std::string output(char_string);
  return output;
}
// ----------------------------------------------------------------------------------------------------
std::string RooContainer::getsysindexName(std::string name,std::string sys_name
				    ,int sys,int direction){
  char* char_string;
  if (direction >= 0)
    char_string = Form("%s_%sUp%.2d_sigma",name.c_str(),sys_name.c_str(),sys);
  else
    char_string = Form("%s_%sDown%.2d_sigma",name.c_str(),sys_name.c_str(),sys);
  std::string output(char_string);
  return output;
}

// ----------------------------------------------------------------------------------------------------
std::string RooContainer::getsysName(std::string name,std::string sys_name){
  char* char_string = Form("%s_%s",name.c_str(),sys_name.c_str());
  std::string output(char_string);
  return output;
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,double xmin,double xmax){

  RooRealVar temp(name.c_str(),name.c_str(),xmin,xmax);
  m_real_var_.insert(std::pair<std::string,RooRealVar >(name,temp));

  if (verbosity_){
  	std::cout << "RooContainer::AddRealVar -- Appended the variable " 
		  << name <<std::endl; 
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,double init,double xmin,double xmax){
  RooRealVar temp(name.c_str(),name.c_str(),init,xmin,xmax);
  m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp));
  if (verbosity_){
            std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	              << name <<std::endl;
  }
  
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::addFormulaVar(std::string name ,std::string formula, std::string var){

  RooFormulaVar temp(name.c_str(),name.c_str(),formula.c_str(),RooArgList(m_real_var_[var]));
  m_form_var_.insert(std::pair<std::string,RooFormulaVar >(name,temp));

  if (verbosity_){
        std::cout << "RooContainer::AddFormulaVar -- Appended the variable " 
	          << name <<std::endl; 
  }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::addGenericPdf(std::string name,std::string formula,std::string obs_name
				,std::vector<std::string> & var, int form
				,double norm_guess, double norm_min, double norm_max ){

   

    RooAbsPdf *temp_1;
    std::map<std::string,RooRealVar>::iterator obs_real_var = m_real_var_.find(obs_name);

    if (obs_real_var == m_real_var_.end()) {

         std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	           << obs_name << std::endl; 
    }
    else{

     // find the map to look in 
     
     bool form_vars = false;
     std::vector<std::string>::iterator it_var = var.begin();

     for (;it_var != var.end();it_var++){ 
        std::map<std::string,RooRealVar>::iterator it_check = m_real_var_.find(*it_var);
	if (it_check == m_real_var_.end()){
	  
           std::map<std::string,RooFormulaVar>::iterator it_f_check = m_form_var_.find(*it_var);
	   if (it_f_check == m_form_var_.end()){
             std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	               << *it_var << std::endl; 
	     return;
	   } else form_vars = true;
        }
     }

     if (form==0){
      RooArgList roo_args;
      roo_args.add((*obs_real_var).second); 
      for (std::vector<std::string>::iterator it_var = var.begin()
	   ;it_var != var.end()
	   ;it_var++
	   ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     if (!form_vars) {
		std::map<std::string,RooRealVar>::iterator real_var =  m_real_var_.find(*it_var);
	        if (real_var != m_real_var_.end())
	          roo_args.add((*real_var).second);
	        else 
      		   std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	     << *it_var << std::endl;
	    } else {
		std::map<std::string,RooFormulaVar>::iterator form_var =  m_form_var_.find(*it_var);
	        if (form_var != m_form_var_.end())
	          roo_args.add((*form_var).second);
	        else 
      		   std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	     << *it_var << std::endl;
  	   }
	}

      std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	        << std::endl;

      temp_1 = new RooGenericPdf(Form("pdf_%s",name.c_str()),name.c_str(),formula.c_str(),roo_args);	
      //v_gen_.push_back(temp_1);

     } else if (form == 1) { //RooExponential  - x,slope
	if (var.size() == 1){
          if (!form_vars)  temp_1 = new RooExponential(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]]);
          else  temp_1 = new RooExponential(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_form_var_[var[0]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 1 arguments for RooExponential, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 2) { //RooGaussian - x,mean,sigma
	if (var.size() == 2){
	  if (!form_vars) temp_1 = new RooGaussian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]);
	  else temp_1 = new RooGaussian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_form_var_[var[0]],m_form_var_[var[1]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 2 arguments for RooGaussian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 3) { //RooBreitWigner - x,centre,width
	if (var.size() == 2){
	  if (!form_vars) temp_1 = new RooBreitWigner(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]);
	  else temp_1 = new RooBreitWigner(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_form_var_[var[0]],m_form_var_[var[1]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 2 arguments for RooBreitWigner, was given: "

		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 4) { //RooCBShape - x,centre,sigma,slope,n
	if (var.size() == 4){
	  if (!form_vars) temp_1 = new RooCBShape(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]
				 ,m_real_var_[var[2]],m_real_var_[var[3]]);
	  else temp_1 = new RooCBShape(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_form_var_[var[0]],m_form_var_[var[1]]
				 ,m_form_var_[var[2]],m_form_var_[var[3]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 4 arguments for RooCBShape, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 5) { //RooVoigtian - x,centre,width,sigma
	if (var.size() == 3){
	  if (!form_vars) temp_1 = new RooVoigtian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]
				  ,m_real_var_[var[2]]);
	  else temp_1 = new RooVoigtian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_form_var_[var[0]],m_form_var_[var[1]]
				  ,m_form_var_[var[2]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 3 arguments for RooVoigtian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 6) { //RooPowerLaw - x,slope
	if (var.size() == 1){
          if (!form_vars)  temp_1 = new RooGenericPdf(Form("pdf_%s",name.c_str()),name.c_str(),"pow(@0,-@1)", 
						      RooArgList( (*obs_real_var).second,m_real_var_[var[0]] ) );
          else  temp_1 = new RooGenericPdf(Form("pdf_%s",name.c_str()),name.c_str(),"pow(@0,-@1)", 
					   RooArgList( (*obs_real_var).second, m_form_var_[var[0]]) );
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 1 arguments for RooPowerLaw, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form>60 && form < 70) { // RooChebychev - x, p1....pform-1

	if ((int)var.size() == form-60){
      	  RooArgList roo_args;
	  if (!form_vars){
	   for (std::vector<std::string>::iterator it_var = var.begin()
	      ;it_var != var.end()
	      ;it_var++
	      ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     std::map<std::string,RooRealVar>::iterator real_var =  m_real_var_.find(*it_var);
	     if (real_var != m_real_var_.end())
	       roo_args.add((*real_var).second);
	     else 
      		std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	  << *it_var << std::endl;
  	   }
	  } else {
	   for (std::vector<std::string>::iterator it_var = var.begin()
	      ;it_var != var.end()
	      ;it_var++
	      ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     std::map<std::string,RooFormulaVar>::iterator form_var =  m_form_var_.find(*it_var);
	     if (form_var != m_form_var_.end())
	       roo_args.add((*form_var).second);
	     else 
      		std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	         	  << *it_var << std::endl;
	   }
	  }
      	   std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	        	  << std::endl;

      	   temp_1 = new RooChebychev(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need "<<form-60 << " arguments for RooChebychev of order " << form-60 <<", was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}

     } else if (form>70) { // RooBernstein (FAST) - x, p1....pform-1
	int norders = form-70;
	if ((int)var.size() == form-70){
      	  RooArgList roo_args;
	  //roo_args.add(RooConst(1.0)); // no need with RooBernsteinFast!
	  if (!form_vars){
	   for (std::vector<std::string>::iterator it_var = var.begin()
	      ;it_var != var.end()
	      ;it_var++
	      ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     std::map<std::string,RooRealVar>::iterator real_var =  m_real_var_.find(*it_var);
	     if (real_var != m_real_var_.end())
	       roo_args.add((*real_var).second);
	     else 
      		std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	  << *it_var << std::endl;
  	   }
	  } else {
	   for (std::vector<std::string>::iterator it_var = var.begin()
	      ;it_var != var.end()
	      ;it_var++
	      ){
	     std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		       << *it_var << std::endl;

	     std::map<std::string,RooFormulaVar>::iterator form_var =  m_form_var_.find(*it_var);
	     if (form_var != m_form_var_.end())
	       roo_args.add((*form_var).second);
	     else 
      		std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	         	  << *it_var << std::endl;
	   }
	  }
      	   std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	        	  << std::endl;

      	   if 	   (norders==1)	temp_1 = new RooBernsteinFast<1>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
      	   else if (norders==2) temp_1 = new RooBernsteinFast<2>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
      	   else if (norders==3) temp_1 = new RooBernsteinFast<3>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
      	   else if (norders==4) temp_1 = new RooBernsteinFast<4>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
           else if (norders==5)	temp_1 = new RooBernsteinFast<5>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
      	   else if (norders==6) temp_1 = new RooBernsteinFast<6>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
      	   else if (norders==7) temp_1 = new RooBernsteinFast<7>(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
	   else temp_1 = new RooBernstein(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args); // switch back to Usual RooBernstein
	   
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need "<<form-70 << " arguments for RooBernstein of order " << form-70 <<", was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}

    	
    } else {
	
       std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Mode " << form
		 << "Understood -- WARNING"
	         << std::endl;
       return;
     }

     m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,temp_1));

     RooRealVar temp_var(Form("pdf_%s_norm",name.c_str()),name.c_str(),norm_guess,norm_min,norm_max);
     m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));

     RooExtendPdf temp(name.c_str(),name.c_str(),*temp_1,m_real_var_[name],"FullObervableRange");
     m_exp_.insert(std::pair<std::string,RooExtendPdf>(name,temp));

     std::cout << "RooContainer::AddGenericPdf -- Made extended PDF " 
	       << name << std::endl;	
    }
}

// ----------------------------------------------------------------------------------------------------
/*
void RooContainer::convolutePdf(std::string name, std::string f_pdf, std::string g_pdf, RooRealVar &observable, double norm_guess){
 
    observable.setBins(1000,"cache");
    
    std::map<std::string,RooAbsPdf*>::iterator it_f = m_gen_.find(f_pdf);
    std::map<std::string,RooAbsPdf*>::iterator it_g = m_gen_.find(g_pdf);
    if (it_f == m_gen_.end() || it_f == m_gen_.end()){
	
      std::cerr << "WARNING -- RooContainer::ConvolutePdf -- One of pdf's f or g wasnt found"
		<< " f: " << f_pdf
		<< " g: " << g_pdf << std::endl;
      return;
    }

    RooFFTConvPdf *convoluted = new RooFFTConvPdf(Form("pdf_%s",name.c_str()),name.c_str(),observable,*(it_f->second),*(it_g->second));

    RooRealVar temp_var(Form("norm_%s",name.c_str()),name.c_str(),norm_guess,0.,1.e6);

    m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));

    RooExtendPdf temp(name.c_str(),name.c_str(),*convoluted,m_real_var_[name]);
    m_exp_.insert(std::pair<std::string,RooExtendPdf>(name,temp));
   
    std::cout << "RooContainer::ConvolvePdf -- Made convoluted PDF " 
	      << name << std::endl;	
    
}
*/
// ----------------------------------------------------------------------------------------------------
std::vector<std::vector<double> > RooContainer::SoverBOptimizedBinning(std::string signalname,std::string bkgname,int nTargetBins,double penaltyScale){

	std::vector<std::vector<double> > return_bins;
	for (int cat=0;cat<ncat;cat++){
	   std::map<std::string,TH1F>::iterator it_ths=m_th1f_.find(getcatName(signalname,cat));
	   std::map<std::string,TH1F>::iterator it_thb=m_th1f_.find(getcatName(bkgname,cat));
	   if (it_ths!=m_th1f_.end() && it_thb!=m_th1f_.end()){
		return_bins.push_back(soverBOptimizedBinning(&(it_ths->second),&(it_thb->second),nTargetBins,penaltyScale));

	   } else {
		std::cerr << "WARNING ! -- RooContainer::SoverBOptimizedBinning -- One of the Following binned datasets not found " << signalname << ", " << bkgname <<std::endl;
	   }
	}

	return return_bins;

}

// ----------------------------------------------------------------------------------------------------
double RooContainer::calculateSig(double s1, double s2, double b1, double b2){
	
	
	double sig =  1.4142*TMath::Sqrt((s1+b1)*TMath::Log((s1+b1)/b1) + (s2+b2)*TMath::Log((s2+b2)/b2) - (s1+s2));
	return sig;
}

double RooContainer::calculateSigMulti(std::vector<double> &s1, std::vector<double> &b1){
	
	int nchannel=s1.size();
	double sterm=0;
	double logterms=0;
	for (int i=0;i<nchannel;i++){
	  logterms+=(s1[i]+b1[i])*TMath::Log((s1[i]+b1[i])/b1[i]);
	  sterm+=s1[i];
	}
	double sig =  1.4142*TMath::Sqrt(logterms - sterm);
	return sig;
}
double RooContainer::calculateSigMulti(double *s1, double *b1, int nchannel){
	
//	int nchannel=s1.size();
	double sterm=0;
	double logterms=0;
	for (int i=0;i<nchannel;i++){
	  logterms+=(s1[i]+b1[i])*TMath::Log((s1[i]+b1[i])/b1[i]);
	  sterm+=s1[i];
	}
	double sig =  1.4142*TMath::Sqrt(logterms - sterm);
	return sig;
}

bool RooContainer::compareLHWide(double s1, double sdiff,double s2,double b1,double bdiff, double b2,double n, std::vector<double> &chanS, std::vector<double> &chanB){


	std::vector<double> tmpS = chanS;
	std::vector<double> tmpB = chanB;
	std::vector<double> tmpS2 = chanS;
	std::vector<double> tmpB2 = chanB;
	
	double chan1 = calculateSig(s1+sdiff,s2,b1+bdiff,b2) ;
	double chan2 = calculateSig(s1,s2+sdiff,b1,b2+bdiff);

	if ((chan1 > chan2) ){
		return true;
	}
	else return false;


}

// ----------------------------------------------------------------------------------------------------
std::vector<std::vector<double> > RooContainer::SignificanceOptimizedBinning(std::string signalname,std::string bkgname,int nTargetBins){

	std::vector<std::vector<double> > return_bins;
	for (int cat=0;cat<ncat;cat++){
	   std::map<std::string,TH1F>::iterator it_ths=m_th1f_.find(getcatName(signalname,cat));
	   std::map<std::string,TH1F>::iterator it_thb=m_th1f_.find(getcatName(bkgname,cat));
	   if (it_ths!=m_th1f_.end() && it_thb!=m_th1f_.end()){
		return_bins.push_back(significanceOptimizedBinning(&(it_ths->second),&(it_thb->second),nTargetBins));

	   } else {
		std::cerr << "WARNING ! -- RooContainer::SignificanceOptimizedBinning -- One of the Following binned datasets not found " << signalname << ", " << bkgname <<std::endl;
	   }
	}

	return return_bins;

}

// Careful, Recursive Function here -------------------------------------------------------------------
void RooContainer::maxSigScan(double *maximumSignificance,int *frozen_counters,int *chosen_counters,TH1F *hs, TH1F *hb, int N,int *counters, int movingCounterIndex){

	if (movingCounterIndex >=N ) std::cout << "Should never have got here!!!" <<std::endl;
	if (movingCounterIndex < 0) return;
	if (counters[movingCounterIndex] < 2) std::cout << "WHAT IS GOING ON?? " <<  movingCounterIndex << " " << counters[movingCounterIndex]<<std::endl;
	int nBins = hs->GetNbinsX();
	double significance_now;
	if (not sweepmode){
	  // N is number of boundaries
	  int m=1;
	  for (int k=movingCounterIndex+1;k<N;k++) {
		counters[k]=counters[movingCounterIndex]+m;
		m++;
	  }
	} else {
	  // N is number of boundaries
	  int m=1;
	  for (int k=movingCounterIndex+1;k<N;k++) {
		int newpoint = (frozen_counters[k] - g_step > 1) ? frozen_counters[k] - g_step:2;
		counters[k]= (newpoint > counters[movingCounterIndex]+m )? newpoint: counters[movingCounterIndex]+m ;
		m++;
	  }
	}
	
	if (counters[movingCounterIndex] < 2) std::cout << "WHAT IS GOING ON?? " <<  movingCounterIndex << " " << counters[movingCounterIndex]<<std::endl;

	if ( movingCounterIndex==N-1) {	
	 if (not sweepmode){
	    for (;counters[N-1]<=nBins;counters[N-1]+=g_step){
		for (int j=0;j<=N-1;j++){
			if (j==0){
			  signalVector1[j] = (hs->Integral(1,counters[j]-1));
			  backgroundVector1[j] = (hb->Integral(1,counters[j]-1));

			} else {
			  signalVector1[j] = (hs->Integral(counters[j-1],counters[j]-1));
			  backgroundVector1[j] = (hb->Integral(counters[j-1],counters[j]-1));

			}
	   	 }
		signalVector1[N]=(hs->Integral(counters[N-1],nBins));
		backgroundVector1[N]=(hb->Integral(counters[N-1],nBins));
		significance_now = calculateSigMulti(signalVector1,backgroundVector1,N+1);

		if (significance_now>*maximumSignificance){
			*maximumSignificance=significance_now;
			for (int j=0;j<N;j++){
				chosen_counters[j]=counters[j];	
				if (chosen_counters[j] < 0) std::cout << "Freak OUT !!! - " << j << "  " << chosen_counters[j] <<std::endl;
			}
		}
	     }
	     maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex-1);
	 } else { // fine scanning
	    int currmax = (nBins < frozen_counters[N-1] + g_step )? nBins : frozen_counters[N-1] + g_step;
	    for (;counters[N-1]<=currmax;counters[N-1]++){
		for (int j=0;j<=N-1;j++){
			if (j==0){
			  signalVector1[j] = (hs->Integral(1,counters[j]-1));
			  backgroundVector1[j] = (hb->Integral(1,counters[j]-1));

			} else {
			  signalVector1[j] = (hs->Integral(counters[j-1],counters[j]-1));
			  backgroundVector1[j] = (hb->Integral(counters[j-1],counters[j]-1));

			}
	   	 }	
		signalVector1[N]=(hs->Integral(counters[N-1],nBins));
		backgroundVector1[N]=(hb->Integral(counters[N-1],nBins));
		significance_now = calculateSigMulti(signalVector1,backgroundVector1,N+1);
		if (significance_now>*maximumSignificance){
			*maximumSignificance=significance_now;
			for (int j=0;j<N;j++){
				chosen_counters[j]=counters[j];
			}
		}
	     }
	     maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex-1);
	 }
	}


	else if (counters[movingCounterIndex]+1 <= nBins-(N-movingCounterIndex)){
		
	 if (not sweepmode){
		if (counters[movingCounterIndex]+g_step <= nBins-(N-movingCounterIndex)){
			counters[movingCounterIndex]+=g_step;
		} else {
			counters[movingCounterIndex]=nBins-(N-movingCounterIndex);
		}
		int m=1;
		for (int k=movingCounterIndex+1;k<N;k++) {
		  counters[k]=counters[movingCounterIndex]+m;
		  m++;
	        }
		maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex+1);
	 } else {
	  	// N is number of boundaries
		if (counters[movingCounterIndex]+1 <= frozen_counters[movingCounterIndex]+g_step){
		  counters[movingCounterIndex]++;
	  	  int m=1;
	  	  for (int k=movingCounterIndex+1;k<N;k++) {
			int newpoint = (frozen_counters[k] - g_step > 1) ? frozen_counters[k] - g_step:2;
			if (newpoint < 0) {std::cout << "Whaaaaaa? "<< newpoint <<std::endl;}
			counters[k]= (newpoint > counters[movingCounterIndex]+m )? newpoint: counters[movingCounterIndex]+m ;
			m++;
	  	  }
		  maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex+1);
		} else {
			if (movingCounterIndex>0){
				maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex-1);
			} else {
				return;
			}
		}
	 }

	}
	
	else { // got to the end,

		if (movingCounterIndex>0){
			maxSigScan(maximumSignificance,frozen_counters,chosen_counters,hs,hb,N,counters,movingCounterIndex-1);
		} else {
			return;
		}
	}


}

// ----------------------------------------------------------------------------------------------------
std::vector<double> RooContainer::significanceOptimizedBinning(TH1F *hs,TH1F *hb,int nTargetBins){

	// Performs Optimized Binning based on a Signal and Background  distributions
	// First runs the optimizedBinning on background and rebins S and B clones, note, always performs 
	// revise_target=false,direction=-1 and use_n_entries=true
	// nTargetBins is used for the flat binning, decision to merge is based on improvement to expected significance
	// Full scan is done for largest significance (wardning, could be very slow for tight constraints)

	int ninitBins = hb->GetNbinsX();
	if (hs->Integral()==0 ||  hb->Integral()==0 || ninitBins < 3) {
		std::vector<double> binEdges;
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(ninitBins+1));
		return binEdges;
	}

	std::vector<double> binEdges = optimizedReverseBinning(hb,nTargetBins,false,true);

	if (binEdges.size() <= 10 ) return binEdges;
	// Just TESTING HERE so remove this line soon!
//	nTargetBins = 150; // this gives us about 144 with the latest thing :)
//	std::vector<double> binEdges = optimizedReverseBinning(hb,nTargetBins,false,false);
	//---------------------------------------------------------------------------------

	int j =0;
	double *arrBins = new double[binEdges.size()];
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		//cout << *it << std::endl;
		arrBins[j]=*it;
		j++;	
	}
	// Create new rebinned histograms (only temporary)
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,"hbnew",arrBins);
	TH1F *hsnew =(TH1F*) hs->Rebin(binEdges.size()-1,"hsnew",arrBins);
	

	// Better smoothing which doesn't use the first and last binsi, performs a fit to the histogram	
	if (hsnew->Integral()!=0 && hbnew->Integral()!=0 && binEdges.size()-1 > 10){
		histogramSmoothingFit(hsnew);
		histogramSmoothingFit(hbnew);
        }


	delete [] arrBins;

	if (hbnew->Integral()==0 || hsnew->Integral()==0) return binEdges;
	if (hbnew->GetNbinsX() <= 10) return binEdges;

	int nNewBins = hbnew->GetNbinsX();

	// Here is the New Algorithm
	int 	*counters, *chosen_counters,*frozen_counters;
	double 	highestMaxSignificance=0;
	int 	chosenN=1;
	int 	*finalCounters ;

	g_step = (int)TMath::Exp(TMath::Log(nNewBins/2)/2);
	if (g_step < 1) g_step=1;
		
	for (int N=2;N<7;N++){				// Refuse to go beyond 7 Bins, will take forever
	  double maximumSignificance=0;
	  counters = new int[N];
	  chosen_counters = new int[N];
	  frozen_counters = new int[N];
	  for (int c=0;c<N;c++) counters[c]=c+2;	// init to starting values
	  for (int c=0;c<N;c++) frozen_counters[c]=c+2; // init to starting values
	  for (int c=0;c<N;c++) chosen_counters[c]=c+2; // init to starting values

	  double diff;
	  clock_t start;

	  sweepmode=0;	// First perform Broad Scan with optimized step size (g_step)
	  std::cout << "Performing Fully optimized Scan"	<<std::endl;
	  start=clock();
	  maxSigScan(&maximumSignificance,frozen_counters,chosen_counters,hsnew,hbnew,N,counters,N-1);

	  sweepmode=1;	// Now do Fine scan after having found rough maximum
	  for (int c=0;c<N;c++) counters[c]=chosen_counters[c]; // init to rough guess
          for (int c=0;c<N;c++) frozen_counters[c]=chosen_counters[c];  // init to rough guess

          // For full scanning, need to move lowest boundary to lowest point now
          int resetpoint = (2>frozen_counters[0]-g_step) ? 2 : frozen_counters[0]-g_step;
          counters[0]=resetpoint;

          maxSigScan(&maximumSignificance,frozen_counters,chosen_counters,hsnew,hbnew,N,counters,N-1);

	  diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
	  std::cout << Form("Finished, time taken = %3.5f",diff)<<std::endl;
	  std::cout << "N Bins, Max Significance -> " << N+1 << " "<<maximumSignificance << std::endl;


	  if ((maximumSignificance-highestMaxSignificance)/highestMaxSignificance > 0.001){
		highestMaxSignificance = maximumSignificance ;
	  	finalCounters= new int[N];
	  	chosenN = N;
     	  	for (int cc=0;cc<N;cc++) finalCounters[cc]=chosen_counters[cc];
	  } else {
		
		break;
	  }
	
	}

        std::vector<double> newbinEdges;
	newbinEdges.push_back(hsnew->GetBinLowEdge(1));
	for (int newguy=0;newguy<chosenN;newguy++){
		//std::cout << "newEdge = " << hsnew->GetBinLowEdge(finalCounters[newguy])<<std::endl;
	 	newbinEdges.push_back(hsnew->GetBinLowEdge(finalCounters[newguy]));
	}
	newbinEdges.push_back(hsnew->GetBinLowEdge(nNewBins+1));

	delete [] finalCounters;
	delete [] counters;
	delete [] chosen_counters;
	
	return newbinEdges;

/*
	 int nNewBins = hbnew->GetNbinsX();
         newbinEdges.push_back(hbnew->GetBinLowEdge(nNewBins+1));
         int i=nNewBins;
 	 std::vector<double> backgroundsofar,signalsofar;
         while (i>1){

                 //std::cout << "At Bin - "<< i <<std::endl;
                 int k = i-1;
                 double highEdge=hbnew->GetBinLowEdge(i);
                 double S = hsnew->GetBinContent(i);
                 double B = hbnew->GetBinContent(i);
                 double Stot =S;
                 double Btot =B;
                 if (B!=0){
                 bool carryOn=true;
 
                 while ( carryOn){
                         if (k>=1){
 
                           double S1 = hsnew->GetBinContent(k);
                           double B1 = nTargetBins;
                           double BB1 = nTargetBins;
                            if (B1==0) {
                                  carryOn=true;
                                  highEdge = hbnew->GetBinLowEdge(k);
                                  k++;
                            } else {
                              //if (compareLH(S1,Stot,B1,Btot,scaler)){
                              // now the think we want to compare to is not S1 and B1 but the rest of S1 and B1
                              double SSumRight = hsnew->Integral(1,k-1);
                              double BSumRight = hbnew->Integral(1,k-1);
 			      //double sigsofar = calculateSigMulti(signalsofar,backgroundsofar);
			   
                              if (compareLHWide(Stot,S1,SSumRight,Btot,B1,BSumRight,1.0, signalsofar,backgroundsofar)){
 
                              //cout << "Merging dude!" << std::endl; 
                              highEdge = hbnew->GetBinLowEdge(k);
 
                              Stot+=S1;
                              Btot+=B1;
 
                              carryOn = true;
                              k--;
                           } else {

				// check if we even care at this point, if we are no longer merging, check that a split is worth is
			      if (signalsofar.size()==0) {carryOn=false;}
			      else{
				std::vector<double> tmpS = signalsofar;
				std::vector<double> tmpB = backgroundsofar;
				std::vector<double> tmpS2 = signalsofar;
				std::vector<double> tmpB2 = backgroundsofar;
				tmpS.push_back(Stot+S1);
				tmpS.push_back(SSumRight);
				tmpB.push_back(Btot+B1);
				tmpB.push_back(BSumRight);
				tmpS2.push_back(Stot+S1+SSumRight);
				tmpB2.push_back(Btot+B1+BSumRight);
				if ((calculateSigMulti(tmpS,tmpB)-calculateSigMulti(tmpS2,tmpB2)) /calculateSigMulti(tmpS2,tmpB2) > 0.0001) carryOn= false;
				else{ carryOn=true; k--;}
			      }
			   }
			}
 
                         } else {
                           highEdge = hbnew->GetBinLowEdge(k+1);
                           carryOn=false;
                         }
                         }
                 }
                 newbinEdges.push_back(highEdge);
		 signalsofar.push_back(Stot);
		 backgroundsofar.push_back(Btot);
                 i=k;
         }




        reverse(newbinEdges.begin(),newbinEdges.end());	


	// now we have new Bin edges to return to the 
*/

}
// ----------------------------------------------------------------------------------------------------
std::vector<double> RooContainer::soverBOptimizedBinning(TH1F *hs,TH1F *hb,int nTargetBins,double penaltyScale){

	// Performs Optimized Binning based on a Signal and Background (S/B) distributions
	// First runs the optimizedBinning on background and rebins S and B clones, note, always performs 
	// revise_target=false,direction=-1 and use_n_entries=true
	// nTargetBins is used for the flat binning, decision to merge based on penaltyScale
	int ninitBins = hb->GetNbinsX();
	if (hs->Integral()==0 ||  hb->Integral()==0 || ninitBins < 2) {
		std::vector<double> binEdges;
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(ninitBins+1));
		return binEdges;
	}

	std::vector<double> binEdges = optimizedReverseBinning(hb,nTargetBins,false,true);

	int j =0;
	double *arrBins = new double[binEdges.size()];
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		//cout << *it << std::endl;
		arrBins[j]=*it;
		j++;	
	}
	// Create new rebinned histograms (only temporary)
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,"hbnew",arrBins);
	TH1F *hsnew =(TH1F*) hs->Rebin(binEdges.size()-1,"hsnew",arrBins);
	

	// Better smoothing which doesn't use the first and last bins	
	if (hsnew->Integral()!=0 && hbnew->Integral()!=0 && binEdges.size()-1 > 10){
		histogramSmoothing(hsnew,1000);
		histogramSmoothing(hbnew,1000);
		//hsnew->Smooth(1000);
		//hbnew->Smooth(1000);
        }
	// Do we really need the background histogram ?  we will be assuming that the first step is nentries per bin

	// Smooth signal new binned histograms, the size of smoothing should be ~1% of the total bins	
	//int nSmooth = (int) 0.01*hsnew->GetNbinsX();
	//hsnew->Smooth(nSmooth);

	delete [] arrBins;

	if (hbnew->Integral()==0 || hsnew->Integral()==0) return binEdges;
	std::vector<double> newbinEdges;
	newbinEdges.push_back(hbnew->GetBinLowEdge(1));
	int nNewBins = hbnew->GetNbinsX();
	int i=1;
	double maxSoB=0;
        for (int j=1;j<=hbnew->GetNbinsX();j++){
                double newMaximum = hsnew->GetBinContent(j)/hbnew->GetBinContent(j) ;
                if (newMaximum>maxSoB) maxSoB=newMaximum;
        }

	while (i<=nNewBins){

		int k = i+1;
		double highEdge=hbnew->GetBinLowEdge(i+1);
		double S = hsnew->GetBinContent(i);
		double B = hbnew->GetBinContent(i);
                double Stot =S;
                double Btot =B;

		if (B!=0){ 
		  bool carryOn=true;

		  while ( carryOn){

		  	double SoB = S/B;

			if (k<=nNewBins){

			  double S1 = hsnew->GetBinContent(k);
			  double B1 = hbnew->GetBinContent(k);
			  if (B1==0) {
				carryOn=true;
			      	highEdge = hbnew->GetBinLowEdge(k+1);
				Stot+=S1;
				k++;
			  }
			  else{

			    double SoB1 = S1/B1;
			    double scaler=penaltyScale;
                            double prob   = TMath::Prob((scaler*(Stot+S1))*(scaler*(Stot+S1))/(Btot+B1),1);
                            double prob1  = TMath::Prob(((scaler*Stot)*(scaler*Stot)/(Btot)) + ((scaler*S1)*(scaler*S1)/(B1)),2);
                  	    double importance = SoB/maxSoB;

//			      if (fabs(SoB-SoB1)/SoB < penaltyScale/importance ){
			      if ( prob<prob1 ){

			      highEdge = hbnew->GetBinLowEdge(k+1);
			      Stot+=S1;
			      Btot+=B1;
			      carryOn = true;
			      k++;
			    } else {carryOn=false;}
			 }

			} else {
			  highEdge = hbnew->GetBinLowEdge(k);
			  carryOn=false;
			}
		  }
		}
	        newbinEdges.push_back(highEdge);
		i=k;
	}

	// now we have new Bin edges to return to the 
	return newbinEdges;

}
// ----------------------------------------------------------------------------------------------------
std::vector<std::vector<double> > RooContainer::RebinConstantEdges(std::string datasetname, int nTargetBins){

	std::vector<std::vector<double> > return_bins;
	for (int cat=0;cat<ncat;cat++){
	   std::map<std::string,TH1F>::iterator it_th=m_th1f_.find(getcatName(datasetname,cat));
	   if (it_th!=m_th1f_.end()){
		return_bins.push_back(rebinConstantEdges(&(it_th->second),nTargetBins));

	   } else {
		std::cerr << "WARNING ! -- RooContainer::OptimizedBinning -- No such binned dataset as " << datasetname << std::endl;
	   }
	}

	return return_bins;

}
// ----------------------------------------------------------------------------------------------------
std::vector<double> RooContainer::rebinConstantEdges(TH1F *hb,int nTargetBins){
	 
	int nBins = hb->GetNbinsX();
	std::vector<double> binEdges;
	
	for (int j=1;j<=nBins;j++) {
		if (j%nTargetBins==0)	binEdges.push_back(hb->GetBinLowEdge(j+1));
	}
	return binEdges;
}
// ----------------------------------------------------------------------------------------------------
std::vector<std::vector<double> > RooContainer::OptimizedBinning(std::string datasetname, int nTargetBins,bool revise_target,bool use_n_entries,int direction){

	std::vector<std::vector<double> > return_bins;
	for (int cat=0;cat<ncat;cat++){
	   std::map<std::string,TH1F>::iterator it_th=m_th1f_.find(getcatName(datasetname,cat));
	   if (it_th!=m_th1f_.end()){
		if (direction!=-1) return_bins.push_back(optimizedBinning(&(it_th->second),nTargetBins,revise_target,use_n_entries));
		else return_bins.push_back(optimizedReverseBinning(&(it_th->second),nTargetBins,revise_target,use_n_entries));

	   } else {
		std::cerr << "WARNING ! -- RooContainer::OptimizedBinning -- No such binned dataset as " << datasetname << std::endl;
	   }
	}

	return return_bins;

}
// ----------------------------------------------------------------------------------------------------
std::vector<double> RooContainer::optimizedBinning(TH1F *hb,int nTargetBins,bool revise_target, bool use_n_target){
	// Return a set of bins which are "smoother" 

	if (revise_target) {
		if (use_n_target){
		   std::cerr << "WARNING -- RooContainer::OptimizedBinning -- Can't use number of Entries as target in revised binning algo " << std::endl; 
		   use_n_target = false;  // geometric algo always use revised number of bins, not number of entries
		
		}
	}

	int nBins = hb->GetNbinsX();
	std::vector<double> binEdges;

	double targetNumbers;
	if (use_n_target) targetNumbers = nTargetBins; 
	else targetNumbers = hb->Integral()/nTargetBins;

	if (hb->Integral() < 2*targetNumbers){
		std::cout << "RooContainer::OptimizedBinning -- Not enough entries in histogram for target numbers calculated - " 
			  << targetNumbers 
			  << ", Returning current bin boundaries "  << std::endl;
		//for (int j=2;j<=nBins+1;j++) binEdges.push_back(hb->GetBinLowEdge(j));
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(nBins+1));
		return binEdges;
	}
	binEdges.push_back(hb->GetBinLowEdge(1));

	double sumBin = 0;
	int i=1;
	while (i<=nBins){
		if (revise_target) targetNumbers = hb->Integral(i,nBins)/nTargetBins;
		sumBin=hb->GetBinContent(i);
		double highEdge=hb->GetBinLowEdge(i+1);

		bool carryOn = sumBin <= targetNumbers;
		while ( carryOn){
			if (i<nBins){
			  sumBin+=hb->GetBinContent(i+1);
			  highEdge = hb->GetBinLowEdge(i+2);
			  carryOn =(sumBin <targetNumbers && i<=nBins);
			  i++;
			} else {
			  highEdge = hb->GetBinLowEdge(i+1);
			  carryOn=false;
			}
		}
	        binEdges.push_back(highEdge);
		i++;
	}
        if (sumBin < 10) binEdges.erase(binEdges.end()-2);
	return binEdges;

}
// ----------------------------------------------------------------------------------------------------
std::vector<double> RooContainer::optimizedReverseBinning(TH1F *hb,int nTargetBins,bool revise_target, bool use_n_target){
	// Return a set of bins which are "smoother" 

	if (revise_target) {
		if (use_n_target){
		   std::cerr << "WARNING -- RooContainer::OptimizedBinning -- Can't use number of Entries as target in revised binning algo " << std::endl; 
		   use_n_target = false;  // geometric algo always use revised number of bins, not number of entries
		
		}
	}

	int nBins = hb->GetNbinsX();
	std::vector<double> binEdges;

	double targetNumbers;
	if (use_n_target) targetNumbers = nTargetBins; 
	else targetNumbers = hb->Integral()/nTargetBins;

	if (hb->Integral() < 2*targetNumbers){
		std::cout << "RooContainer::OptimizedBinning -- Not enough entries in histogram for target numbers calculated - " 
			  << targetNumbers 
			  << ", Returning current bin boundaries "  << std::endl;
		//for (int j=2;j<=nBins+1;j++) binEdges.push_back(hb->GetBinLowEdge(j));
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(nBins+1));
		return binEdges;
	}
	binEdges.push_back(hb->GetBinLowEdge(nBins+1));

	std::cout << "RooContainer::optimizedBinning -- Performing Reverse Optimize Binning" <<std::endl;
	double sumBin = 0;
	int i=nBins;
	while (i>=1){
		if (revise_target) targetNumbers = hb->Integral(1,i)/nTargetBins;
		sumBin=hb->GetBinContent(i);
		double highEdge=hb->GetBinLowEdge(i);

		bool carryOn = sumBin <= targetNumbers;
		while ( carryOn){
			if (i>1){
			  sumBin+=hb->GetBinContent(i-1);
			  highEdge = hb->GetBinLowEdge(i-1);
			  carryOn =(sumBin <targetNumbers && i>=1);
			  i--;
			} else {
			  highEdge = hb->GetBinLowEdge(i);
			  carryOn=false;
			}
		}
	        binEdges.push_back(highEdge);
		i--;
	}
        if (sumBin < 10) binEdges.erase(binEdges.end()-2);
	reverse(binEdges.begin(),binEdges.end());
	return binEdges;

}
// ----------------------------------------------------------------------------------------------------
void RooContainer::RebinBinnedDataset(std::string new_name,std::string name,std::vector<double>  catBinEdges, bool systematics){

	for (int cat=0;cat<ncat;cat++){
	  std::string catName = getcatName(name,cat);
	  std::string catNewName = getcatName(new_name,cat);
	  std::map<std::string,TH1F>::iterator it = m_th1f_.find(catName);
	  if (it!=m_th1f_.end())
		rebinBinnedDataset(catNewName,catName,&(it->second),catBinEdges);
	  else {
		std::cerr << "WARNING -- RooContainer::RebinBinnedDataset -- No Such Binned Dataset as "
			  << getcatName(name,cat)
			  << std::endl;
	  }

	  if (systematics){
		for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
		  for (int sys=1;sys<=nsigmas;sys++){
		     std::string sysDName = getsysindexName(catName,it_sys->first,sys,-1);
		     std::string sysUName = getsysindexName(catName,it_sys->first,sys,1);
		     std::string sysDNewName = getsysindexName(catNewName,it_sys->first,sys,-1);
		     std::string sysUNewName = getsysindexName(catNewName,it_sys->first,sys,1);
		     std::map<std::string,TH1F>::iterator itD = m_th1f_.find(sysDName);
		     std::map<std::string,TH1F>::iterator itU = m_th1f_.find(sysUName);
		     if (itD!=m_th1f_.end()) {
			rebinBinnedDataset(sysDNewName,sysDName,&(itD->second),catBinEdges);
			rebinBinnedDataset(sysUNewName,sysUName,&(itU->second),catBinEdges);
		     }	
		  }
	        }
	  }

	}
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::RebinBinnedDataset(std::string new_name,std::string name,std::vector <std::vector<double> > catBinEdges, bool systematics){

	for (int cat=0;cat<ncat;cat++){
	  std::string catName = getcatName(name,cat);
	  std::string catNewName = getcatName(new_name,cat);
	  std::map<std::string,TH1F>::iterator it = m_th1f_.find(catName);
	  if (it!=m_th1f_.end())
		rebinBinnedDataset(catNewName,catName,&(it->second),catBinEdges[cat]);
	  else {
		std::cerr << "WARNING -- RooContainer::RebinBinnedDataset -- No Such Binned Dataset as "
			  << getcatName(name,cat)
			  << std::endl;
	  }

	  if (systematics){
		for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
		  for (int sys=1;sys<=nsigmas;sys++){
		     std::string sysDName = getsysindexName(catName,it_sys->first,sys,-1);
		     std::string sysUName = getsysindexName(catName,it_sys->first,sys,1);
		     std::string sysDNewName = getsysindexName(catNewName,it_sys->first,sys,-1);
		     std::string sysUNewName = getsysindexName(catNewName,it_sys->first,sys,1);
		     std::map<std::string,TH1F>::iterator itD = m_th1f_.find(sysDName);
		     std::map<std::string,TH1F>::iterator itU = m_th1f_.find(sysUName);
		     if (itD!=m_th1f_.end()) {
			rebinBinnedDataset(sysDNewName,sysDName,&(itD->second),catBinEdges[cat]);
			rebinBinnedDataset(sysUNewName,sysUName,&(itU->second),catBinEdges[cat]);
		     }	
		  }
	        }
	  }

	}
}
// ----------------------------------------------------------------------------------------------------

void RooContainer::rebinBinnedDataset(std::string new_name,std::string name, TH1F *hb,std::vector<double> binEdges){

	double *arrBins = new double[binEdges.size()];
	int j=0;
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		arrBins[j]=*it;
		j++;
		
	}
	//const char *h_name = (const char *) hb->GetName;
	//const char *title  = (const char *) hb->GetTitle;
	
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,hb->GetName(),arrBins);
	hbnew->SetName(Form("th1f_%s",new_name.c_str()));
	//cout << "title for new re-binned histogream - " << hb->GetTitle()<<endl; 
	hbnew->SetTitle(hb->GetTitle());

	// Just a quick test, mask the last "channel"
	//hbnew->SetBinContent(hbnew->GetNbinsX(),0);
	//cout << "DONT DO THIS IN MAIN PROGRAM ----- LINE 1563 rebin setting last bin to 0" <<endl;
	delete [] arrBins;
	m_th1f_[new_name]=*hbnew;
	
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::composePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula,bool use_extended){

    RooArgList roo_args;
    RooArgList roo_funs;
    RooArgList normalisationTerms;
    TString expression = "@0";

    std::vector <RooRealVar*> norms;
    double check_sum = 0.;
    int arg_count    = 1;

    for (std::vector<std::string>::iterator it_fun = formula.begin()
	;it_fun != formula.end()
	;it_fun++
	){
	  if (use_extended){
	    std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(*it_fun);

	    if (exp != m_exp_.end()){
	      roo_funs.add((*exp).second);
	      roo_args.add(m_real_var_[(*it_fun)]);
              norms.push_back(&m_real_var_[(*it_fun)]);
	      normalisationTerms.add(m_real_var_[(*it_fun)]);
              if (arg_count < (int)formula.size()) expression+=Form("+@%d",arg_count);
	      std::cout << "RooContainer::ComposePdf -- Including Extended Pdf " 
		        << *it_fun << std::endl;
	      arg_count++;

	    } else {
      	      std::cerr << "WARNING -- RooContainer::ComposePdf -- No Pdf Found Named " 
	              << *it_fun << std::endl;
	      return;
	    }
	  }

	  else {
	    std::map<std::string,RooAbsPdf*>::iterator exp = m_gen_.find(*it_fun);

	    if (exp != m_gen_.end()){
	      roo_funs.add(*((*exp).second));

              if (arg_count < (int)formula.size()){

	        RooRealVar *functionNormValue = &m_real_var_[(*it_fun)];
	        addRealVar(Form("frac_%s",(*it_fun).c_str()),functionNormValue->getVal(),functionNormValue->getMin(),functionNormValue->getMax());
	        RooRealVar *fractionVar = &m_real_var_[Form("frac_%s",(*it_fun).c_str())];
	        roo_args.add(*fractionVar);
	        check_sum += fractionVar->getVal();
		//fractionVar->removeRange();
	        arg_count++;
	      }
	      std::cout << "RooContainer::ComposePdf -- Including Pdf " 
		        << *it_fun << std::endl;

	    } else {
      	      std::cerr << "WARNING -- RooContainer::ComposePdf -- No Pdf Found Named " 
	              << *it_fun << std::endl;
	      return;
	    }
	  }
	}	

    if (!use_extended){
	if (check_sum > 1.){
      	      std::cerr << "WARNING -- RooContainer::ComposePdf -- Fractions of Pdf components sum to greater than 1 !!! will not make Composed Pdf" 
	                << std::endl;
	      return;
	} else {

    	  RooAddPdf *temp = new RooAddPdf(Form("pdf_%s",name.c_str()),composition.c_str(),roo_funs,roo_args,true);
	  //temp->fixCoefRange("FullObservableRange");
     	  RooRealVar temp_var(Form("pdf_%s_norm",name.c_str()),name.c_str(),100.,0.0,1e6);
     	  m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));
          m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,temp));

     	  RooExtendPdf temp_extend(name.c_str(),name.c_str(),*temp,m_real_var_[name]);
          m_exp_.insert(std::pair<std::string,RooExtendPdf>(name,temp_extend));
        }
    }

    else{
      RooAddPdf temp(name.c_str(),composition.c_str(),roo_funs,roo_args);
      m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,&temp));
      m_pdf_.insert(std::pair<std::string,RooAddPdf>(name,temp));
      std::cerr << "EXPRESSSION - " << expression.Data()<<std::endl;;
      RooFormulaVar tmpSum(Form("pdf_%s_norm",name.c_str()),name.c_str(),expression.Data(),normalisationTerms);
      m_form_var_.insert(std::pair<std::string,RooFormulaVar >(name,tmpSum));
      m_comp_pdf_norm_[name] = norms;
    }

    std::cout << "RooContainer::ComposePdf -- Created Composed PDF from Extended pdf's called " 
	      << name << std::endl;			       
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::createDataSet(std::string name,std::string data_name,int nbins,double x1,double x2){

    std::map<std::string,RooRealVar>::const_iterator test=m_real_var_.find(name);
    if (test != m_real_var_.end()){ 

      bins_[data_name] = nbins;
      m_data_var_ptr_[data_name] = &(m_real_var_[name]);
      m_weight_var_ptr_[data_name] = &(m_real_var_[getweightName(name)]);
      
      double xmin = (test->second).getMin();
      double xmax = (test->second).getMax();
      double r1,r2;
      int number_of_bins;

      if (x1 < -990 || x2 < -990 || x1 < xmin || x2 > xmax){
	r1 = xmin;
        r2 = xmax;
      } else {
	r1 = x1;
        r2 = x2;
        //m_data_var_ptr_[data_name]->setRange(Form("%s_%f_%f",data_name.c_str(),r1,r2),r1,r2);
        //RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet((*test).second,*m_weight_var_ptr_[data_name]),WeightVar(getweightName(name).c_str()),CutRange(Form("%s_%f_%f",data_name.c_str(),r1,r2)));
        //data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));
      }

      m_var_min_.insert(std::pair<std::string, double >(data_name,r1));
      m_var_max_.insert(std::pair<std::string, double >(data_name,r2));
           
      RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet((*test).second,*m_weight_var_ptr_[data_name]),getweightName(name).c_str());
      data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));
      data_obs_names_.insert(std::pair<std::string,std::string>(data_name,name));


      number_of_bins = nbins;
      
      TH1F tmp_hist(Form("th1f_%s",data_name.c_str()),name.c_str(),number_of_bins,r1,r2);
      tmp_hist.Sumw2();
      if (blind_data) {tmp_hist.SetLineColor(0); tmp_hist.SetMarkerColor(0);} // Invisible points for bliding
      tmp_hist.GetYaxis()->SetTitle(Form("Events / (%.3f)",tmp_hist.GetBinWidth(1)));
      m_th1f_[data_name] = tmp_hist;

      if (verbosity_){
          std::cout << "RooContainer::CreateDataSet -- Created RooDataSet from " << name 
	       << " with name " << data_name <<std::endl;
          std::cout << "RooContainer::CreateDataSet -- Created TH1F from " << data_name << std::endl;
      }

    } else {
      std::cerr << "WARNING -- RooContainer::CreateDataSet -- No RealVar found Named "
                << name
		<< " CRASH Expected soon!!! -- WARNING"
		<< std::endl;
    }	
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::fitToData(std::string name_func, std::string name_data, std::string  name_var
			    ,double x1, double x2, double x3, double x4){

    bool use_composed_pdf = false;
    int n_pars;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To Dataset " 
	      << name_var
	      << std::endl; 

    std::map<std::string,RooDataSet>::iterator it_data = data_.find(name_data);
    if (it_data == data_.end()){
      std::cerr << "WARNING -- RooContainer::FitToData -- No DataSet Found Named "
 	 	<< name_data << std::endl;
      return;
    }

    RooFitResult *fit_result = new RooFitResult;
    RooAbsPdf *pdf_ptr;
    // Look in the composed pdf before checking the standards
    std::map<std::string ,RooAddPdf>::iterator it_pdf = m_pdf_.find(name_func);
    std::map<std::string ,RooExtendPdf>::iterator it_exp = m_exp_.find(name_func);

    int bins = bins_[name_data];
    double x_min = m_var_min_[name_data];
    double x_max = m_var_max_[name_data];
    int mode=0;

    RooRealVar *real_var = m_data_var_ptr_[name_data];

    if (it_pdf != m_pdf_.end()){
      use_composed_pdf = true;
      pdf_ptr = &(it_pdf->second);
     if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){ // Full Range Fit
     //   real_var->setRange("FullObservableRange",x_min,x_max);
        real_var->setRange("rnge",x_min,x_max);
        fit_result = (it_pdf->second).fitTo(it_data->second,Range("rnge"),RooFit::Extended(true),RooFit::Save(true),RooFit::Strategy(1));
        mode = 0;
     } else if (x3 < -990 && x4 < -990){ // Single window fit
      if (x1 < x_min || x2 > x_max){
     //   real_var->setRange("FullObservableRange",x_min,x_max);
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!! " 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second,RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
      //  real_var->setRange("FullObservableRange",x_min,x_max);
        real_var->setRange("rnge",x1,x2);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge"),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
        mode = 1;
      }
     } else {
      if (x1 < x_min || x4 > x_max){ // Excluded window fit
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second,RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
      //  real_var->setRange("FullObservableRange",x2,x3);
        real_var->setRange("rnge1",x1,x2);
        real_var->setRange("rnge2",x3,x4);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge1,rnge2"),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
	// Not sure what to do with combined PDF in terms of getting Nevents correct with multiple ranges
	// idea is to sum over all components and then scale by sum / sum entries
	double nEntries = (it_data->second).sumEntries(0,"rnge1") + (it_data->second).sumEntries(0,"rnge2");
	double sumComp = 0;
//	for (std::vector<RooRealVar*>::iterator it_comp_norm = m_comp_pdf_norm_[name_func].begin();it_comp_norm != m_comp_pdf_norm_[name_func].end();it_comp_norm++)
//	  {
//		sumComp += (*it_comp_norm)->getVal();
//	  }
	for (std::vector<RooRealVar*>::iterator it_comp_norm = m_comp_pdf_norm_[name_func].begin();it_comp_norm != m_comp_pdf_norm_[name_func].end();it_comp_norm++)
	  {
//		(*it_comp_norm)->setVal(nEntries/sumComp);
//		(*it_comp_norm)->setError(TMath::Sqrt((*it_comp_norm)->getVal()));
		(*it_comp_norm)->setVal(2*(*it_comp_norm)->getVal());
		(*it_comp_norm)->setError(2*(*it_comp_norm)->getError());
	  }
	
        mode = 2;
      }
     }
    }

    else {
      if (it_exp != m_exp_.end()){
      pdf_ptr = &(it_exp->second);
       if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){
     //   real_var->setRange("FullObservableRange",x_min,x_max);
        real_var->setRange("rnge",x_min,x_max);
        fit_result = (it_exp->second).fitTo(it_data->second,Range("rnge"),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
        mode = 0;
      } else if (x3 < -990 && x4 < -990){ // Single window fit
        if (x1 < x_min || x2 > x_max){
     //     real_var->setRange("FullObservableRange",x_min,x_max);
          std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo(it_data->second,RooFit::Extended(true),RooFit::Save(true),RooFit::Strategy(1));
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
     //    real_var->setRange("FullObservableRange",x_min,x_max);
         real_var->setRange("rnge",x1,x2);
         fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge"),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
         mode = 1;
        }	
       } else {
        if (x1 < x_min || x4 > x_max){
     //     real_var->setRange("FullObservableRange",x_min,x_max);
          std::cerr << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo((it_data->second),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
        //  real_var->setRange("FullObservableRange",x1,x4);
          real_var->setRange("rnge1",x1,x2);
          real_var->setRange("rnge2",x3,x4);
          fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge1,rnge2"),RooFit::Save(true),RooFit::Extended(true),RooFit::Strategy(1));

	  // This Fix just counts events in the sideband and uses that for Normalisation purposes
	  double newNorm =(it_data->second).sumEntries(0,"rnge1") + (it_data->second).sumEntries(0,"rnge2");
	  m_real_var_[name_func].setVal(newNorm);
	  m_real_var_[name_func].setError(TMath::Sqrt(newNorm));

	  // This Fix is to take account of there being 2 NLL for the 2 range fit
//	  m_real_var_[name_func].setVal(2*m_real_var_[name_func].getVal());
//	  m_real_var_[name_func].setError(2*m_real_var_[name_func].getError());
          mode = 2;
        }
       }
      }
      else
	std::cerr << "WARNING -- RooContainer::FitToData -- No Pdf Found Named "
	 	  << name_func << std::endl;
    }

    // keep hold of outer range from latest fit
    if (mode==0 or mode==1){
  	RooAbsReal* integral = pdf_ptr->createIntegral(*real_var,NormSet(*real_var),Range("rnge"));
      	latestFitRangeIntegral_[name_func] = integral;
    } else {
  	RooAbsReal* integral1 = pdf_ptr->createIntegral(*real_var,NormSet(*real_var),Range("rnge1"));
  	RooAbsReal* integral2 = pdf_ptr->createIntegral(*real_var,NormSet(*real_var),Range("rnge2"));
//  	RooAbsReal* integral1 = pdf_ptr->createIntegral(*real_var,Range("rnge1"));
//  	RooAbsReal* integral2 = pdf_ptr->createIntegral(*real_var,Range("rnge2"));
	RooFormulaVar *integral = new RooFormulaVar("sum","sum","@0+@1",RooArgSet(*integral1,*integral2));
	
        real_var->setRange("rngeiSignal",x2,x3);
	std::cout << "Number of Events in signal region = " << (it_data->second).sumEntries(0,"rngeSignal") <<std::endl;
	DUMP_[Form("partial_1_%s",name_func.c_str())]=integral1;
	DUMP_[Form("partial_2_%s",name_func.c_str())]=integral2;

      	latestFitRangeIntegral_[name_func] = integral;
    }

    RooPlot *xframe = (*real_var).frame(x_min,x_max);

    const RooCmdArg blind = (blind_data) ? RooFit::Invisible() : RooFit::MarkerColor(kBlack);
    if (bins > 0) (it_data->second).plotOn(xframe,Binning(bins),blind);
    else  (it_data->second).plotOn(xframe,blind);

    if (use_composed_pdf){
      (it_pdf->second).plotOn(xframe,LineColor(4));
      int npdfs = (it_pdf->second).pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	(it_pdf->second).plotOn(xframe
			,Components(*((it_pdf->second).pdfList().at(i)))
			,LineColor(i+1)
			,LineStyle(kDashed));
	(it_pdf->second).paramOn(xframe,RooFit::Layout(0.15,0.86,0.86));
      }
      pdf_saves_.push_back(&(it_pdf->second));
      n_pars=(it_pdf->second).getParameters(it_data->second)->getSize() -1 ;
    }

    else {
	if (mode==2){
	  (it_exp->second).plotOn(xframe,LineColor(4),RooFit::Range("rnge1,rnge2"),RooFit::Normalization(m_real_var_[name_func].getVal(),RooAbsReal::NumEvent));
	} else {
	  (it_exp->second).plotOn(xframe,LineColor(4),RooFit::Range("rnge"),RooFit::Normalization(m_real_var_[name_func].getVal(),RooAbsReal::NumEvent));
	}
	//if (mode==2){
         // real_var->setRange("rnge",x2,x3);
	 // (it_exp->second).plotOn(xframe,LineColor(2),LineStyle(3),Range("rnge")); //RooFit is so annoying that multiple goes at this ruins the normalisation
	//}
	(it_exp->second).paramOn(xframe,RooFit::Layout(0.25,0.95,0.86));
    	pdf_saves_.push_back(&(it_exp->second));
        //n_pars=(it_exp->second).getParameters(it_data->second)->getSize() -1 ;
        n_pars=(it_exp->second).getParameters(it_data->second)->getSize();
    }

    double chi_square=xframe->chiSquare(n_pars);
    if (mode == 0) {     // full range fit
      xframe->SetName(Form("%s_%s",name_func.c_str(),name_data.c_str()));
    } else if (mode == 1){ // single window fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2));
    } else {	        // sideband fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2,x3,x4));
    }

    if (fit_result->covQual() != 3) {
      std::cerr << " WARNING!!! -- RooContainer::FitToData -- Fit result covariance quality POOR from fit of " 
	   	<< name_func << " to "
	   	<< name_data << std::endl;
    }

    fit_result->printValue(std::cout);

    // -------------------------------------------------------------------------------------------
    // Make a plot of the fit
    // Dont want to let the plots pop up!
    gROOT->SetBatch(true);
    gROOT->SetStyle("Plain");
    //TLatex *text = new TLatex();
    TCanvas *can = new TCanvas(Form("plot_%s",xframe->GetName()),xframe->GetName(),1200,900) ;    
  
    can->cd(); 
    xframe->Draw();
    //text->SetNDC();
    //text->DrawLatex(0.11,0.15,Form("#chi^{2} / n d.o.f = %.3f",chi));
    fit_canvases_.push_back(can);
    fit_results_[name_func]=(fit_result); // keep a record of the last fit of a pdf
    // -------------------------------------------------------------------------------------------
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::fitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name
	     ,double x1,double x2,double x3,double x4){
  if (! fit_systematics){
    std::cerr << "WARNING!!! -- RooContainer::FitToSystematicSet -- Systematic Pdfs not created ! (Need do MakeSystematicPdfs() ) "<<
		 " No fit will be performed" << std::endl;
    return;
  }

  if (! save_systematics_data){
    std::cerr << "WARNING!!! -- RooContainer::FitToSystematicSet -- Systematic dataset was not saved ! (Need do SaveSystematicsData() ) "<<
		 " No fit will be performed" << std::endl;
    return;
  }

  for (int sys=1;sys<=nsigmas;sys++){
    fitToData(getsysindexName(name_func,sys_name,sys,-1),getsysindexName(name_var,sys_name,sys,-1),name_var,x1,x2,x3,x4);
    fitToData(getsysindexName(name_func,sys_name,sys, 1),getsysindexName(name_var,sys_name,sys, 1),name_var,x1,x2,x3,x4);
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::generateBinnedPdf(int catNo,std::string hist_nocat_name,std::string hist_name,std::string pdf_name,std::string obs_name,RooRealVar *obs,RooDataSet &dataset,int effect_type,int nbins,int mode,double x1,double x2){

  double r1,r2;
  double xmin = obs->getMin();
  double xmax = obs->getMax();
  bool multi_pdf = false;

  if (effect_type > 1 || effect_type < -1){

    std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  You have put in effect type  " << effect_type
              << " must choose from -1 : systematics from fit effect signal, 0: systematics from fit effects signal + background or 1 systematic effects background" << std::endl;  
    return;
  }

  if (x1 < -990 || x2 < -990 || x1 < xmin || x2 > xmax){
    r1 = xmin;
    r2 = xmax;    
	 
  } else {
    r1 = x1; 
    r2 = x2;
  }

  TH1F *temp_hist;
  RooAbsPdf 	 *pdf_ptr;
  std::map<std::string,RooFitResult*>::iterator res_itr = fit_results_.find(pdf_name);
  if (res_itr == fit_results_.end()){
       std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  No Fit Result Found!!! Did you fit? "
	         << pdf_name << std::endl;  
       return;
  }

  RooFitResult *res_ptr = res_itr->second;

  std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(pdf_name);

  if (exp != m_exp_.end()) {
    temp_hist = (TH1F*)((*exp).second).createHistogram(Form("th1f_%s",hist_name.c_str()),*obs,Binning(nbins,r1,r2));
    pdf_ptr = &m_exp_[pdf_name];
  } else {
    std::map<std::string,RooAddPdf>::iterator pdf  = m_pdf_.find(pdf_name);
    if (pdf != m_pdf_.end()) {
      multi_pdf = true;
      pdf_ptr = &m_pdf_[pdf_name];
      temp_hist = (TH1F*)((*pdf).second).createHistogram(Form("th1f_%s",hist_name.c_str()),*obs,Binning(nbins,r1,r2));
    }
    else {
       std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  No Pdf named "
	         << pdf_name << std::endl;  
       return;
    } 
  }


  // Need to normalize the histogram to integral under the curve
  bool external_fit;
  if (mode == 0) external_fit=true;
  if (mode == 1) external_fit=false;

  double normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf,external_fit);
  temp_hist->SetTitle(obs_name.c_str());  // Used later to keep track of the realvar associated, must be a better way to do this ?

 //  TH1F *temp_hist_up = (TH1F*)temp_hist->Clone();
 // TH1F *temp_hist_dn = (TH1F*)temp_hist->Clone();

  temp_hist->Scale(normalisation/temp_hist->Integral());
  temp_hist->SetName(Form("th1f_%s",hist_name.c_str()));
  m_th1f_.insert(std::pair<std::string,TH1F>(hist_name,*temp_hist));

  // ---------------------------------------------------------------------------------------------------------------------
  // for proper treatment of the variations of the systematic errors from the fitting in binned pdf's we need to construct
  // the set of uncorrelated parameters, 
  std::cout << "RooContainer::GenerateBinnedPdf -- Getting Latest fit result " << std::endl;
  TMatrixD cov = res_ptr->covarianceMatrix();
  TVectorD eval;
  TMatrixD evec   = cov.EigenVectors(eval);
  std::cout << "Matrix of EigenVectors from "<< pdf_name << std::endl;
  evec.Print();
  std::cout << "Eigenvalues (squares of errors)"<< std::endl;
  eval.Print();
  // ---------------------------------------------------------------------------------------------------------------------

  // now we know that the eigenvalues of the covariance matrix must scale each parameter as given by the new paraemeters
  int n_par = eval.GetNoElements();
  std::cout << "RooContainer::GenerateBinnedPdf -- Number of Parameters from the Fit -- " << n_par  << std::endl;
  RooArgSet* rooParameters = pdf_ptr->getParameters(dataset);

  // fill vector with original parameters
  std::vector <double> original_values;
  getArgSetParameters(rooParameters,original_values);

  // now loop over rows of the eigenvector Matrix
  std::vector<double> new_values;
  for (int par=0;par<n_par;par++){
    
    // this row in evec is the scales for the parameters
    std::cout << "Systematic from parameter "<< par << std::endl;

    double sigmaUnit = sigmaRange/nsigmas;
    double err = sigmaUnit*TMath::Sqrt(eval[par]);
    // If mode was 1 then the same data as was used in the fit is used in setting the limit  (ie need underconstrain)
    // if mode was 0 then we are assuming this fit is an EXTERNAL constraint to the errors from the fit must be included
    if (mode) err*=3.;  // cant go too high for fear of <ve pdfs

    std::string sys_name = (std::string)Form("%s_p%d",hist_name.c_str(),par);

    // push a new systematic into the study if this is the first category.
    std::string s_name = (std::string)Form("%s_p%d",hist_nocat_name.c_str(),par);
    if (catNo==0){
      systematics_.insert(std::pair<std::string,int>(s_name,effect_type));
    }

    for (int sys=1;sys<=nsigmas;sys++){

      new_values.clear(); // make sure its empty before we fill it
      for (int i=0;i<n_par;i++){
        new_values.push_back(original_values[i] + evec[par][i]*err);	
      }

      setArgSetParameters(rooParameters,new_values);
      normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf,external_fit);

      std::string hist_sys_name_up = getsysindexName(hist_name,s_name,sys,1);
      temp_hist = (TH1F*) pdf_ptr->createHistogram(Form("th1f_%s",hist_sys_name_up.c_str()),*obs,Binning(nbins,r1,r2));
      temp_hist->SetTitle(obs_name.c_str());  // Used later to keep track of the realvar associated, must be a better way to do this ?
      temp_hist->Scale(normalisation/temp_hist->Integral());
      temp_hist->SetName(Form("th1f_%s",hist_sys_name_up.c_str()));
      m_th1f_.insert(std::pair<std::string,TH1F>(hist_sys_name_up,*temp_hist));

      new_values.clear();  // need to clear again
      for (int i=0;i<n_par;i++){
	  new_values.push_back(original_values[i] - evec[par][i]*err);	
      }

      setArgSetParameters(rooParameters,new_values);
      normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf,external_fit);

      std::string hist_sys_name_dn = getsysindexName(hist_name,s_name,sys,-1);
      temp_hist = (TH1F*) pdf_ptr->createHistogram(Form("th1f_%s",hist_sys_name_dn.c_str()),*obs,Binning(nbins,r1,r2));
      temp_hist->SetTitle(obs_name.c_str());  // Used later to keep track of the realvar associated, must be a better way to do this ?
      temp_hist->Scale(normalisation/temp_hist->Integral());
      temp_hist->SetName(Form("th1f_%s",hist_sys_name_dn.c_str()));
      m_th1f_.insert(std::pair<std::string,TH1F>(hist_sys_name_dn,*temp_hist));
 
      // after each parameter we reset the originals back
      setArgSetParameters(rooParameters,original_values);
   }
  } 

  std::cout << "RooContainer::GenerateBinnedPdf -- Generated TH1F named  "
	    << hist_name << " From Pdf "
	    << pdf_name  << ", Normalised Histogram to " << normalisation << std::endl;      
}

// ---------------------------------------------------------------------------------------------------------------------

void RooContainer::getArgSetParameters(RooArgSet* params,std::vector<double> &val){

  val.clear();
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      //std::cout << "RooContainer:: -- Getting values from parameter -- " << std::endl;
      //rrv->Print();
      val.push_back(rrv->getVal());
  }
}

void RooContainer::setArgSetParameters(RooArgSet* params,std::vector<double> &val){

  int i = 0;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
      // each var must be shifted
      RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
      //std::cout << "RooContainer:: -- Setting values from parameter -- " << std::endl;
      //rrv->Print();
      *rrv = val[i];
      i++;
  }
}

std::pair<double,double> RooContainer::getNormalisationAndErrorFromFit(std::string pdf_name,std::string hist_name,RooAbsPdf *pdf_ptr,RooRealVar *obs,double r1,double r2,bool multi_pdf,bool external_range){

  double normalisation = 0;
  RooAbsReal* normalisationVar;
  if (multi_pdf){

	normalisationVar = &m_form_var_[pdf_name];
  }
  else{
    normalisationVar = &m_real_var_[pdf_name];
  }

  normalisation = normalisationVar->getVal();

  std::cout << "RooContainer::getNormalisationFromFit - Calculating Integral from Pdf " << pdf_name << " in range " << r1 << " -> " << r2 <<std::endl;
  obs->setRange("rngeNorm",r1,r2);
  RooAbsReal* integral = pdf_ptr->createIntegral(*obs,NormSet(*obs),Range("rngeNorm"));
  RooFormulaVar normIntVar("normIntVar","normIntVar","@0*@1/@2",RooArgSet(*normalisationVar,*integral,*latestFitRangeIntegral_[pdf_name]));

  double integralValue = integral->getVal();
  double integralError = integral->getPropagatedError(*fit_results_[pdf_name]);

  double result    = normIntVar.getVal();
  double fullError = normIntVar.getPropagatedError(*fit_results_[pdf_name]);

//  std::cout << "	Pdf Integral and Error From Latest Fit - " << integralValue << "+/-" << integralError << std::endl;
//  std::cout << "	Full (extended term) Norm - " << normalisationVar->getVal()  << std::endl;
//  std::cout << "	Integral Of latest Fit - " << latestFitRangeIntegral_[pdf_name]->getVal()  << std::endl;
//  std::cout << "	Pdf Get Norm - " << pdf_ptr->getNorm()  << std::endl;
//  std::cout << "	Normalisation In Range and Error - " << result << "+/-" << fullError << std::endl;
  return std::pair<double,double>(result,fullError);

}

double RooContainer::getNormalisationFromFit(std::string pdf_name,std::string hist_name,RooAbsPdf *pdf_ptr,RooRealVar *obs,double r1,double r2,bool multi_pdf,bool external_range){

  double normalisation = 0;
  RooAbsReal* normalisationVar;
  if (multi_pdf){

	normalisationVar = &m_form_var_[pdf_name];
  }
  else{
    normalisationVar = &m_real_var_[pdf_name];
  }

  normalisation = normalisationVar->getVal();

  std::cout << "RooContainer::getNormalisationFromFit - Calculating Integral from Pdf " << pdf_name << " in range " << r1 << " -> " << r2 <<std::endl;
  obs->setRange("rngeNorm",r1,r2);
  RooAbsReal* integral = pdf_ptr->createIntegral(*obs,NormSet(*obs),Range("rngeNorm"));
  RooFormulaVar normIntVar("normIntVar","normIntVar","@0*@1/@2",RooArgSet(*normalisationVar,*integral,*latestFitRangeIntegral_[pdf_name]));
  double result    = normIntVar.getVal();

  return result;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooDataHist(std::string name, TH1F *hist){

  // Having some problems with the RooDataHists right now in the MVA code, stick to the TH1's, also importing is very slow
  // Makes Pdf from Histogram, ie bin Integral is used not bin content

  if (save_roodatahists){

       RooDataHist tmp(Form("roohist_%s",name.c_str()),name.c_str(),RooArgList(m_real_var_[hist->GetTitle()]),hist);
       ws.import(tmp);
  }

  hist->Write();
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooPlot(RooPlot *plot,double chi){
  // Dont want to let the plots pop up!
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  // ---------------------------------
  TLatex *text = new TLatex();
  TCanvas *can = new TCanvas(Form("plot_%s",plot->GetName()),plot->GetName(),1200,900) ;    
  
  can->cd(); 
  plot->Draw();
  text->SetNDC();
  //text->DrawLatex(0.11,0.15,Form("#chi^{2} / n d.o.f = %.3f",chi));
  can->Write();
  delete can;
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::removeDuplicateElements(std::vector<RooAbsPdf*> &k){
// Useful function to remove duplicate pointers from stackoverflow
  std::vector<RooAbsPdf*>::iterator r,w;
  std::set<RooAbsPdf*> tmpset ;

  for(r=k.begin(),w=k.begin();r!=k.end();++r){
     if (tmpset.insert(*r).second ){ *w++ = *r;}
  }
        k.erase( w , k.end() );
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::histogramSmoothing(TH1F* h, int n){
   // Nothing too special, a function which will smooth a histogram but ignore the first and last
   // bins, useful for the "flat-binning" approach! 
   if (h->GetNbinsX()>3){
     int nbin = h->GetNbinsX();
     TH1F *h2 = new TH1F(Form("hn%s",h->GetName()),Form("hn%s",h->GetName()),nbin-2,0,1);
     for (int i=1;i<=nbin-2;i++){
           h2->SetBinContent(i,h->GetBinContent(i+1));
     }
     h2->Smooth(n);
     for (int i=2;i<=nbin-1;i++){
           h->SetBinContent(i,h2->GetBinContent(i-1));
     }

   }

   return;

   
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::histogramSmoothingFit(TH1F* h){
   // Nothing too special, a function which will smooth a histogram but ignore the first and last
   // bins, useful for the "flat-binning" approach! 
	float originalIntegral=h->Integral();
	if (h->GetNbinsX()>3){
	  int nbin = h->GetNbinsX();
	  TH1F *h2 = new TH1F(Form("hn%s",h->GetName()),Form("hn%s",h->GetName()),nbin-2,0,1);
	  for (int i=1;i<=nbin-2;i++){
		h2->SetBinContent(i,h->GetBinContent(i+1));
          }
	  h->Fit("pol9","F","",h->GetBinLowEdge(2),h->GetBinLowEdge(h->GetNbinsX()));
	  //h2->Smooth(n);
	  for (int i=2;i<=nbin-1;i++){
		if (h->GetFunction("pol9")->Eval(h->GetBinCenter(i-1)) < 0){
			h->SetBinContent(i,0.5*(h->GetBinContent(i-1)+h->GetBinContent(i+1)));
			
		} else {
			h->SetBinContent(i,h->GetFunction("pol9")->Eval(h->GetBinCenter(i)));
		}
          }
	
	}
	h->Scale(originalIntegral/h->Integral());
	return;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::makeSystematics(std::string observable,std::string s_name, int effect){

  for (std::map<std::string,int>::iterator isys=systematics_.begin()
      ;isys!=systematics_.end();++isys){
  
    if (isys->second != effect && isys->second != 0) continue;

    const std::string & sys_name = isys->first; 
    std::map<std::string,RooDataSet>::iterator test_it = data_.find(s_name);

    if (test_it != data_.end()){

      int bins   = bins_[s_name];

      std::vector<RooDataSet*> v_sys_up;
      std::vector<RooDataSet*> v_sys_dn;

      std::vector<TH1F*> v_th1f_up;
      std::vector<TH1F*> v_th1f_dn;

  
      for (int sys=1;sys<=nsigmas;sys++){

        std::string name_up = getsysindexName(s_name,sys_name,sys,1);
        std::string name_dn = getsysindexName(s_name,sys_name,sys,-1);
    
        TH1F *hist = &(m_th1f_[s_name]);

        double x1 = hist->GetBinLowEdge(1); 
        double x2 = hist->GetBinLowEdge(hist->GetNbinsX()+1); 

        createDataSet(observable,name_up,bins,x1,x2);
        createDataSet(observable,name_dn,bins,x1,x2);

        v_sys_up.push_back(&(data_[name_up]));
        v_sys_dn.push_back(&(data_[name_dn]));

        v_th1f_up.push_back(&(m_th1f_[name_up]));
        v_th1f_dn.push_back(&(m_th1f_[name_dn]));
      }

    std::string map_name = getsysName(s_name,sys_name);

    data_up_.insert(std::pair<std::string,std::vector<RooDataSet*> >(map_name,v_sys_up));
    data_dn_.insert(std::pair<std::string,std::vector<RooDataSet*> >(map_name,v_sys_dn));

    m_th1f_up_.insert(std::pair<std::string,std::vector<TH1F*> >(map_name,v_th1f_up));
    m_th1f_dn_.insert(std::pair<std::string,std::vector<TH1F*> >(map_name,v_th1f_dn));
  
    }
    else
      std::cerr << "RooContainer::MakeSystematics -- Cannot Create Systematics Set, "
	        << "No DATASET Found Named " << s_name << std::endl;
  }
}

void RooContainer::setAllParametersConstant() {
 
   std::map<std::string,RooRealVar>::iterator rrv = m_real_var_.begin();
   std::map<std::string,RooRealVar>::iterator rrv_end = m_real_var_.end();
   for (;rrv!=rrv_end;++rrv){
	(rrv->second).setConstant(true);
   }
//   TIterator* iter(params->createIterator());
//   for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
//    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
//    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
//   }
}
 

// ----------------------------------------------------------------------------------------------------
// A few Extra functions are desired so that RooWorkspaces can be merged:
// CategoryLooper Functions:

void RooContainer::AppendDataSet(std::string name, RooDataSet *extraData){

  // Find the DataSet
  std::map<std::string, RooDataSet>::iterator it_var  = data_.find(name);
    if (it_var == data_.end()) 
      std::cerr << "WARNING -- RooContainer::InputDataPoint -- No DataSet named "<< name << std::endl;
    else {

      (it_var->second).append(*extraData);

    }
}

void RooContainer::AppendTH1F(std::string name, TH1F *extraHist){
  
    std::map<std::string, TH1F>::iterator it_var  = m_th1f_.find(name);
    if (it_var == m_th1f_.end()) 
      std::cerr << "WARNING -- RooContainer::InputDataPoint -- No DataSet named "<< name << std::endl;
    else (it_var->second).Add(extraHist);

}

//Getty type functions
std::vector<std::string> RooContainer::GetDataSetNames(){

	std::vector<std::string> ret;
	for (std::map<std::string,RooDataSet>::iterator it_dset = data_.begin()
	    ;it_dset!=data_.end()
	    ;it_dset++){

		ret.push_back(it_dset->first);
	}
	return ret;

}
std::vector<std::string> RooContainer::GetTH1FNames(){

	std::vector<std::string> ret;
	for (std::map<std::string,TH1F>::iterator it_dset = m_th1f_.begin()
	    ;it_dset!=m_th1f_.end()
	    ;it_dset++){

		ret.push_back(it_dset->first);
	}
	return ret;

}

//EOF

