/****************************************************** 
   RooContainer.cc  
   Original Author - Nicholas Wardle - Imperial College
*******************************************************/
   
#include "RooContainer.h"

using namespace RooFit;

RooContainer::RooContainer(int n, int s):ncat(n),nsigmas(s),make_systematics(false){}
// ----------------------------------------------------------------------------------------------------
void RooContainer::SetNCategories(int n){
   ncat = n;
}
void RooContainer::AddGlobalSystematic(std::string name,double val_sig, double val_bkg){
  global_systematics_[name] = std::pair<double,double>(val_sig,val_bkg);
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
  ws.import(temp);
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddObservable(std::string name,double xmin,double xmax){
    addRealVar(name,xmin,xmax);
    addRealVar(getweightName(name),0.,1.0e6);
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddRealVar(std::string name,double init,double xmin, double xmax){
  for (int cat=0;cat<ncat;cat++){
    
    addRealVar(getcatName(name,cat),init,xmin,xmax);
    if (make_systematics){
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
    if (make_systematics){
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

    if (make_systematics){
      
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

    if (make_systematics){
      
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
	
    if (make_systematics){	
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
	
    if (make_systematics){	
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
	      << observable << endl;
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

  std::map<RooPlot*,double>::iterator it;
  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
      
       writeRooPlot((*it).first,(*it).second);
  }
  
  std::cout << "RooContainer::Save -- Saving Pdfs "
            << std::endl;

  // check first and remove duplicates of the pointers. 
  removeDuplicateElements(pdf_saves_);

  std::vector<RooAbsPdf*>::iterator it_pdf;

  for(it_pdf  = pdf_saves_.begin()
     ;it_pdf != pdf_saves_.end()
     ;it_pdf++ ){
     
       ws.import(**it_pdf);
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

// ----------------------------------------------------------------------------------------------------
void RooContainer::InputDataPoint(std::string var_name, int cat, double x, double w){
 
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
void RooContainer::InputSystematicSet(std::string s_name, std::string sys_name, std::vector<int> cats
				      ,std::vector<double> x, std::vector<double> weights){

  if( weights.empty() ) { weights.resize(x.size(),1.); }
  
  assert( x.size() == cats.size() );
  assert( x.size() == weights.size() );

  if (x.size() != 2*nsigmas){
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
	  int ishift = istep < nsigmas ? istep - nsigmas : istep - nsigmas + 1;

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
	    	*ptr = val;
	    	data_set.add(RooArgSet(*ptr),w);
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
     file.open(Form("cms-hgg-datacard_%s_parBKG.txt",filename.c_str()));
     file << "CMS-HGG DataCard for Binned Limit Setting with RooDataHist+Parameterised Background\n";
   }
   else{
     file.open(Form("cms-hgg-datacard_%s.txt",filename.c_str()));
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
          std::string bcatName = getcatName(bkg_name,cat);
          std::string dcatName = getcatName(data_name,cat);
	  std::string rngeName = getcatName("datacard_",cat);
	  double r1 = m_var_min_[dcatName] ;
	  double r2 = m_var_max_[dcatName] ;
	  double norm = getNormalisationFromFit(bcatName,rngeName,m_gen_[bcatName],obs,r1,r2,multi_pdf);
	  file << " -1 " << norm << " ";   
	}
   }
   file << "\n---------------------------------------------\n";

   // Now write the systematics lines:
   // first Global Systematics
   for (std::map<std::string,std::pair<double,double> >::iterator it_g_sys = global_systematics_.begin();it_g_sys!=global_systematics_.end(); it_g_sys++){ 
     file << it_g_sys->first << " lnN  ";
     for (int cat=0;cat<ncat;cat++) file << (it_g_sys->second).first << " " << (it_g_sys->second).second << " ";	
     file <<endl;
   }

  
    
   double sigmaUnitInv = 1./(sigmaRange/nsigmas);
   for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){ 
    
     file << it_sys->first << " shape  ";
     if (it_sys->second == 0) // both background and signal effected
      for (int cat=0;cat<ncat;cat++) file << sigmaUnitInv<< " " << sigmaUnitInv<< " ";
     else if (it_sys->second == -1) // signal effected
      for (int cat=0;cat<ncat;cat++) file << sigmaUnitInv <<" 0 ";
     else if (it_sys->second == 1) // background effected
      for (int cat=0;cat<ncat;cat++) file << " 0 "<<sigmaUnitInv;
     file << endl;
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
  m_real_var_.insert(pair<std::string,RooRealVar >(name,temp));

  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl; 
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,double init,double xmin,double xmax){
  RooRealVar temp(name.c_str(),name.c_str(),init,xmin,xmax);
 // temp.removeRange();
  temp.setRange(-10.,10.);
  m_real_var_.insert(pair<std::string,RooRealVar>(name,temp));
  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl;
  
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::addFormulaVar(std::string name ,std::string formula, std::string var){

  RooFormulaVar temp(name.c_str(),name.c_str(),formula.c_str(),RooArgList(m_real_var_[var]));
  m_form_var_.insert(pair<std::string,RooFormulaVar >(name,temp));

  std::cout << "RooContainer::AddFormulaVar -- Appended the variable " 
	    << name <<std::endl; 
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
     } else if (form>60 && form < 70) { // RooChebychev - x, p1....pform-1

	if (var.size() == form-60){
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

     } else if (form>70) { // RooBernstein - x, p1....pform-1

	if (var.size() == form-70){
      	  RooArgList roo_args;
	  roo_args.add(RooConst(1.0));
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

      	   temp_1 = new RooBernstein(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,roo_args);
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

     RooRealVar temp_var(Form("%s_norm",name.c_str()),name.c_str(),norm_guess,norm_min,norm_max);
     m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));

     RooExtendPdf temp(name.c_str(),name.c_str(),*temp_1,m_real_var_[name]);
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

    m_real_var_.insert(pair<std::string,RooRealVar>(name,temp_var));

    RooExtendPdf temp(name.c_str(),name.c_str(),*convoluted,m_real_var_[name]);
    m_exp_.insert(pair<std::string,RooExtendPdf>(name,temp));
   
    std::cout << "RooContainer::ConvolvePdf -- Made convoluted PDF " 
	      << name << std::endl;	
    
}
*/
// ----------------------------------------------------------------------------------------------------
void RooContainer::composePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula,bool use_extended){

    RooArgList roo_args;
    RooArgList roo_funs;
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
	      std::cout << "RooContainer::ComposePdf -- Including Extended Pdf " 
		        << *it_fun << std::endl;

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

              if (arg_count < formula.size()){

	        roo_args.add(m_real_var_[(*it_fun)]);
	        check_sum += m_real_var_[*it_fun].getVal();
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

    	  RooAddPdf *temp = new RooAddPdf(Form("pdf_%s",name.c_str()),composition.c_str(),roo_funs,roo_args);
     	  RooRealVar temp_var(Form("norm_%s",name.c_str()),name.c_str(),100.,0.0,1e6);
     	  m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));
          m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,temp));

     	  RooExtendPdf temp_extend(name.c_str(),name.c_str(),*temp,m_real_var_[name]);
          m_exp_.insert(pair<std::string,RooExtendPdf>(name,temp_extend));
        }
    }

    else{
      RooAddPdf temp(name.c_str(),composition.c_str(),roo_funs,roo_args);
      m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,&temp));
      m_pdf_.insert(pair<std::string,RooAddPdf>(name,temp));
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

      m_var_min_.insert(pair<std::string, double >(data_name,r1));
      m_var_max_.insert(pair<std::string, double >(data_name,r2));
           
      RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet((*test).second,*m_weight_var_ptr_[data_name]),getweightName(name).c_str());
      data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));
      data_obs_names_.insert(std::pair<std::string,std::string>(data_name,name));


      number_of_bins = nbins;
      
      TH1F tmp_hist(Form("th1f_%s",data_name.c_str()),name.c_str(),number_of_bins,r1,r2);
      tmp_hist.GetYaxis()->SetTitle(Form("Events / (%.3f)",tmp_hist.GetBinWidth(1)));
      m_th1f_[data_name] = tmp_hist;

      cout << "RooContainer::CreateDataSet -- Created RooDataSet from " << name 
	   << " with name " << data_name <<endl;
      cout << "RooContainer::CreateDataSet -- Created TH1F from " << data_name << endl;

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

     if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){ // Full Range Fit
        real_var->setRange("rnge",x_min,x_max);
        fit_result = (it_pdf->second).fitTo(it_data->second,Range("rnge"),RooFit::Save(true));
        mode = 0;
     } else if (x3 < -990 && x4 < -990){ // Single window fit
      if (x1 < x_min || x2 > x_max){
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!! " 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second,RooFit::Save(true),"r");
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
        real_var->setRange("rnge",x1,x2);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge"),RooFit::Save(true),"r");
        mode = 1;
      }
     } else {
      if (x1 < x_min || x4 > x_max){ // Excluded window fit
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second,RooFit::Save(true),"rq");
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
        real_var->setRange("rnge1",x1,x2);
        real_var->setRange("rnge2",x3,x4);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge1,rnge2"),RooFit::Save(true),"r");
        mode = 2;
      }
     }
    }

    else {
      if (it_exp != m_exp_.end()){
       if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){
        real_var->setRange("rnge",x_min,x_max);
        fit_result = (it_exp->second).fitTo(it_data->second,Range("rnge"),RooFit::Save(true),"r");
        mode = 0;
      } else if (x3 < -990 && x4 < -990){ // Single window fit
        if (x1 < x_min || x2 > x_max){
          std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo(it_data->second,RooFit::Save(true));
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
         real_var->setRange("rnge",x1,x2);
         fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge"),RooFit::Save(true),"r");
         mode = 1;
        }	
       } else {
        if (x1 < x_min || x4 > x_max){
          std::cerr << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range, Normalisation may be wrong!!!" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo((it_data->second),RooFit::Save(true),"rq");
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
          real_var->setRange("rnge1",x1,x2);
          real_var->setRange("rnge2",x3,x4);
          fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge1,rnge2"),RooFit::Save(true),"rq");
          mode = 2;
        }
       }
      }
      else
	std::cerr << "WARNING -- RooContainer::FitToData -- No Pdf Found Named "
	 	  << name_func << std::endl;
    }

    RooPlot *xframe = (*real_var).frame(x_min,x_max);

    if (bins > 0) (it_data->second).plotOn(xframe,Binning(bins));
    else  (it_data->second).plotOn(xframe);

    if (use_composed_pdf){
      (it_pdf->second).plotOn(xframe,LineColor(4));
      int npdfs = (it_pdf->second).pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	(it_pdf->second).plotOn(xframe
			,Components(*((it_pdf->second).pdfList().at(i)))
			,LineColor(i+1)
			,LineStyle(kDashed));
	(it_pdf->second).paramOn(xframe);
      }
      pdf_saves_.push_back(&(it_pdf->second));
      n_pars=(it_pdf->second).getParameters(it_data->second)->getSize() -1 ;
    }

    else {
	(it_exp->second).plotOn(xframe,LineColor(4));
	(it_exp->second).paramOn(xframe);
    	pdf_saves_.push_back(&(it_exp->second));
        n_pars=(it_exp->second).getParameters(it_data->second)->getSize() -1 ;
    }

    double chi_square = xframe->chiSquare(n_pars);
    if (mode == 0)      // full range fit
      xframe->SetName(Form("%s_%s",name_func.c_str(),name_data.c_str()));
    else if (mode == 1) // single window fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2));
    else	        // sideband fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2,x3,x4));

    fit_res_.insert(std::pair<RooPlot*,double>(xframe,chi_square));
    fit_results_[name_func]=(fit_result); // keep a record of the last fit of a pdf
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::fitToSystematicSet(std::string name_func,std::string name_var
	     ,std::string sys_name
	     ,double x1,double x2,double x3,double x4){
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
              << " must choose from -1 : systematics from fit effect signal, 0: systematics from fit effects signal + background or 1 systematic effects background" << endl;  
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
	         << pdf_name << endl;  
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
	         << pdf_name << endl;  
       return;
    } 
  }


  // Need to normalize the histogram to integral under the curve
  // 2 modes - full observable range from fit and histogram (default, no range arguments)
  // 0 	     - assumes same data as fit to will be used to set limit -> make up some bullshit constraints weaker than the data
  // 1       - assumes External constraints to data being input, ie actual fit errors will be used

  double normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf);
  temp_hist->SetTitle(obs_name.c_str());  // Used later to keep track of the realvar associated, must be a better way to do this ?

  TH1F *temp_hist_up = (TH1F*)temp_hist->Clone();
  TH1F *temp_hist_dn = (TH1F*)temp_hist->Clone();
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
    cout << "Systematic from parameter "<< par << endl;

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
      normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf);

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
      normalisation = getNormalisationFromFit(pdf_name,hist_name,pdf_ptr,obs,r1,r2,multi_pdf);

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

double RooContainer::getNormalisationFromFit(std::string pdf_name,std::string hist_name,RooAbsPdf *pdf_ptr,RooRealVar *obs,double r1,double r2,bool multi_pdf){

  double normalisation = 0;
  if (multi_pdf){

     std::map<std::string,std::vector<RooRealVar*> >::iterator it = m_comp_pdf_norm_.find(pdf_name);
     // get the values which correspond to normalisations of each component pdf
     for ( std::vector<RooRealVar*>::iterator it_r = (it->second).begin() ; it_r!=(it->second).end();it_r++){
	normalisation += (*it_r)->getVal();
     }
  }

  else{
    normalisation = m_real_var_[pdf_name].getVal();
  }

  obs->setRange(hist_name.c_str(),r1,r2);
  RooAbsReal* integral = pdf_ptr->createIntegral(*obs,NormSet(*obs),Range(hist_name.c_str()));
  normalisation *= integral->getVal();


  return normalisation;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooDataHist(std::string name, TH1F *hist){

  RooDataHist tmp(Form("roohist_%s",name.c_str()),name.c_str(),RooArgList(m_real_var_[hist->GetTitle()]),hist);

  ws.import(tmp);
  hist->Write();
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooPlot(RooPlot *plot,double chi){
  // Dont want to let the plots pop up!
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  // ---------------------------------
  TLatex *text = new TLatex();
  TCanvas *can = new TCanvas(Form("plot_%s",plot->GetName()),plot->GetName(),900,600) ;    
  
  can->cd(); 
  plot->Draw();
  text->SetNDC();
  text->DrawLatex(0.8,0.5,Form("#Chi^{2} = %.3f",chi));
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
//EOF

