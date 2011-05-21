/****************************************************** 
   RooContainer.cc  
   Original Author - Nicholas Wardle - Imperial College
*******************************************************/
   
#include "RooContainer.h"

using namespace RooFit;

RooContainer::RooContainer(int n):ncat(n),nsigmas(30),make_systematics(false){}
// ----------------------------------------------------------------------------------------------------
void RooContainer::SetNCategories(int n){
   ncat = n;
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematicStudy(std::vector<std::string> sys_names){
   make_systematics = true;
   for (it_sys = sys_names.begin();it_sys!=sys_names.end();it_sys++)
     systematics_.push_back(*it_sys);
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::AddObservable(std::string name,float xmin,float xmax){
    addRealVar(name,xmin,xmax);
    addRealVar(getweightName(name),0.,1.0e6);
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::AddRealVar(std::string name,float init, float vmin, float vmax){
  for (int cat=0;cat<ncat;cat++){
    
addRealVar(getcatName(name,cat),init,vmin,vmax);
    if (make_systematics){
	for (it_sys=systematics_.begin(); it_sys!=systematics_.end();it_sys++){
	  for (int sys=1;sys<=nsigmas;sys++)
	    addRealVar(getsysindexName(getcatName(name,cat),*it_sys,sys,-1),init,vmin,vmax);
	  for (int sys=1;sys<=nsigmas;sys++)
	    addRealVar(getsysindexName(getcatName(name,cat),*it_sys,sys,1),init,vmin,vmax);
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
          cat_var.push_back(getsysindexName(getcatName(*it,cat),*it_sys,sys,-1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),*it_sys,sys,-1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
       } 
       for (int sys=1;sys<=nsigmas;sys++){
        std::vector<std::string> cat_var;
        for (std::vector<std::string>::iterator it=var.begin()
	  ;it!=var.end()
	  ;++it){
          cat_var.push_back(getsysindexName(getcatName(*it,cat),*it_sys,sys,1));
        } 
        addGenericPdf(getsysindexName(getcatName(name,cat),*it_sys,sys,1),formula,obs_name,cat_var,form,norm_guess,norm_min,norm_max);
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
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),*it_sys,sys,-1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),*it_sys,sys,-1),composition,cat_formula,use_extended);
      }
      for (int sys=1;sys<=nsigmas;sys++){
       std::vector<std::string> cat_formula;
       for (std::vector<std::string>::iterator it=formula.begin()
	   ;it!=formula.end()
	   ;++it){
          cat_formula.push_back(getsysindexName(getcatName(*it,cat),*it_sys,sys,1));
       }  
       composePdf(getsysindexName(getcatName(name,cat),*it_sys,sys,1),composition,cat_formula,use_extended);
      }
     }
   }
 }
}

// ----------------------------------------------------------------------------------------------------
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
       convolutePdf(getsysindexName(getcatName(name,cat),*it_sys,sys,-1),getsysindexName(getcatName(f_pdf,cat),*it_sys,sys,-1)
		   ,getsysindexName(getcatName(g_pdf,cat),*it_sys,sys,-1),(*obs).second,norm_guess);
      }
      for (int sys=1;sys<=nsigmas;sys++){
       convolutePdf(getsysindexName(getcatName(name,cat),*it_sys,sys,1),getsysindexName(getcatName(f_pdf,cat),*it_sys,sys,1)
		   ,getsysindexName(getcatName(g_pdf,cat),*it_sys,sys,1),(*obs).second,norm_guess);
      }
     }
   }
 }
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::CreateDataSet(std::string name,std::string data_name,int nbins, double x1, double x2){
  for (int cat=0;cat<ncat;cat++){
    std::string cat_name = getcatName(data_name,cat);
    createDataSet(name,cat_name,nbins,x1,x2);  
  }
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::GenerateBinnedPdf(std::string hist_name,std::string pdf_name,std::string obs_name,int nbins,double x1,double x2){

  std::map<std::string,RooRealVar>::iterator obs_var = m_real_var_.find(obs_name);
  if (obs_var == m_real_var_.end()) {
   std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  No Observable named "
	     << obs_name << std::endl;
   return;
  }
  for (int cat=0;cat<ncat;cat++){
    generateBinnedPdf(getcatName(hist_name,cat),getcatName(pdf_name,cat),obs_name,obs_var->second,nbins,x1,x2);  
  }
  
}
// ----------------------------------------------------------------------------------------------------
void RooContainer::MakeSystematics(std::string observable,std::string s_name, std::string sys_name){

  std::map<std::string,RooRealVar>::iterator obs = m_real_var_.find(observable);
  if (obs == m_real_var_.end()){
    std::cerr << "RooContainer::MakeSystematics -- WARNING!!! -- No Observable named "
	      << observable << endl;
    return;
  }

  bool in_sys = false;  
  for (it_sys=systematics_.begin();it_sys!=systematics_.end();it_sys++)
    if (sys_name == *it_sys) { in_sys=true;break;}
  if (in_sys){
    for (int cat=0;cat<ncat;cat++){
      makeSystematics(observable,getcatName(s_name,cat),sys_name);
    }
  }
  else
    std::cerr << "RooContainer::MakeSystematics -- WARNING!!! -- No systematic study available named "
	      << sys_name << endl;
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

  std::map<RooPlot*,RooFitResult*>::iterator it;
  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
      
       writeRooPlot((*it).first);
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

  std::cout << "RooContainer::Save -- Saving Histos "
            << std::endl;

  std::map<std::string,TH1F>::iterator it_d = m_th1f_.begin();
  std::map<std::string,TH1F>::iterator it_e = m_th1f_.end();

  for(;it_d != it_e; it_d++){
     writeRooDataHist((*it_d).first,&((*it_d).second));
  }
  ws.SetName("cms_hgg_workspace");
  ws.Write();
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::InputDataPoint(std::string var_name, int cat, float x, float w){
 
  if (cat>-1 && cat<ncat){
    std::string name = getcatName(var_name,cat);
    std::map<std::string, RooDataSet>::iterator it_var  = data_.find(name);
    if (it_var == data_.end()) 
      std::cerr << "WARNING -- RooContainer::InputDataPoint -- No DataSet named "<< name << std::endl;
    else{
      float min_x = m_data_var_ptr_[name]->getMin();
      float max_x = m_data_var_ptr_[name]->getMax();

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
void RooContainer::InputSystematicSet(std::string s_name, std::string sys_name,int cat
		  ,std::vector<float> x, float w){

  if (x.size() != 2*nsigmas){
    std::cerr << "WARNING -- RooContainer::InputSystematicSet -- Size of vector must be equal to "
	      << 2*nsigmas << " -- WARNING" << std::endl;
    return;	
  }
 
 else{

   if (cat>-1 && cat<ncat){

    std::string cat_name = getcatName(s_name,cat);
    std::string name = getsysName(cat_name,sys_name);

    std::map<std::string,RooDataSet>::iterator it_var  = data_.find(cat_name);

    if (it_var == data_.end()) 
      std::cerr << "WARNING -- RooContainer::InpusSystematicSet -- No DataSet named "<< cat_name << std::endl;
  
    else {

      // Safe to use this since creation of systematic set guaranteed to have come from an already existing dataset
      RooRealVar *ptr = m_data_var_ptr_[cat_name];

      std::vector<RooDataSet*>::iterator data_set_up = data_up_[name].begin();
      std::vector<RooDataSet*>::iterator data_set_dn = data_dn_[name].begin();
      std::vector<TH1F*>::iterator th1f_set_up = m_th1f_up_[name].begin();
      std::vector<TH1F*>::iterator th1f_set_dn = m_th1f_dn_[name].begin();

      std::vector<float>::iterator val = x.begin();
 
      float min_x = ptr->getMin();
      float max_x = ptr->getMax();
      
      // Loop over the first nsigmas elements as the -1 -> dn sys
      for (;data_set_dn != data_dn_[name].end()
	   ;data_set_dn++,th1f_set_dn++,val++){

        if (*val > min_x && *val < max_x){
           *ptr = *val;
           (*data_set_dn)->add(RooArgSet(*ptr),w);
           (*th1f_set_dn)->Fill(*val,w);
        }
      }

      // Loop over the second nsigmas elements as the +1 -> up sys
      for (;data_set_up != data_up_[name].end()
	   ;data_set_up++,th1f_set_up++,val++){

        if (*val > min_x && *val < max_x){
           *ptr = *val;
           (*data_set_up)->add(RooArgSet(*ptr),w);
           (*th1f_set_up)->Fill(*val,w);
        }
      }
    }
   }
   else {
    std::cerr << "WARNING -- RooContainer::InputDataPoint -- No Category Number " << cat 
              << ", category must be from 0 to " << ncat-1
	      << std::endl;
   }
  }
}

void RooContainer::WriteDataCard(std::string filename,std::string data_name
				,std::string sig_name, std::string bkg_name){

   std::vector<double> data_norms,sig_norms,bkg_norms;
   data_norms.resize(ncat); sig_norms.resize(ncat); bkg_norms.resize(ncat);

   for (int cat=0;cat<ncat;cat++){
     std::map<std::string,TH1F>::iterator it_data = m_th1f_.find(getcatName(data_name,cat));
     std::map<std::string,TH1F>::iterator it_sig  = m_th1f_.find(getcatName(sig_name,cat));
     std::map<std::string,TH1F>::iterator it_bkg  = m_th1f_.find(getcatName(bkg_name,cat));

     if (it_data !=m_th1f_.end() && it_sig !=m_th1f_.end() && it_bkg !=m_th1f_.end() ){

       data_norms[cat] = (it_data->second).Integral();
       sig_norms[cat]  = (it_sig->second).Integral();
       bkg_norms[cat]  = (it_bkg->second).Integral();   
     } else {

	std::cerr << "WARNING -- RooContainer::WriteDataCard -- Cannot find one of "
		  << data_name << " " << sig_name << " " << " " << bkg_name
		  << " DataCard will NOT be Writen -- WARNING!" << std::endl;
	return;
     }
   }

   ofstream file (Form("cms-hgg-%s-datacard.txt",sig_name.c_str()));
   file << "CMS-HGG DataCard\n";
   file << "Run with: combine cms-hgg-datacard -M Routine -D "<< data_name; 
   file << "\n---------------------------------------------\n";
   file << "imax *\n";
   file << "jmax *\n";
   file << "kmax *\n";
   file << "---------------------------------------------\n";
   file << "shapes * * "<<filename<<" th1f_$PROCESS_$CHANNEL\n";
   file << "---------------------------------------------\n";

   // Write the data info	 
   file << "bin	           ";
   for (int cat=0;cat<ncat;cat++) file << "cat"<<cat << " ";
   file << "\nobservation  ";
   for (int cat=0;cat<ncat;cat++) file << data_norms[cat] << " ";
   file << "\n---------------------------------------------\n";
   // Write the model info	 
   file << "bin	       ";
   for (int cat=0;cat<ncat;cat++) file << "cat"<<cat << "  cat" <<cat<<" ";
   file << "\nprocess  ";
   for (int cat=0;cat<ncat;cat++) file << sig_name<<" "<<bkg_name<< " ";
   file << "\nprocess  ";
   for (int cat=0;cat<ncat;cat++) file << " -1   1";
   file << "\nrate     ";
   for (int cat=0;cat<ncat;cat++) file << sig_norms[cat]<<" "<<bkg_norms[cat]<<" ";
   file << "\n---------------------------------------------\n";

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
    char_string = Form("%s_%s-up%.2d_sigma",name.c_str(),sys_name.c_str(),sys);
  else
    char_string = Form("%s_%s-down%.2d_sigma",name.c_str(),sys_name.c_str(),sys);
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
void RooContainer::addRealVar(std::string name ,float xmin,float xmax){

  RooRealVar temp(name.c_str(),name.c_str(),xmin,xmax);

  m_real_var_.insert(pair<std::string,RooRealVar >(name,temp));
  m_var_min_.insert(pair<std::string, float >(name,xmin));
  m_var_max_.insert(pair<std::string, float >(name,xmax));

  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl; 
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::addRealVar(std::string name ,float init, float vmin, float vmax){
  RooRealVar temp(name.c_str(),name.c_str(),init,vmin,vmax);
  m_real_var_.insert(pair<std::string,RooRealVar>(name,temp));
  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
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
     std::vector<std::string>::iterator it_var = var.begin();

     for (;it_var != var.end();it_var++){ 
        std::map<std::string,RooRealVar>::iterator it_check = m_real_var_.find(*it_var);
	if (it_check == m_real_var_.end()){
           std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	             << *it_var << std::endl; 
	   return;
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

	     std::map<std::string,RooRealVar>::iterator real_var =  m_real_var_.find(*it_var);
	     if (real_var != m_real_var_.end())
	       roo_args.add((*real_var).second);
	     else 
      		std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Variable Found Named " 
	        	  << *it_var << std::endl;
  	   }

      std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	        << std::endl;

      temp_1 = new RooGenericPdf(Form("pdf_%s",name.c_str()),name.c_str(),formula.c_str(),roo_args);	
      //v_gen_.push_back(temp_1);

     } else if (form == 1) { //RooExponential  - x,slope
	if (var.size() == 1){
	  temp_1 = new RooExponential(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 1 arguments for RooExponential, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 2) { //RooGaussian - x,mean,sigma
	if (var.size() == 2){
	  temp_1 = new RooGaussian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 2 arguments for RooGaussian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 3) { //RooBreitWigner - x,centre,width
	if (var.size() == 2){
	  temp_1 = new RooBreitWigner(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 2 arguments for RooBreitWigner, was given: "

		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 4) { //RooCBShape - x,centre,sigma,slope,n
	if (var.size() == 4){
	  temp_1 = new RooCBShape(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]
				 ,m_real_var_[var[2]],m_real_var_[var[3]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 4 arguments for RooCBShape, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     } else if (form == 5) { //RooVoigtian - x,centre,width,sigma
	if (var.size() == 3){
	  temp_1 = new RooVoigtian(Form("pdf_%s",name.c_str()),name.c_str(),(*obs_real_var).second,m_real_var_[var[0]],m_real_var_[var[1]]
				  ,m_real_var_[var[2]]);
          //v_gen_.push_back(temp_1);
	} else {
		
          std::cerr << "WARNING -- RooContainer::AddGenericPdf -- Need 3 arguments for RooVoigtian, was given: "
		    << var.size() << " -- WARNING"
	            << std::endl;
	  return;
	}
     }

     else {
	
       std::cerr << "WARNING -- RooContainer::AddGenericPdf -- No Mode " << form
		 << "Understood -- WARNING"
	         << std::endl;
       return;
     }

     m_gen_.insert(std::pair<std::string,RooAbsPdf*>(name,temp_1));

     RooRealVar temp_var(Form("norm_%s",name.c_str()),name.c_str(),norm_guess,0.0,1e6);
     m_real_var_.insert(std::pair<std::string,RooRealVar>(name,temp_var));

     RooExtendPdf temp(name.c_str(),name.c_str(),*temp_1,m_real_var_[name]);
     m_exp_.insert(std::pair<std::string,RooExtendPdf>(name,temp));

     std::cout << "RooContainer::AddGenericPdf -- Made extended PDF " 
	       << name << std::endl;	
    }
}

// ----------------------------------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------------------------------
void RooContainer::composePdf(std::string name, std::string  composition
			     ,std::vector<std::string> & formula,bool use_extended){

    RooArgList roo_args;
    RooArgList roo_funs;
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
      m_pdf_.insert(pair<std::string,RooAddPdf>(name,temp));
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
        RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet((*test).second,*m_weight_var_ptr_[data_name]),getweightName(name).c_str());
        data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));

      } else {
	r1 = x1;
        r2 = x2;
        m_data_var_ptr_[data_name]->setRange(Form("%s_%f_%f",data_name.c_str(),r1,r2),r1,r2);
        RooDataSet data_tmp(data_name.c_str(),data_name.c_str(),RooArgSet((*test).second,*m_weight_var_ptr_[data_name]),WeightVar(getweightName(name).c_str()),CutRange(Form("%s_%f_%f",data_name.c_str(),r1,r2)));
        data_.insert(std::pair<std::string,RooDataSet>(data_name,data_tmp));
      }


      if (nbins < -1)
        number_of_bins = (int) r2-r1;
      else
	number_of_bins = nbins;
      
      TH1F tmp_hist(Form("th1f_%s",data_name.c_str()),name.c_str(),number_of_bins,r1,r2);
      tmp_hist.GetYaxis()->SetTitle(Form("Events / (*.43f)",tmp_hist.GetBinWidth(1)));
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
    double chi_square;
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

    RooFitResult *fit_result;
    // Look in the composed pdf before checking the standards
    std::map<std::string ,RooAddPdf>::iterator it_pdf = m_pdf_.find(name_func);
    std::map<std::string ,RooExtendPdf>::iterator it_exp = m_exp_.find(name_func);

    int bins = bins_[name_data];
    float x_min = m_data_var_ptr_[name_data]->getMin();
    float x_max = m_data_var_ptr_[name_data]->getMax();
    int mode=0;

    RooRealVar *real_var = m_data_var_ptr_[name_data];

    if (it_pdf != m_pdf_.end()){
      use_composed_pdf = true;

     if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){ // Full Range Fit
        fit_result = (it_pdf->second).fitTo(it_data->second);
        mode = 0;
     } else if (x3 < -990, x4 < -990){ // Single window fit
      if (x1 < x_min || x2 > x_max){
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range" 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second);
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
        real_var->setRange("rnge",x1,x2);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge"));
        mode = 1;
      }
     } else {
      if (x1 < x_min || x4 > x_max){ // Excluded window fit
        std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range" 
		  << std::endl;
        fit_result = (it_pdf->second).fitTo(it_data->second);
        std::cerr << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
        mode = 0;
      } else {
        real_var->setRange("rnge1",x1,x2);
        real_var->setRange("rnge2",x3,x4);
        fit_result = (it_pdf->second).fitTo((it_data->second),Range("rnge1,rnge2"));
        mode = 2;
      }
     }
    }

    else {
      if (it_exp != m_exp_.end()){
       if (x1 < -990. && x2 < -990.  && x3 < -990. && x4 < -990.){
        fit_result = (it_exp->second).fitTo(it_data->second);
        mode = 0;
      } else if (x3 < -990, x4 < -990){ // Single window fit
        if (x1 < x_min || x2 > x_max){
          std::cerr << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo(it_data->second);
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
         real_var->setRange("rnge",x1,x2);
         fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge"));
         mode = 1;
        }	
       } else {
        if (x1 < x_min || x4 > x_max){
          std::cerr << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range! Will Fit to full range" 
		    << std::endl;
          fit_result = (it_exp->second).fitTo((it_data->second));
          std::cerr << " Fitted To Full Range -- WARNING!!" 
		    << std::endl;
          mode = 0;
        } else {
          real_var->setRange("rnge1",x1,x2);
          real_var->setRange("rnge2",x3,x4);
          fit_result = (it_exp->second).fitTo((it_data->second),Range("rnge1,rnge2"));
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
    }

    else {
	(it_exp->second).plotOn(xframe,LineColor(4));
	(it_exp->second).paramOn(xframe);
    	pdf_saves_.push_back(&(it_exp->second));
    }

    if (mode == 0)      // full range fit
      xframe->SetName(Form("%s_%s",name_func.c_str(),name_data.c_str()));
    else if (mode == 1) // single window fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2));
    else	        // sideband fit
      xframe->SetName(Form("%s_%s_%.1f_%.1f_%.1f_%.1f",name_func.c_str(),name_data.c_str(),x1,x2,x3,x4));

    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
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
void RooContainer::generateBinnedPdf(std::string hist_name,std::string pdf_name,std::string obs_name,RooRealVar &obs,int nbins,double x1,double x2){

  double r1,r2;
  double xmin = obs.getMin();
  double xmax = obs.getMax();
  bool normalise_full_range = true;

  if (x1 < -990 || x2 < -990 || x1 < xmin || x2 > xmax){
    r1 = xmin;
    r2 = xmax;    
	 
  } else {
    r1 = x1; 
    r2 = x2;
    normalise_full_range = false;
  }

  TH1F *temp_hist;
  RooAbsPdf *pdf_ptr;
  std::map<std::string,RooExtendPdf>::iterator exp = m_exp_.find(pdf_name);
  if (exp != m_exp_.end()) {
    temp_hist = (TH1F*)((*exp).second).createHistogram(Form("th1f_%s",hist_name.c_str()),obs,Binning(nbins,r1,r2));
    pdf_ptr = &(exp->second);
  } else {
    std::map<std::string,RooAddPdf>::iterator pdf  = m_pdf_.find(pdf_name);
    if (pdf != m_pdf_.end()) {
      temp_hist = (TH1F*)((*pdf).second).createHistogram(Form("th1f_%s",hist_name.c_str()),obs,Binning(nbins,r1,r2));
      pdf_ptr = &(pdf->second);
    }
    else {
       std::cerr << "WARNING!!! -- RooContainer::GenerateBinnedPdf --  No Pdf named "
	         << pdf_name << endl;  
       return;
    } 
  }

  // Need to normalize the histogram to integral under the curve
  double normalisation = m_real_var_[pdf_name].getVal();
  if (!normalise_full_range) {

    obs.setRange(hist_name.c_str(),x1,x2);
    RooAbsReal*integral = pdf_ptr->createIntegral(obs,NormSet(obs),Range(hist_name.c_str()));
    normalisation *= integral->getVal();
  }

  temp_hist->Scale(normalisation/temp_hist->Integral());
  temp_hist->SetName(Form("th1f_%s",hist_name.c_str()));

  temp_hist->SetTitle(obs_name.c_str());  // Used later to keep track of the realvar associated, must be a better way to do this ?

  m_th1f_.insert(std::pair<std::string,TH1F>(hist_name,*temp_hist));
  std::cout << "RooContainer::GenerateBinnedPdf -- Generated TH1F named  "
	    << hist_name << " From Pdf "
	    << pdf_name  << std::endl;      
}

// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooDataHist(std::string name, TH1F *hist){

  RooDataHist tmp(Form("roohist_%s",name.c_str()),name.c_str(),RooArgList(m_real_var_[hist->GetTitle()]),hist);

  ws.import(tmp);
  hist->Write();
}


// ----------------------------------------------------------------------------------------------------
void RooContainer::writeRooPlot(RooPlot *plot){
  // Dont want to let the plots pop up!
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  // ---------------------------------
  TCanvas *can = new TCanvas(Form("plot_%s",plot->GetName()),plot->GetName(),900,600) ;    
  can->cd(); plot->Draw();
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
void RooContainer::makeSystematics(std::string observable,std::string s_name, std::string sys_name){

  std::map<std::string,RooDataSet>::iterator test_it = data_.find(s_name);

  if (test_it != data_.end()){

    std::vector<RooDataSet*> v_sys_up;
    std::vector<RooDataSet*> v_sys_dn;

    std::vector<TH1F*> v_th1f_up;
    std::vector<TH1F*> v_th1f_dn;

    int bins   = bins_[observable];
  
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
// ----------------------------------------------------------------------------------------------------
//EOF

