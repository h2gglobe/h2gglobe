#include "RooContainer.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include <cmath>
#include <iostream>

using namespace RooFit;


void RooContainer::AddRealVar(const char* name ,float xmin,float xmax){
  RooRealVar temp(name,name,xmin,xmax);

  m_real_var_.insert(pair<const char*,RooRealVar>(name,temp));
  m_var_min_.insert(pair<const char*, float>(name,xmin));
  m_var_max_.insert(pair<const char*, float>(name,xmax));
  
}

void RooContainer::AddRealVar(const char* name ,float init, float vmin, float vmax){
  RooRealVar temp(name,name,init,vmin,vmax);
  m_real_var_.insert(pair<const char*,RooRealVar>(name,temp));
  std::cout << "RooContainer::AddRealVar -- Appended the variable " 
	    << name <<std::endl;
  
}

void RooContainer::AddGenericPdf(const char* name,const char* formula
				,std::vector<const char*> & var
				,double norm_guess ){

    RooArgList roo_args;

    for (std::vector<const char*>::iterator it_var = var.begin()
	;it_var != var.end()
	;it_var++
	){
	  std::cout << "RooContainer::AddGenericPdf -- Adding Parameter " 
		    << *it_var << std::endl;
	  roo_args.add(m_real_var_[*it_var]);
	}



    std::cout << "RooContainer::AddGenericPdf -- Added all variables" 
	      << std::endl;
    RooGenericPdf temp_1(Form("comp_%s",name),name,formula,roo_args);	
    m_gen_.insert(std::pair<const char *, RooGenericPdf>(name,temp_1));

    RooRealVar temp_var(Form("norm_%s",name),name,norm_guess,0.0,10000);
    m_real_var_.insert(pair<const char*,RooRealVar>(name,temp_var));

    RooExtendPdf  temp(name,name,m_gen_[name],m_real_var_[name]);

    std::cout << "RooContainer::AddGenericPdf -- Made extended PDF " 
	      << name << std::endl;			       
    m_exp_.insert(pair<const char*,RooExtendPdf>(name,temp));
}

void RooContainer::ComposePdf(const char* name, const char * composition
			     ,std::vector<const char*> & formula){

    RooArgList roo_args;
    RooArgList roo_funs;

    for (std::vector<const char*>::iterator it_fun = formula.begin()
	;it_fun != formula.end()
	;it_fun++
	){
	  std::cout << "RooContainer::ComposePdf -- Including Function " 
		    << *it_fun << std::endl;
	  roo_funs.add(m_exp_[*it_fun]);
	  roo_args.add(m_real_var_[(*it_fun)]);
	}

    RooAddPdf temp(name,composition,roo_funs,roo_args);
    //RooAddPdf temp(name,composition,roo_funs);
    std::cout << "RooContainer::ComposePdf -- Created Composed PDF " 
	      << name << std::endl;			       
    m_pdf_.insert(pair<const char*,RooAddPdf>(name,temp));
}


void RooContainer::CreateDataSet(const char *name){
 
    data_[name] = new RooDataSet("data","data",RooArgSet(m_real_var_[name]) );

}

void RooContainer::FitToData(const char* name_func, const char * name_var){

    bool use_composed_pdf = false;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To data " 
	      << name_var
	      << std::endl; 

    RooFitResult *fit_result;
    // Look in the composed pdf before checking the standards
    std::map<const char *,RooAddPdf>::const_iterator it_pdf_ = m_pdf_.find(name_func);
    if (it_pdf_ != m_pdf_.end()){
      fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]));
      use_composed_pdf = true;
    }
    else {
      fit_result = m_exp_[name_func].fitTo(*(data_[name_var]));
    }

    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];
 
    RooPlot *xframe = m_real_var_[name_var].frame(x_min,x_max);
    data_[name_var]->plotOn(xframe);

    if (use_composed_pdf){
      m_pdf_[name_func].plotOn(xframe,LineColor(4));
      int npdfs = m_pdf_[name_func].pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	m_pdf_[name_func].plotOn(xframe
				,Components(*(m_pdf_[name_func].pdfList().at(i)))
				,LineColor(i+1)
				,LineStyle(kDashed));
	m_pdf_[name_func].paramOn(xframe);
      } 
    }
    else {
	m_exp_[name_func].plotOn(xframe,LineColor(4));
	m_exp_[name_func].paramOn(xframe);
    }
    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
    //xframe->Draw();
}

void RooContainer::FitToData(const char* name_func, const char * name_var
			    ,double x1, double x2, double x3, double x4){

    
    float x_min = m_var_min_[name_var];
    float x_max = m_var_max_[name_var];

    bool use_composed_pdf = false;
    std::cout << "RooContainer::FitToData -- Fitting function " 
	      << name_func 
	      << " To data " 
	      << name_var
	      << std::endl; 

    RooFitResult *fit_result;
    // Look in the composed pdf before checking the standards
    std::map<const char *,RooAddPdf>::const_iterator it_pdf_ = m_pdf_.find(name_func);
    if (it_pdf_ != m_pdf_.end()){

      if (x1 < x_min || x4 > x_max){
        std::cout << "RooContainer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		  << std::endl;
        fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]));
        std::cout << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
      } else {
        m_real_var_[name_var].setRange("rnge1",x1,x2);
        m_real_var_[name_var].setRange("rnge2",x3,x4);
        fit_result = m_pdf_[name_func].fitTo(*(data_[name_var]),Range("rnge1,rnge2"));
      }
      use_composed_pdf = true;
    }
    else {
     if (x1 < x_min || x4 > x_max){
        std::cout << "RooContianer::FitToData -- WARNING!! Ranges outside of DataSet Range!" 
		  << std::endl;
        fit_result = m_exp_[name_func].fitTo(*(data_[name_var]));
        std::cout << " Fitted To Full Range -- WARNING!!" 
		  << std::endl;
      } else {
        m_real_var_[name_var].setRange("rnge1",x1,x2);
        m_real_var_[name_var].setRange("rnge2",x3,x4);
        fit_result = m_exp_[name_func].fitTo(*(data_[name_var]),Range("rnge1,rnge2"));
      }
    }

    RooPlot *xframe = m_real_var_[name_var].frame(x_min,x_max);
    data_[name_var]->plotOn(xframe);

    if (use_composed_pdf){
      m_pdf_[name_func].plotOn(xframe,LineColor(4));
      int npdfs = m_pdf_[name_func].pdfList().getSize();
    
      for (int i=0;i<npdfs;i++){
	m_pdf_[name_func].plotOn(xframe
			,Components(*(m_pdf_[name_func].pdfList().at(i)))
			,LineColor(i+1)
			,LineStyle(kDashed));
	m_pdf_[name_func].paramOn(xframe);
      }
    }

    else {
	m_exp_[name_func].plotOn(xframe,LineColor(4));
	m_exp_[name_func].paramOn(xframe);
    }
    fit_res_.insert(std::pair<RooPlot*,RooFitResult*>(xframe,fit_result));
    //xframe->Draw();
}

void RooContainer::SetRealVar(const char * name, float x){

  
  std::map<const char*, RooRealVar>::const_iterator it_var  = m_real_var_.find(name);

    float min_x = m_var_min_[name];
    float max_x = m_var_max_[name];

    if (x > min_x && x < max_x){
      m_real_var_[name] = x;
      data_[name]->add(RooArgSet(m_real_var_[name]));
    }

}

void RooContainer::Save(){
 
  std::map<RooPlot*,RooFitResult*>::const_iterator it;

  std::cout << "RooContainer::Save -- Saving To File "
            << std::endl;

  for(it  = fit_res_.begin()
     ;it != fit_res_.end()
     ;it++ ){
       it->first->Write();
       //it->second->Write();
  }
}

