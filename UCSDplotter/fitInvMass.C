#include "TROOT.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TProcessID.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"

#include <vector>
#include <iostream>


std::vector<TCanvas*> canvases;
std::vector<std::vector<RooRealVar*> > vars;
std::vector<RooRealVar*> norms;
std::vector<RooVoigtian*> voigts;
std::vector<RooExponential*> exps;
std::vector<RooAddPdf*> pdfs;
std::vector<RooExtendPdf*> ePdfs;

void saveFile(TFile* f, const char* OutFile, Int_t nCat, RooWorkspace* rw) {

  TFile* out = new TFile(OutFile, "recreate");
  f->cd();
  
  TIter nextTopLevelKey;
  nextTopLevelKey=(f->GetListOfKeys());
  TKey *keyTopLevel;
  while (keyTopLevel = (TKey*)nextTopLevelKey()) {
    
    TString name(keyTopLevel->GetName());
    //std::cout << name << std::endl;
    TString className(keyTopLevel->GetClassName());
    //std::cout << className << std::endl;
    
    if (className.CompareTo("TH1F") == 0) {
      TH1F* th1 = (TH1F*) f->Get(name);
      out->cd();
      th1->Write();
      f->cd();
    }
  
    if (className.CompareTo("TMacro") == 0) {
      TMacro* macro = (TMacro*) f->Get(name);
      out->cd();
      macro->Write();
      f->cd();
    }
  
    //if (className.CompareTo("TProcessID") == 0) {
    //  TProcessID* pid = (TProcessID*) f->Get(name);
    //  out->cd();
    //  pid->Write();
    //  f->cd();
    //}
  }
  
  for (int i=0; i<nCat; i++) {
    for (unsigned int j=0; j<vars[i].size(); j++) 
      rw->import(*(vars[i][j]));
    
    rw->import(*(norms[i]));
    rw->import(*(voigts[i]));
    rw->import(*(exps[i]));
    rw->import(*(pdfs[i]), RooFit::RecycleConflictNodes());
    rw->import(*(ePdfs[i]), RooFit::RecycleConflictNodes());
  }

  out->cd();
  rw->Write();
  for (unsigned int i=0; i<canvases.size(); i++) 
    canvases[i]->Write();
  out->Close();
}

int fitInvMass(const char* InputFile, const char* OutFile, Int_t nCat=4) {
  
  TFile* f = new TFile(InputFile);
  RooWorkspace *w = (RooWorkspace*)(f->Get("cms_hgg_workspace"));
  
  char a[100];
  
  for (int i=0; i<nCat; i++) {
    std::vector<RooRealVar*> temp;
    sprintf(a, "CMS_hll_voigtexp_0_8TeV_cat%d", i);
    temp.push_back(new RooRealVar(a, a, 91.186, 80, 100));
    sprintf(a, "CMS_hll_voigtexp_1_8TeV_cat%d", i);
    temp.push_back(new RooRealVar(a, a, 2.495, 0., 10.));
    sprintf(a, "CMS_hll_voigtexp_2_8TeV_cat%d", i);
    temp.push_back(new RooRealVar(a, a, 5., 0., 10.));
    sprintf(a, "CMS_hll_voigtexp_3_8TeV_cat%d", i);
    temp.push_back(new RooRealVar(a, a, -.001, -1., -1e-6));
    sprintf(a, "CMS_hll_voigtexp_4_8TeV_cat%d", i);
    temp.push_back(new RooRealVar(a, a, .5, .0, 1.0));
    vars.push_back(temp);
    
    sprintf(a, "pdf_data_voigtexp_model_8TeV_cat%d_norm", i);
    norms.push_back(new RooRealVar(a, a, 1e3, 0, 1e8));

    sprintf(a, "voigt_cat%d", i);
    voigts.push_back(new RooVoigtian(a, a, *(w->var("CMS_hll_mass")), *vars[i][0], *vars[i][1], *vars[i][2]));
    sprintf(a, "exp_cat%d", i);
    exps.push_back(new RooExponential(a, a, *(w->var("CMS_hll_mass")), *vars[i][3]));
    sprintf(a, "pdf_data_voigtexp_model_8TeV_cat%d", i);
    pdfs.push_back(new RooAddPdf(a, a, RooArgList(*(voigts[i]), *(exps[i])), *vars[i][4], true));
    sprintf(a, "data_voigtexp_model_8TeV_cat%d", i);
    ePdfs.push_back(new RooExtendPdf(a, a, *(pdfs[i]), *norms[i]));
  }

  for(int i=0; i<nCat; i++) {
    sprintf(a, "data_mass_cat%d", i);
    RooDataSet* ds = (RooDataSet*)w->data(a);
    RooDataSet* pippo = (RooDataSet*)ds->reduce("CMS_hll_mass < 120 || CMS_hll_mass > 130");
    RooPlot* frame = (RooPlot*)w->var("CMS_hll_mass")->frame();//110, 180);
    w->var("CMS_hll_mass")->setRange("range1", 88, 94);
    voigts[i]->fitTo(*ds, RooFit::Range("range1"), RooFit::Extended(false), RooFit::Save(true), RooFit::Strategy(1), RooFit::NumCPU(8));
    w->var("CMS_hll_mass")->setRange("range2", 130, 180);
    exps[i]->fitTo(*ds, RooFit::Range("range2"), RooFit::Extended(false), RooFit::Save(true), RooFit::Strategy(1));

    //vars[i][3]->setConstant(1);
    w->var("CMS_hll_mass")->setRange("range3", 110, 120);
    w->var("CMS_hll_mass")->setRange("range4", 130, 180);
    ePdfs[i]->setNormRange("range3,range4");
    ePdfs[i]->fitTo(*ds, RooFit::Range("range3,range4"), RooFit::Extended(true), RooFit::Save(true), RooFit::Strategy(1), RooFit::NumCPU(8));
    //vars[i][3]->setConstant(0);

    
    //////bool blind_data = true;
    //////const RooCmdArg blind = (blind_data) ? RooFit::Invisible() : RooFit::MarkerColor(kBlack);
    pippo->plotOn(frame);//, RooFit::Invisible());
    voigts[i]->plotOn(frame, RooFit::LineStyle(kDashed));
    exps[i]->plotOn(frame, RooFit::LineStyle(kDashed));
    ePdfs[i]->plotOn(frame);
    ePdfs[i]->paramOn(frame, RooFit::Layout(0.25,0.95,0.86));

    //////pdf_saves_.push_back(&(it_pdf->second));
    //////n_pars=(it_pdf->second).getParameters(it_data->second)->getSize() -1 ;
    //////(it_exp->second).plotOn(xframe,LineColor(4),RooFit::Range("rnge1,rnge2"),RooFit::Normalization(m_real_var_[name_func].getVal(),RooAbsReal::NumEvent));
    //////fit_result->printValue(std::cout);

    gROOT->SetBatch(true);
    gROOT->SetStyle("Plain");
    sprintf(a,"plot_%s_%s", ePdfs[i]->GetName(), ds->GetName());
    canvases.push_back(new TCanvas(a, a, 1200, 900));    
    canvases[i]->cd(); 
    frame->Draw();
  }

  saveFile(f, OutFile, nCat, w);
  return 0;
}
