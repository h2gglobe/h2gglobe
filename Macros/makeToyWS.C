#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooGenericPdf.h"

using namespace std;
using namespace RooFit;

void printTree(TTree *tree){
  
  float bdtoutput;
  float mgg;
  float weight;
  int vbf;

  tree->SetBranchAddress("bdtoutput",&bdtoutput);
  tree->SetBranchAddress("CMG_hgg_mass",&mgg);
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("vbf",&vbf);

  for (int e=0; e<tree->GetEntries(); e++){
    tree->GetEntry(e);
    cout << Form("M: %3.1f -- B: %1.3f -- W: %1.3f -- V: %d",mgg,bdtoutput,weight,vbf) << endl;
  }

}

void makeToyWS(string inFileName, string outFileName){

  gROOT->SetBatch();

  TFile *outFile = new TFile(outFileName.c_str(),"RECREATE");
  RooWorkspace *outWS = new RooWorkspace("fits_workspace");

  TFile *treeFile = TFile::Open(inFileName.c_str());
  TTree *dataTree = (TTree*)treeFile->Get("dataTree");
  TTree *sigTree = (TTree*)treeFile->Get("sigTree");
  //printTree(sigTree);
  //printTree(dataTree);

  RooRealVar *mass = new RooRealVar("CMS_hgg_mass","CMS_hgg_mass",100,180);
  RooRealVar *bdtoutput = new RooRealVar("bdtoutput","bdtoutput",-1,1);
  RooRealVar *vbf = new RooRealVar("vbf","vbf",0,1);
  RooRealVar *weight = new RooRealVar("weight","weight",0,1);

  RooDataSet *data_novbf = new RooDataSet("data_bdt_novbf","data_bdt",RooArgSet(*bdtoutput,*mass,*vbf),Import(*dataTree),Cut("vbf==0"));
  RooDataSet *data_vbf = new RooDataSet("data_bdt_vbf","data_bdt",RooArgSet(*bdtoutput,*mass,*vbf),Import(*dataTree),Cut("vbf==1"));

  RooKeysPdf *data_pdf = new RooKeysPdf("data_pdf","data_pdf",*bdtoutput,*data_novbf);

  RooDataSet *data_formass_all = new RooDataSet("data_bdt_cut_all","data_bdt",RooArgSet(*bdtoutput,*mass,*vbf),Import(*dataTree),Cut("bdtoutput>=0.05"));
  RooRealVar *r1 = new RooRealVar("r1","r1",-8.,-50.,0.); 
  RooRealVar *r2 = new RooRealVar("r2","r2",-1.,-50.,0.); 
  RooRealVar *f1 = new RooRealVar("f1","f1",0.5,0.,1.); 
  RooGenericPdf *fit  = new RooGenericPdf("data_pow_model","data_pow_model","(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)",RooArgList(*mass,*r1,*r2,*f1));
  fit->fitTo(*data_formass_all);

  TCanvas *canv = new TCanvas();
  RooPlot *datB = bdtoutput->frame();
  data_novbf->plotOn(datB,Binning(100));
  data_pdf->plotOn(datB);
  datB->Draw();
  canv->SaveAs("datB.pdf");

  RooPlot *datM = mass->frame();
  data_novbf->plotOn(datM,Binning(160));
  datM->Draw();
  canv->SaveAs("datM.pdf");

  RooDataSet *sig_novbf = new RooDataSet("sig_bdt_novbf","sig_bdt",RooArgSet(*bdtoutput,*mass,*vbf,*weight),Import(*sigTree),Cut("vbf==0"),WeightVar(*weight));
  RooDataSet *sig_vbf = new RooDataSet("sig_bdt_vbf","sig_bdt",RooArgSet(*bdtoutput,*mass,*vbf,*weight),Import(*sigTree),Cut("vbf==1"),WeightVar(*weight));
  
  bdtoutput->setBins(50);
  mass->setBins(360);
  
  RooDataHist *sig_hist = new RooDataHist("sig_bdt_hist","sig_bdt",RooArgSet(*bdtoutput,*mass),*sig_novbf);
  RooHistPdf *sig_pdf = new RooHistPdf("sig_pdf","sig_pdf",RooArgSet(*bdtoutput,*mass),*sig_hist);

  RooPlot *sigB = bdtoutput->frame();
  sig_novbf->plotOn(sigB,Binning(50));
  sig_pdf->plotOn(sigB);
  sigB->Draw();
  canv->SaveAs("sigB.pdf");

  RooPlot *sigM = mass->frame();
  sig_novbf->plotOn(sigM,Binning(160));
  sig_pdf->plotOn(sigM);
  sigM->Draw();
  canv->SaveAs("sigM.pdf");
  

  outWS->import(*data_novbf);
  outWS->import(*data_vbf);
  outWS->import(*data_pdf);
  outWS->import(*data_formass_all);
  outWS->import(*fit);
  outWS->import(*sig_novbf);
  outWS->import(*sig_vbf);
  outWS->import(*sig_pdf);

  outFile->cd();
  outWS->Write();
  sigTree->Write();
  dataTree->Write();

  outFile->Close();
  treeFile->Close();



}
