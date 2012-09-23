#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooAbsArg.h"
#include "RooArgList.h"
#include "RooPlot.h"

using namespace std;
using namespace RooFit;

void BeamspotParametrization(string filename_newBS, string filename_oldBS){
 
  gROOT->SetBatch();

  TFile *newFile = TFile::Open(filename_newBS.c_str());
  TFile *oldFile = TFile::Open(filename_oldBS.c_str());
  TFile *outFile = new TFile("BSplots/bs_fits.root","RECREATE");
 
  // BS variable
  RooRealVar *deltaZ = new RooRealVar("deltaZgenToChosen","deltaZgenToChosen",-25.,25.);
  RooRealVar *bdtoutput = new RooRealVar("bdtoutput","bdtoutput",-1,1.);
  
  // get trees
  map<string,TTree*> trees;
  trees.insert(pair<string,TTree*>("ggh_m124_newBS",(TTree*)newFile->Get("ggh_m124_8TeV_newBS")));
  //trees.insert(pair<string,TTree*>("ggh_m125_newBS",(TTree*)newFile->Get("ggh_m125_8TeV_newBS")));
  //trees.insert(pair<string,TTree*>("ggh_m126_newBS",(TTree*)newFile->Get("ggh_m126_8TeV_newBS")));
  //trees.insert(pair<string,TTree*>("ggh_m124_oldBS",(TTree*)oldFile->Get("ggh_m124_8TeV")));
  //trees.insert(pair<string,TTree*>("ggh_m125_oldBS",(TTree*)oldFile->Get("ggh_m125_8TeV")));
  //trees.insert(pair<string,TTree*>("ggh_m123_oldBS",(TTree*)oldFile->Get("ggh_m123_8TeV")));

  // set up triple gaussian
  // norm terms
  RooRealVar *norm0 = new RooRealVar("n0","n0",2e4,0,1e6);
  RooRealVar *norm1 = new RooRealVar("n1","n1",500,0,1e6);
  RooRealVar *norm2 = new RooRealVar("n2","n2",500,0,1e6);

  // mean terms
  RooRealVar *mean0 = new RooRealVar("m0","m0",0.,-0.02,0.02);
  RooRealVar *mean1 = new RooRealVar("m1","m1",0.,-1.,1.);
  RooRealVar *mean2 = new RooRealVar("m2","m2",0.,-10.,10.);

  // sigma terms
  RooRealVar *sigma0 = new RooRealVar("s0","s0",0.003,0.,0.02);
  RooRealVar *sigma1 = new RooRealVar("s1","s1",0.5,0.,5.);
  RooRealVar *sigma2 = new RooRealVar("s2","s2",6.,0.,10.);

  RooGaussian *gaus0 = new RooGaussian("g0","g0",*deltaZ,*mean0,*sigma0);
  RooGaussian *gaus1 = new RooGaussian("g1","g1",*deltaZ,*mean1,*sigma1);
  RooGaussian *gaus2 = new RooGaussian("g2","g2",*deltaZ,*mean2,*sigma2);
  
  RooExtendPdf *ext0 = new RooExtendPdf("e0","e0",*gaus0,*norm0);
  RooExtendPdf *ext1 = new RooExtendPdf("e1","e1",*gaus1,*norm1);
  RooExtendPdf *ext2 = new RooExtendPdf("e2","e2",*gaus2,*norm2);
  
  RooAddPdf *gausPdf = new RooAddPdf("gausPdf","gausPdf",RooArgList(*ext0,*ext1,*ext2));

  map<string,RooFitResult*> fitResults;

  deltaZ->setRange("narrow",-0.02,0.02);
  deltaZ->setRange("wide_low",-25.,-2.);
  deltaZ->setRange("wide_high",2.,25.);
  deltaZ->setRange("mid_low",-2.,-0.02);
  deltaZ->setRange("mid_high",0.02,2.);

  TCanvas *canv = new TCanvas();
  // loop trees and fit gaussian
  for (map<string,TTree*>::iterator mapIt=trees.begin(); mapIt!=trees.end(); mapIt++){
   
    cout << (mapIt->second)->GetName() << " " << (mapIt->second)->GetEntries() << endl;
    RooDataSet *data = new RooDataSet((mapIt->first).c_str(),(mapIt->first).c_str(),RooArgSet(*deltaZ,*bdtoutput),Import(*(mapIt->second)),Cut("bdtoutput>=0.05"));
    // first fit narrow gaussian
    ext0->fitTo(*data,Range("narrow"),Save());
    mean0->setRange(mean0->getVal(),mean0->getVal());
    sigma0->setRange(sigma0->getVal(),sigma0->getVal());
    // then fit mid gaussian
    ext1->fitTo(*data,Range("mid_low,mid_high"),Save());
    mean1->setRange(mean1->getVal(),mean1->getVal());
    sigma1->setRange(sigma1->getVal(),sigma1->getVal());
    // then fit wide gaussian
    ext2->fitTo(*data,Range("wide_low,wide_high"),Save());
    mean2->setRange(mean2->getVal(),mean2->getVal());
    sigma2->setRange(sigma2->getVal(),sigma2->getVal());
    
    RooPlot *plotN = deltaZ->frame();
    data->plotOn(plotN,Binning(200));
    ext0->plotOn(plotN,Range("narrow"),NormRange("narrow"));
    plotN->Draw();
    plotN->GetXaxis()->SetRangeUser(-0.05,0.05);
    canv->Print(Form("BSplots/nar_%s.pdf",(mapIt->first).c_str()));
    
    RooPlot *plotM = deltaZ->frame();
    data->plotOn(plotM,Binning(200));
    ext1->plotOn(plotM,Range("mid_low,mid_high"),NormRange("mid_low,mid_high"));
    plotM->Draw();
    plotM->GetYaxis()->SetRangeUser(0,500);
    canv->Print(Form("BSplots/mid_%s.pdf",(mapIt->first).c_str()));
    
    RooPlot *plotW = deltaZ->frame();
    data->plotOn(plotW,Binning(200));
    ext2->plotOn(plotW,Range("wide_low,wide_high"),NormRange("wide_low,wide_high"));
    plotW->Draw();
    plotW->GetYaxis()->SetRangeUser(0,500);
    canv->Print(Form("BSplots/wide_%s.pdf",(mapIt->first).c_str()));
   
    // then fit all
    RooFitResult *fitRes = gausPdf->fitTo(*data,Save());
    
    fitResults.insert(pair<string,RooFitResult*>(mapIt->first,fitRes));

    RooPlot *plot = deltaZ->frame();
    data->plotOn(plotW,Binning(200));
    gausPdf->plotOn(plot);
    plot->GetYaxis()->SetRangeUser(0,500);
    plot->Draw();
    canv->SaveAs(Form("BSplots/%s.C",(mapIt->first).c_str()));
    canv->Print(Form("BSplots/%s.pdf",(mapIt->first).c_str()));
    canv->Write();
    
    system(Form("cp BSplots/%s.pdf ~/www/",(mapIt->first).c_str()));
    system(Form("cp BSplots/nar_%s.pdf ~/www/",(mapIt->first).c_str()));
    system(Form("cp BSplots/mid_%s.pdf ~/www/",(mapIt->first).c_str()));
    system(Form("cp BSplots/wide_%s.pdf ~/www/",(mapIt->first).c_str()));
  }

  for (map<string,RooFitResult*>::iterator iter=fitResults.begin(); iter!=fitResults.end(); iter++){
   
    cout << "--- FIT RESULTS --- " << iter->first << endl;
    (iter->second)->Print();
  }
 
  newFile->Close();
  oldFile->Close();
  outFile->Close();

}
