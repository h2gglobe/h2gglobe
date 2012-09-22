#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
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
  RooRealVar *g0_p0 = new RooRealVar("g0_p0","g0_p0",21294,0.,1e6);
  RooRealVar *g0_p1 = new RooRealVar("g0_p1","g0_p1",0.);//,-0.02,0.02);
  RooRealVar *g0_p2 = new RooRealVar("g0_p2","g0_p2",0.003);//,0.,0.01);
  RooRealVar *g1_p0 = new RooRealVar("g1_p0","g1_p0",500,0.,1e3);
  RooRealVar *g1_p1 = new RooRealVar("g1_p1","g1_p1",0.015,-1.,1.);
  RooRealVar *g1_p2 = new RooRealVar("g1_p2","g1_p2",0.5,0.,1.);
  RooRealVar *g2_p0 = new RooRealVar("g2_p0","g2_p0",250.,0.,1e5);
  RooRealVar *g2_p1 = new RooRealVar("g2_p1","g2_p1",-0.08,-20.,20.);
  RooRealVar *g2_p2 = new RooRealVar("g2_p2","g2_p2",6.,0.,10.);

  RooGaussian *g0 = new RooGaussian("g0","g0",*deltaZ,*g0_p1,*g0_p2);
  RooGaussian *g1 = new RooGaussian("g1","g1",*deltaZ,*g1_p1,*g1_p2);
  RooGaussian *g2 = new RooGaussian("g2","g2",*deltaZ,*g2_p1,*g2_p2);

  RooAddPdf *gausPdf = new RooAddPdf("gausPdf","gausPdf",RooArgList(*g0,*g1,*g2),RooArgList(*g0_p0,*g1_p0,*g2_p0));

  map<string,RooFitResult*> fitResults;

  // loop trees and fit gaussian
  for (map<string,TTree*>::iterator mapIt=trees.begin(); mapIt!=trees.end(); mapIt++){
   
    cout << (mapIt->second)->GetName() << " " << (mapIt->second)->GetEntries() << endl;
    RooDataSet *data = new RooDataSet((mapIt->first).c_str(),(mapIt->first).c_str(),RooArgSet(*deltaZ,*bdtoutput),Import(*(mapIt->second)),Cut("bdtoutput>=0.05"));
    RooFitResult *fitRes = gausPdf->fitTo(*data,Save());
    fitResults.insert(pair<string,RooFitResult*>(mapIt->first,fitRes));

    TCanvas *canv = new TCanvas();
    RooPlot *plot = deltaZ->frame();
    data->plotOn(plot,Binning(200));
    gausPdf->plotOn(plot);
    plot->Draw();
    canv->SaveAs(Form("BSplots/%s.C",(mapIt->first).c_str()));
    canv->Print(Form("BSplots/%s.pdf",(mapIt->first).c_str()));
    canv->Write();
    delete canv;
    delete data;
  }

  for (map<string,RooFitResult*>::iterator iter=fitResults.begin(); iter!=fitResults.end(); iter++){
   
    cout << "--- FIT RESULTS --- " << iter->first << endl;
    (iter->second)->Print();
  }
 
  newFile->Close();
  oldFile->Close();
  outFile->Close();

}
