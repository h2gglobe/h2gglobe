#include "TCanvas.h"
#include "TGraph.h"

#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooSpline1D.h"

#include "../interface/Packager.h"

using namespace std;
using namespace RooFit;

Packager::Packager(RooWorkspace *ws, bool splitVH, int nCats, int mhLow, int mhHigh, bool is2011):
  outWS(ws),
  splitVH_(splitVH),
  nCats_(nCats),
  mhLow_(mhLow),
  mhHigh_(mhHigh),
	is2011_(is2011)
{
  procs.push_back("ggh"); 
  procs.push_back("vbf"); 
  if (splitVH_){
    procs.push_back("wh");
    procs.push_back("zh");
  }
  else {
    procs.push_back("wzh");
  }
  procs.push_back("tth");
	if (is2011) sqrts_=7;
	else sqrts_=8;
  normalization = new Normalization_8TeV(is2011);
}

Packager::~Packager(){}

void Packager::packageOutput(){

  vector<string> expectedObjectsNotFound;

  // sum datasets first
  for (int mh=mhLow_; mh<=mhHigh_; mh+=5){
    RooDataSet *allDataThisMass = 0;
    for (int cat=0; cat<nCats_; cat++) {
      RooDataSet *allDataThisCat = NULL;
      for (vector<string>::iterator proc=procs.begin(); proc!=procs.end(); proc++){
        RooDataSet *tempData = (RooDataSet*)outWS->data(Form("sig_%s_mass_m%d_cat%d",proc->c_str(),mh,cat));
        if (!tempData) {
          cerr << "WARNING -- dataset: " << Form("sig_%s_mass_m%d_cat%d",proc->c_str(),mh,cat) << " not found. It will be skipped" << endl;
          expectedObjectsNotFound.push_back(Form("sig_%s_mass_m%d_cat%d",proc->c_str(),mh,cat));
          continue;
        }
        if (cat==0) allDataThisMass = (RooDataSet*)tempData->Clone(Form("sig_mass_m%d_AllCats",mh));
        else allDataThisMass->append(*tempData);
        if (*proc=="ggh") allDataThisCat = (RooDataSet*)tempData->Clone(Form("sig_mass_m%d_cat%d",mh,cat));
        else allDataThisCat->append(*tempData);
      }
      if (!allDataThisCat) {
        cerr << "WARNING -- allData for cat " << cat << " is NULL. Probably because the relevant datasets couldn't be found. Skipping.. " << endl;
        continue;
      }
      outWS->import(*allDataThisCat);
    }
    if (!allDataThisMass) {
      cerr << "WARNING -- allData for mass " << mh << " is NULL. Probably because the relevant datasets couldn't be found. Skipping.. " << endl;
      continue;
    }
    outWS->import(*allDataThisMass);
  }

  // now create pdf sums (these don't the relative amounts as just used for plotting so can use ThisLum versions)
  RooArgList *sumPdfs = new RooArgList();
  RooArgList *runningNormSum = new RooArgList();
  for (int cat=0; cat<nCats_; cat++){
    RooArgList *sumPdfsThisCat = new RooArgList();
    for (vector<string>::iterator proc=procs.begin(); proc!=procs.end(); proc++){
      
      // sum eA
      RooSpline1D *norm = (RooSpline1D*)outWS->function(Form("hggpdfsmrel_%dTeV_%s_cat%d_norm",sqrts_,proc->c_str(),cat));
      if (!norm) {
        cerr << "WARNING -- ea: " << Form("hggpdfsmrel_%dTeV_%s_cat%d_norm",sqrts_,proc->c_str(),cat) << "not found. It will be skipped" << endl;
      }
      else {
        runningNormSum->add(*norm);
      }
      
      // sum pdf
      RooExtendPdf *tempPdf = (RooExtendPdf*)outWS->pdf(Form("extendhggpdfsmrel_%dTeV_%s_cat%dThisLumi",sqrts_,proc->c_str(),cat));
      if (!tempPdf) {
        cerr << "WARNING -- pdf: " << Form("extendhggpdfsmrel_%dTeV_%s_cat%d",sqrts_,proc->c_str(),cat) << " not found. It will be skipped" << endl;
        expectedObjectsNotFound.push_back(Form("extendhggpdfsmrel_%dTeV_%s_cat%d",sqrts_,proc->c_str(),cat));
        continue;
      }
      sumPdfsThisCat->add(*tempPdf);
      sumPdfs->add(*tempPdf);
    }
    if (sumPdfsThisCat->getSize()==0){
      cerr << "WARNING -- sumPdfs for cat " << cat << " is EMPTY. Probably because the relevant pdfs couldn't be found. Skipping.. " << endl;
      continue;
    }
		// Dont put sqrts here as combine never uses this (but our plotting scripts do)
    RooAddPdf *sumPdfsPerCat = new RooAddPdf(Form("sigpdfrelcat%d_allProcs",cat),Form("sigpdfrelcat%d_allProcs",cat),*sumPdfsThisCat);
    outWS->import(*sumPdfsPerCat,RecycleConflictNodes());
  }
  if (sumPdfs->getSize()==0){
    cerr << "WARNING -- sumAllPdfs is EMPTY. Probably because the relevant pdfs couldn't be found. Skipping.. " << endl;
  }
  else {
		// Dont put sqrts here as combine never uses this (but our plotting scripts do)
    RooAddPdf *sumPdfsAllCats = new RooAddPdf("sigpdfrelAllCats_allProcs","sigpdfrelAllCats_allProcs",*sumPdfs);
    outWS->import(*sumPdfsAllCats,RecycleConflictNodes());
  }

  if (runningNormSum->getSize()==0){
    cerr << "WARNING -- runningNormSum is EMPTY. Probably because the relevant normalizations couldn't be found. Skipping.. " << endl;
  }
  else {
    RooAddition *normSum = new RooAddition("normSum","normSum",*runningNormSum);
    outWS->import(*normSum,RecycleConflictNodes());
    
    RooRealVar *MH = (RooRealVar*)outWS->var("MH");
    RooRealVar *intLumi = (RooRealVar*)outWS->var("IntLumi");
    RooAddition *norm = (RooAddition*)outWS->function("normSum");
    TGraph *effAccGraph = new TGraph();
    TGraph *expEventsGraph = new TGraph();
    int p=0;
    for (double mh=mhLow_; mh<mhHigh_+0.5; mh+=1){
      MH->setVal(mh);
      expEventsGraph->SetPoint(p,mh,intLumi->getVal()*norm->getVal());
      effAccGraph->SetPoint(p,mh,norm->getVal()/(normalization->GetXsection(mh)*normalization->GetBR(mh)));
      p++;
    }
    TCanvas *canv = new TCanvas();
    effAccGraph->SetLineWidth(3);
    effAccGraph->GetXaxis()->SetTitle("m_{H} (GeV)");
    effAccGraph->GetYaxis()->SetTitle("efficiency #times acceptance");
    effAccGraph->Draw("AL");
    canv->Print("plots/effAccCheck.pdf");
    canv->Print("plots/effAccCheck.png");
    expEventsGraph->SetLineWidth(3);
    expEventsGraph->GetXaxis()->SetTitle("m_{H} (GeV)");
    expEventsGraph->GetYaxis()->SetTitle(Form("Expected Events for %4.1ffb^{-1}",intLumi->getVal()/1000.));
    expEventsGraph->Draw("AL");
    canv->Print("plots/expEventsCheck.pdf");
    canv->Print("plots/expEventsCheck.png");

  }
}
