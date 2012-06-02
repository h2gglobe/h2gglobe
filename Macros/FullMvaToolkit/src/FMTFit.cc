#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"

#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"

#include "../interface/FMTFit.h"

using namespace std;
using namespace RooFit;

FMTFit::FMTFit(TFile *tFile, double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int,vector<double> > > AllBinEdges, bool verbose):
	
	FMTBase(intLumi,is2011,mHMinimum, mHMaximum, mHStep, massMin, massMax, nDataBins, signalRegionWidth, sidebandWidth, numberOfSidebands, numberOfSidebandsForAlgos, numberOfSidebandGaps, massSidebandMin, massSidebandMax, includeVBF, nVBFCategories, includeLEP, nLEPCategories, systematics, rederiveOptimizedBinEdges, AllBinEdges, verbose),
	blind_(false),
	plot_(true)
{
  gROOT->SetStyle("Plain");
  r1 = new RooRealVar("r1","r1",-8.,-50.,0.); 
  r2 = new RooRealVar("r2","r2",-1.,-50.,0.); 
  f1 = new RooRealVar("f1","f1",0.5,0.,1.); 
  nBkgInSigReg = new RooRealVar("nbis","nbis",10,0,100000);
	outWS = (RooWorkspace*)tFile->Get("cms_hgg_workspace");
  mass_var = (RooRealVar*)outWS->var("CMS_hgg_mass");
	
	// get data and combine all cats
	cout << "Looking for datasets....." << endl;
	data = (RooDataSet*)((RooDataSet*)outWS->data("data_mass_cat0"))->Clone("data_mass");
	if (getincludeVBF()) data->append(*((RooDataSet*)outWS->data("data_mass_cat1")));
	if (getincludeLEP()) data->append(*((RooDataSet*)outWS->data("data_mass_cat2")));
	if (!outWS->data("data_mass")) outWS->import(*data);

}

FMTFit::~FMTFit(){
	delete r1;
	delete r2;
	delete f1;
	delete nBkgInSigReg;
}

pair<double,double> FMTFit::FitPow(double mass){
  
  gROOT->SetBatch();
  if (mass<getmHMinimum() || mass>getmHMaximum()) {
    cout << Form("%3.1f invalid mass. Not in range [",mass) << getmHMinimum() << "," << getmHMaximum() << "]" << endl;
    exit(1);
  }

	// set up fit function
	r1->SetName(Form("r1_%3.1f",mass));
  r2->SetName(Form("r2_%3.1f",mass));
  f1->SetName(Form("f1_%3.1f",mass));
  fit  = new RooGenericPdf(Form("data_pow_model_%3.1f",mass),"data_pow_model","(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)",RooArgList(*mass_var,*r1,*r2,*f1));

	// set up fit region
	double mLow = getmassMin();
	double mHigh = getmassMax();
  double sidebL = mass*(1-getsignalRegionWidth());
  double sidebH = mass*(1+getsignalRegionWidth());
  mass_var->setRange(Form("rangeLow_m%3.1f",mass),mLow,sidebL);
  mass_var->setRange(Form("rangeHig_m%3.1f",mass),sidebH,mHigh);
  mass_var->setRange(Form("sigReg_m%3.1f",mass),sidebL,sidebH);

	data->Print();
	cout << data->GetName() << " " << data->numEntries() << endl;
	
  RooFitResult *fitRes = fit->fitTo(*data,Range(Form("rangeLow_m%3.1f,rangeHig_m%3.1f",mass,mass)),Save(true),Strategy(1));

	// make plot
	if (plot_) Plot(mass);
  
	// integral in sig region
  RooAbsReal *integral = fit->createIntegral(*mass_var,NormSet(*mass_var),Range(Form("sigReg_m%3.1f",mass)));
  // comb integral in two sideband regions 
  RooAbsReal *sidebandInt = fit->createIntegral(*mass_var,NormSet(*mass_var),Range(Form("rangeLow_m%3.1f,rangeHig_m%3.1f",mass,mass)));
  double tempEv = data->sumEntries(Form("CMS_hgg_mass>=%3.1f && CMS_hgg_mass<%3.1f",mLow,sidebL))+data->sumEntries(Form("CMS_hgg_mass>%3.1f && CMS_hgg_mass<=%3.1f",sidebH,mHigh));
  RooConstVar sidebandNevents(Form("sbEvents_m%3.1f",mass),Form("sbEvents_m%3.1f",mass),tempEv);
  RooFormulaVar normIntVar(Form("normIntVar_m%3.1f",mass),Form("normIntVar_m%3.1f",mass),"@0*@1/@2",RooArgSet(sidebandNevents,*integral,*sidebandInt));

  double result     = normIntVar.getVal();
  double fullError  = normIntVar.getPropagatedError(*fitRes);
  
	// if fit isn't in workspace then import else it will have been automatically updated (just make sure the workspace is written somewhere)
  // bit messy but gets the job done

  nBkgInSigReg = outWS->var(Form("NBkgInSignal_mH%3.1f",mass));
  if (nBkgInSigReg!=NULL){
    nBkgInSigReg = (RooRealVar*)outWS->var(Form("NBkgInSignal_mH%3.1f",mass));
    nBkgInSigReg->setVal(result);
    nBkgInSigReg->setError(fullError);
		cout << "---- var " << nBkgInSigReg->GetName() << " already found - updating ------" << endl;
    outWS->Write(outWS->GetName(),TObject::kWriteDelete);
  }
  else {
    RooRealVar *temp = new RooRealVar(Form("NBkgInSignal_mH%3.1f",mass),"t",10,0,1e6);
    temp->setVal(result);
    temp->setError(fullError);
    outWS->import(*temp);
		cout << "---- writing var " << temp->GetName() << " to workspace -----" << endl;
    outWS->Write(outWS->GetName(),TObject::kWriteDelete);
    delete temp;
  }
	if (verbose_) outWS->allVars().Print();
	delete fit;

  return pair<double,double>(result,fullError);
}

void FMTFit::Plot(double mass){
    
    RooPlot *frame = mass_var->frame(Title(Form("Mass fit for %3.1f",mass)));
    frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    mass_var->setRange("unblindReg",150,getmassMax());
		if (blind_) {
      data->plotOn(frame,Binning(getnDataBins()),CutRange("unblindReg"));
      data->plotOn(frame,Binning(getnDataBins()),Invisible());
		}
    else data->plotOn(frame,Binning(getnDataBins()));
    fit->plotOn(frame,Range(Form("rangeLow_m%3.1f,rangeHig_m%3.1f",mass,mass)),NormRange(Form("rangeLow_m%3.1f,rangeHig_m%3.1f",mass,mass)));
    frame->SetMinimum(0.0001);
    TCanvas *c1 = new TCanvas();
		frame->Draw();
		// make signal reg and sideband boxes
		TLine l1, l2;
    l1.SetLineColor(kRed);
    l1.SetLineWidth(2);
    l2.SetLineColor(kBlue-2);
    l2.SetLineWidth(2);
    l2.SetLineStyle(9);
    TBox b1, b2;
    b1.SetFillColor(kRed-9);
    b2.SetFillColor(kBlue-9);
    b1.SetFillStyle(3003);
    b2.SetFillStyle(3003);
    b1.IsTransparent();
    b2.IsTransparent();
		double sidebL = mass*(1-getsignalRegionWidth());
		double sidebH = mass*(1+getsignalRegionWidth());
    l1.DrawLine(sidebL,frame->GetMinimum(),sidebL,frame->GetMaximum());
    l1.DrawLine(sidebH,frame->GetMinimum(),sidebH,frame->GetMaximum());
    b1.DrawBox(sidebL,frame->GetMinimum(),sidebH,frame->GetMaximum());
    vector<double> lowEdges = getLowerSidebandEdges(mass);
    vector<double> highEdges = getUpperSidebandEdges(mass);
    for (int i=0; i<lowEdges.size(); i++) {
      l2.DrawLine(lowEdges[i],frame->GetMinimum(),lowEdges[i],frame->GetMaximum());
      if (i>0) b2.DrawBox(lowEdges[i-1],frame->GetMinimum(),lowEdges[i],frame->GetMaximum());
    }
    for (int i=0; i<highEdges.size(); i++) {
      l2.DrawLine(highEdges[i],frame->GetMinimum(),highEdges[i],frame->GetMaximum()); 
      if (i>0) b2.DrawBox(highEdges[i-1],frame->GetMinimum(),highEdges[i],frame->GetMaximum());
    }
		TLatex *text = new TLatex();
		text->SetTextSize(0.04);
		text->SetNDC();
		text->DrawLatex(0.68,0.85,"CMS preliminary");
    text->DrawLatex(0.75,0.78,"#sqrt{s} = 8 TeV");
		text->DrawLatex(0.73,0.71,Form("#int L = %1.1f fb^{-1}",getintLumi()));
    if (blind_) text->DrawLatex(0.67,0.64,"Blinded: [100,150]");
    c1->SaveAs(Form("plots/pdf/fit_m%3.1f.pdf",mass));
    c1->SaveAs(Form("plots/png/fit_m%3.1f.png",mass));
}

void FMTFit::redoFit(double mass){
  
  pair<double,double> dummyVar = FitPow(mass);
}

bool FMTFit::getblind(){
	return blind_;
}

bool FMTFit::getplot(){
	return plot_;
}

void FMTFit::setblind(bool blind){
	blind_=blind;
}

void FMTFit::setplot(bool plot){
	plot_=plot;
}



