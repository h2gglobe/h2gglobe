// Original Author - Doug Berry

#include "TH1F.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooMsgService.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "th1fmorph.C"
#include "Normalization_8TeV.h"
//#include "Normalization_ff.C"

using namespace std;
using namespace RooFit;

Normalization_8TeV *normalizer = new Normalization_8TeV();

template <class type> string makestring(type value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

void dofit(double fitmass, vector <TString> InterpolationList, TFile* SourceFile, TFile* OutputFile, RooWorkspace* WorkSpace, int noverwritemasses=0, double *overwritemasses=NULL, int debug=0) {

  if (floor(fitmass)-fitmass<0.0000001 && floor(fitmass)-fitmass>0) fitmass=floor(fitmass);
  if (fitmass-ceil(fitmass)>-0.0000001 && fitmass-ceil(fitmass)<0) fitmass=ceil(fitmass);
  
  if (fitmass>150 || fitmass<110) {
    cout << "Warning!!!!!!!!!!! You must have an input mass between 110 and 150 GeV!" << endl << "Skipping !!!!" << endl;
    return;
  }

  bool OverWriteMass = false;
  if (overwritemasses!=NULL) {
    for (int i=0; i<noverwritemasses; i++) {
      if (fabs(overwritemasses[i]-fitmass)<0.0000001) OverWriteMass=true;
    }
  }
  
  double Masses[] = {105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0}; //Mass points used for interpolation
  size_t NumMasses = sizeof(Masses)/sizeof(Masses[0]);
  double lowerbound = 0;
  double upperbound = 0;
  for (size_t i=0; i<NumMasses; i++) {
    if (fitmass>Masses[i] && fitmass<Masses[i+1]) {
      lowerbound = Masses[i];
      upperbound = Masses[i+1];
    } else if (fitmass==Masses[0] || fitmass==Masses[NumMasses-1]) {
      lowerbound = Masses[i];
      upperbound = Masses[i];
    } else if (fitmass==Masses[i]) {
      lowerbound = Masses[i-1];
      upperbound = Masses[i+1];
    }
  }

  TString MassString = makestring(fitmass);
  //MassString.ReplaceAll(".","_");
  TString LowerBoundString = makestring(lowerbound);
  LowerBoundString.ReplaceAll(".0","");
  TString UpperBoundString = makestring(upperbound);
  UpperBoundString.ReplaceAll(".0","");
  RooRealVar RooRealMass = *(WorkSpace->var("CMS_hgg_mass"));
  
  cout << "Calculating mass point at " << fitmass << endl;
  
  for (unsigned int k=0; k < InterpolationList.size(); k++) {

    TString LowerHistName = InterpolationList[k];
    LowerHistName.ReplaceAll("110",LowerBoundString);
    TString UpperHistName = InterpolationList[k];
    UpperHistName.ReplaceAll("110",UpperBoundString);
    TString HistName = InterpolationList[k];
    HistName.ReplaceAll("110",MassString);
    
    TString HistTitle(Form("Interpolated Mass at %fGeV",fitmass));

    TH1F* LowerHist = (TH1F*) SourceFile->Get(LowerHistName.Data());
    TH1F* UpperHist = (TH1F*) SourceFile->Get(UpperHistName.Data());
    if (debug>=1) cout << "Calculating mass point at " << fitmass << "GeV with histograms " << LowerHistName << " and " << UpperHistName << endl;

    double Normalization = normalizer->GetNorm(lowerbound, LowerHist, upperbound, UpperHist, fitmass);

    TH1F* MCHist = (TH1F*) SourceFile->Get(HistName.Data());
    TH1F* InterpolatedHist = (TH1F*) th1fmorph((Char_t*) HistName.Data(),(Char_t*) HistTitle.Data(),LowerHist,UpperHist,lowerbound,upperbound,fitmass,Normalization,0);

    if (MCHist!=NULL && !OverWriteMass) {
      TString ResidualHistName = HistName;
      ResidualHistName += "_Residual";
      TH1F* ResidualHist = (TH1F*) InterpolatedHist->Clone(ResidualHistName.Data());
      ResidualHist->Add(MCHist,-1);
      OutputFile->WriteTObject(ResidualHist);
      ResidualHistName.ReplaceAll("th1f_","");
      RooDataHist RooDataResidual(Form("roohist_%s",ResidualHistName.Data()),ResidualHistName.Data(),RooRealMass,ResidualHist);
      WorkSpace->import(RooDataResidual);
    }

    if (MCHist==NULL || OverWriteMass) {
      OutputFile->WriteTObject(InterpolatedHist,InterpolatedHist->GetName());
      HistName.ReplaceAll("th1f_","");
      RooDataHist RooDataInterpolated(Form("roohist_%s",HistName.Data()),HistName.Data(),RooRealMass,InterpolatedHist);
      if (WorkSpace->data(Form("roohist_%s",HistName.Data()))==NULL) WorkSpace->import(RooDataInterpolated);
      else {
        if (debug>=1) cout << "Before: " << ((RooDataHist*) WorkSpace->data(Form("roohist_%s",HistName.Data())))->sumEntries() << endl;
        ((RooDataHist*) WorkSpace->data(Form("roohist_%s",HistName.Data())))->reset();
        ((RooDataHist*) WorkSpace->data(Form("roohist_%s",HistName.Data())))->add(RooDataInterpolated);
        if (debug>=1) cout << "After: " << ((RooDataHist*) WorkSpace->data(Form("roohist_%s",HistName.Data())))->sumEntries() << endl;
      }
    } else {
      HistName += "_Interpolated";
      InterpolatedHist->SetName(HistName.Data());
      OutputFile->WriteTObject(InterpolatedHist);
      HistName.ReplaceAll("th1f_","");
      RooDataHist RooDataInterpolated(Form("roohist_%s",HistName.Data()),HistName.Data(),RooRealMass,InterpolatedHist);
      WorkSpace->import(RooDataInterpolated);
      
    }
  }
  cout << "Done " << fitmass << endl;

}

void InterpolateMass(double fitmass, TString SourceFileName="CMS-HGG.root", double overwritemass=0) {

  TString FileNameMod = Form("_%s.root",makestring(fitmass).c_str());
  TString FileName = SourceFileName;
  FileName.ReplaceAll(".root",FileNameMod);

  if (fitmass>=150 || fitmass<=110) {
    cout << "Warning!!!!!!!!!!! You must have an input mass between 110 and 150 GeV!" << endl << "Exiting Program!!!!" << endl;
    exit(1);
  }

  cout << "Opening file: " << SourceFileName << endl;
  TFile* SourceFile = new TFile(SourceFileName.Data());
  TList* HistList = SourceFile->GetListOfKeys();
  RooWorkspace * WorkSpace = (RooWorkspace*) SourceFile->Get("cms_hgg_workspace");
  TFile* OutputFile = new TFile(FileName.Data(),"RECREATE");
  OutputFile->cd();
  vector<TString> InterpolationList;
  for (Int_t j=0; j<HistList->GetSize(); ++j) {
    TString HistName(HistList->At(j)->GetName());
    if (HistName.Contains("110_")) InterpolationList.push_back(HistName);
    if (HistName.Contains("th1f")) {
      TH1F* temphist = (TH1F*) SourceFile->Get(HistName.Data());
      OutputFile->WriteTObject(temphist);
    }
    if (HistName.Contains("plot_data_pol_model")) {
      TCanvas* tempcan = (TCanvas*) SourceFile->Get(HistName.Data());
      OutputFile->WriteTObject(tempcan);
    }
  }
  dofit(fitmass, InterpolationList, SourceFile, OutputFile, WorkSpace, 1, &overwritemass);

  cout << "Writing data to disk..." << endl;
  WorkSpace->Write();
  OutputFile->Close();
  delete WorkSpace;
  delete SourceFile;
  delete OutputFile;
  cout << "Done!" << endl;
}

void InterpolateMassPoints(int nmasses, double * masses, TString SourceFileName="CMS-HGG",
                           TString FileName="", int noverwritemasses=0, double *overwritemasses=NULL) 
{
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if( FileName == "" ) { FileName = SourceFileName + "_interpolated.root"; }

  SourceFileName += ".root";
  TFile* SourceFile = new TFile(SourceFileName);
  TList* HistList = SourceFile->GetListOfKeys();
  RooWorkspace * WorkSpace = (RooWorkspace*) SourceFile->Get("cms_hgg_workspace");
  TFile* OutputFile = new TFile(FileName.Data(),"RECREATE");
  OutputFile->cd();
  vector<TString> InterpolationList;

  for (Int_t j=0; j<HistList->GetSize(); ++j) {

    TString HistName(HistList->At(j)->GetName());
    if (HistName.Contains("110_")) InterpolationList.push_back(HistName);
    if (HistName.Contains("th1f")) {
      TH1F* temphist = (TH1F*) SourceFile->Get(HistName.Data());
      TString temphistname = temphist->GetName();
      OutputFile->WriteTObject(temphist);
    }
    if (HistName.Contains("plot_data_pol_model")) {
      TCanvas* tempcan = (TCanvas*) SourceFile->Get(HistName.Data());
      OutputFile->WriteTObject(tempcan);
    }
    
  }
  
  for(int imass=0; imass<nmasses; ++imass) {
    dofit(masses[imass], InterpolationList, SourceFile, OutputFile, WorkSpace, noverwritemasses, overwritemasses);
  }
  
  cout << "Writing data to disk..." << endl;
  WorkSpace->Write();
  OutputFile->Close();
  delete WorkSpace;
  delete SourceFile;
  delete OutputFile;
  cout << "Done!" << endl;

}

void InterpolateMassRange(double Min, double Max, double Step, TString SourceFileName="CMS-HGG", int noverwritemasses=0, double *overwritemasses=NULL) 
{
  TString FileName = "";
  std::vector<double> masses;
  for (double fitmass=Min; fitmass<=Max; fitmass+=Step) {
    masses.push_back(fitmass);
  }

  InterpolateMassPoints(masses.size(), &masses[0], SourceFileName,FileName, noverwritemasses, overwritemasses); 
  
}


#ifndef __CINT__

int main(int argc, char ** argv)
{
	InterpolateMassRange(atof(argv[1]),atof(argv[2]),atof(argv[3]));
}

#endif
