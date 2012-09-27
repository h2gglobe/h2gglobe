#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TSpline.h"

#include <iostream>

void cdf(TH1F* hcdf, TH1F* h) {
  for (int i=0; i<h->GetNbinsX(); i++)
    hcdf->SetBinContent(i+1, h->Integral(1, i+1));
}

float findBin(float x, TH1F* h) {

  int nBins = h->GetNbinsX();
  for (int i=0; i<nBins; i++) {
    if (h->GetBinContent(i+1) >= x){
      return (i+1);
    }
  }
  
  return -1;
}

Float_t scalingFunction(Float_t x, TH1F* h1, TH1F* h2, TRandom3* r) {
  
  Int_t iBin1 = h1->FindBin(x);
  Double_t x1[3], y1[3];
  x1[1] = h1->GetBinCenter(iBin1);
  y1[1] = h1->GetBinContent(iBin1);
  x1[2] = h1->GetBinCenter(iBin1+1);
  y1[2] = h1->GetBinContent(iBin1+1);
  x1[0] = h1->GetBinCenter(iBin1-1);
  y1[0] = h1->GetBinContent(iBin1-1);
  TSpline3* spline1 = new TSpline3("spline1", x1, y1, 3);
  Float_t integral = spline1->Eval(x);
  //h1->Draw();
  //h1->SetLineColor(kRed);
  //spline1->Draw("SAME");
  //spline1->SetLineColor(kRed);
  //h2->Draw("SAME");
  //h2->SetLineColor(kGreen);
  //std::cout << iBin << std::endl;
  //Float_t integral = h1->GetBinContent(iBin);
  //std::cout << integral << std::endl;
  //Float_t newval = findBin(integral, h2);// * h2->GetBinWidth(1);
  Double_t x2[3], y2[3];
  int iBin2 = findBin(integral, h2)-1;
  //std::cout << iBin2 << " " << x << std::endl;
  x2[1] = h2->GetBinCenter(iBin2);
  y2[1] = h2->GetBinContent(iBin2);
  x2[2] = h2->GetBinCenter(iBin2+1);
  y2[2] = h2->GetBinContent(iBin2+1);
  x2[0] = h2->GetBinCenter(iBin2-1);
  y2[0] = h2->GetBinContent(iBin2-1);
  //TSpline3* spline2 = new TSpline3("spline2", x2, y2, 3);
  TSpline3* spline3 = new TSpline3("spline3", y2, x2, 3);
  //std::cout << integral << " " << spline3->Eval(integral) << std::endl;
  //spline3->Draw("SAME");
  //spline3->SetLineColor(kGreen);
  Float_t newval = spline3->Eval(integral);
  //std::cout << x << " " << newval <<  " " << iBin2 << std::endl;
  //newval = newval;// + r->Uniform(-h2->GetBinWidth(1), h2->GetBinWidth(1));
  //std::cout << x << " " << newval * h2->GetBinWidth(1) << std::endl;
  return newval;//* h2->GetBinWidth(1);
}


void scaling() {
  
  //ROOT.gStyle.SetOptStat(0)

  TRandom3* random = new TRandom3();

  TFile* f = new TFile("out.root");
  TH1F* data = (TH1F*)f->Get("data");
  TH1F* mc = (TH1F*)f->Get("mc");

  data->Scale(1./data->Integral());
  mc->Scale(1./mc->Integral());

  TH1F* hd, *hm;
  char a[100];
  sprintf(a, "%s_cdf", data->GetName());
  hd = (TH1F*)data->Clone(a);
  sprintf(a, "%s_cdf", mc->GetName());
  hm = (TH1F*)mc->Clone(a);

  cdf(hd, data);
  cdf(hm, mc);

  //hd.SetLineColor(ROOT.kRed)
  //hd.Draw()
  //hm.Draw("SAME")

  Float_t itype, vtxprob;
  
  TFile* fo = new TFile("opttree.root");
  TTree* t = (TTree*)fo->Get("opttree");
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("vtxprob", 1);
  t->SetBranchStatus("itype", 1);
  t->SetBranchAddress("vtxprob", &vtxprob);
  t->SetBranchAddress("itype", &itype);

  TH1F* scaled = new TH1F("scaled", "scaled", 1000, 0, 1);

  Int_t entries = t->GetEntries();
  for (int z=1600000; z<entries; z++) {
    if (z%10000 == 1)
      std::cout << z << std::endl;
    t->GetEntry(z);

    if (itype == 41) {
      Float_t newval = scalingFunction(vtxprob, hm, hd, random);
      scaled->Fill(newval);
    }
  }

  scaled->Scale(1./scaled->Integral());

  scaled->Rebin(2);
  //scaled->Smooth(4);
  scaled->Draw();
  scaled->SetLineColor(kGreen);
  data->SetLineColor(kBlue);
  data->Rebin(2);
  data->Draw("SAME");
  data->SetLineColor(kRed);
  mc->Rebin(2);
  mc->Draw("SAME");
}


//TSpline3()
//findBin(integral, h2)
