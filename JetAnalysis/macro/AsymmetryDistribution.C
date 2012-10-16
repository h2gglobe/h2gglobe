#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

#include "TROOT.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TDirectory.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooHistPdf.h"
#include "RooStats/RooStatsUtils.h"

#define ASYM_UP 1.0
#define ASYM_DOWN -1.0
#define EXPECTED_EVENTS 48

using namespace RooFit;

TCut base  = "diphotonMVA>-0.05 && category>=4 && (sampleType<0 || (mass>130)) && sampleType>-400";
TCut tight = "leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>30 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>500 && TMath::Abs(deltaPhiJJGamGam)>2.6";
TCut loose = "leadPtOverM>0.5 && subleadPtOverM>0.3 && leadJPt>30 && subleadJPt>20 && TMath::Abs(deltaEtaJJ)>3 && TMath::Abs(Zep)<2.5 && MJJ>250 && TMath::Abs(deltaPhiJJGamGam)>2.6";
TCut cat4  = "diphotonMVA>-0.05 && category==4";
TCut cat5  = "diphotonMVA>-0.05 && category==5";

typedef Int_t Yield_t;

struct AsymHyp
{
  TString name;
  TH1F* dist;
  std::map<Yield_t, std::vector<std::vector<Double_t> > > toys;
  std::map<Yield_t, std::vector<Double_t> > asym;
};

class AsymTester
{
  public:
    AsymTester();
    ~AsymTester();

    void AddHyp(TString name, TH1F* dist);
    void AddYield(Yield_t y, TString name);

    void SetToys(Int_t t);

    void GenIt();
    TCanvas* PlotYield(Yield_t);
    void PrintLLR(Yield_t);
    void SetThreshold(Double_t t);

  private:
    std::vector<std::pair<Yield_t, TString> > yields_;
    std::vector<AsymHyp> hypothesis_;

    Int_t nToys_;
    Double_t threshold_;

    void GenToys();
    void GenAsym();
    void OptimizeThreshold();

    std::vector<std::vector<Double_t> > GenToys(UInt_t hyp, UInt_t yield);
    std::vector<Double_t> GenAsym(UInt_t hyp, UInt_t yield);

    TH1F* GetHist(UInt_t hyp, UInt_t yield);

    Double_t GetMedian(UInt_t hyp, UInt_t yield);
    Double_t GetLLR(Yield_t id, char op = ' ');

  protected:
};

template <class T>
TString my_toa(T i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

AsymTester::AsymTester():
  nToys_(2000)
{
}

AsymTester::~AsymTester()
{
  for(UInt_t i = 0; i < hypothesis_.size(); i++)
  {
    delete hypothesis_[i].dist;
    hypothesis_[i].dist = NULL;
  }
}

void AsymTester::AddHyp(TString name, TH1F* dist)
{
  AsymHyp temp;
  temp.name = name;
  temp.dist = (TH1F*)dist->Clone();

  hypothesis_.push_back(temp);
}

void AsymTester::AddYield(Yield_t y, TString name)
{
  std::pair<Int_t, TString> temp;
  temp.first = y;
  temp.second = name;
  yields_.push_back(temp);
}

void AsymTester::SetToys(Int_t t)
{
  nToys_ = t;
}

void AsymTester::OptimizeThreshold()
{
  threshold_ = TMath::Pi()/2;

  if(hypothesis_.size() == 0)
  {
    threshold_ = TMath::Pi()/2;
    return;
  }

  Int_t nBins = hypothesis_[0].dist->GetNbinsX();
  Double_t Min = hypothesis_[0].dist->GetXaxis()->GetXmin();
  Double_t Max = hypothesis_[0].dist->GetXaxis()->GetXmax();
  Double_t Step = (Max-Min)/nBins;

  Double_t maxThreshold = 0.0, maxLLR = 1.0;

  for(Int_t i = nBins/3; i <= 2*nBins/3; i++)
  {
    threshold_ = Min + i*Step;
    GenAsym();
    Double_t LLR = GetLLR(0);
    if(LLR < maxLLR)
    {
      maxLLR = LLR;
      maxThreshold = threshold_;
    }
  }

  threshold_ = maxThreshold;

  std::cout << "Max LLR p-value: " << maxLLR << " -> Z-val: " << RooStats::PValueToSignificance(maxLLR) << std::endl << "threshold: " << threshold_ << std::endl;
}

void AsymTester::GenIt()
{
  GenToys();
  OptimizeThreshold();
  GenAsym();
}

void AsymTester::GenToys()
{
  for(UInt_t i = 0; i < hypothesis_.size(); i++)
  {
    for(UInt_t j = 0; j < yields_.size(); j++)
    {
      std::vector<std::vector<Double_t> > temp = GenToys(i, j);
      hypothesis_[i].toys[yields_[j].first] = temp;
    }
  }
}

std::vector<std::vector<Double_t> > AsymTester::GenToys(UInt_t hyp, UInt_t yield)
{
  TAxis* axis = hypothesis_[hyp].dist->GetXaxis();
  RooRealVar x("x", "x", axis->GetXmin(), axis->GetXmax());
  RooDataHist hist("hist", "hist", x, hypothesis_[hyp].dist);
  RooHistPdf pdf("pdf", "pdf", x, hist);

  std::vector<Double_t> toy;
  std::vector<std::vector<Double_t> > toys;

  for(Int_t i = 0; i < nToys_; i++)
  {
    RooDataSet* toyDS = pdf.generate(x, yields_[yield].first, AutoBinned(false));

    for(Int_t j = 0; j < yields_[yield].first; j++)
    {
      const RooArgSet* ev = toyDS->get(j);
      toy.push_back(ev->getRealValue(x.GetName()));
      //std::cout << toy[j] << std::endl;
      for(Int_t k = j; k > 0; k--)
      {
        if(toy[k] < toy[k-1])
          std::swap(toy[k], toy[k-1]);
        else
          break;
      }
    }
    toys.push_back(toy);
    toy.clear();
    toy.resize(0);

    delete toyDS;
  }
  return toys;
}

void AsymTester::GenAsym()
{
  for(UInt_t i = 0; i < hypothesis_.size(); i++)
  {
    for(UInt_t j = 0; j < yields_.size(); j++)
    {
      std::vector<Double_t> temp = GenAsym(i, j);
      hypothesis_[i].asym[yields_[j].first] = temp;
    }
  }
}

std::vector<Double_t> AsymTester::GenAsym(UInt_t hyp, UInt_t yield)
{
  std::vector<Double_t> asym;

  for(UInt_t i = 0; i < hypothesis_[hyp].toys[yields_[yield].first].size(); i++)
  {
    Int_t Low = 0;
    Int_t High = 0;
    for(UInt_t j = 0; j < hypothesis_[hyp].toys[yields_[yield].first][i].size(); j++)
    {
      if(hypothesis_[hyp].toys[yields_[yield].first][i][j] < threshold_)
        Low++;
      else
        High++;
    }

    asym.push_back(((Double_t)(High-Low))/((Double_t)(Low+High)));

    for(UInt_t j = i; j > 0; j--)
    {
      if(asym[j] < asym[j-1])
      {
        std::swap(asym[j], asym[j-1]);
        std::swap(hypothesis_[hyp].toys[yields_[yield].first][j], hypothesis_[hyp].toys[yields_[yield].first][j-1]);
      }
      else
        break;
    }
  }

  return asym;
}

TCanvas* AsymTester::PlotYield(Yield_t yield)
{
  Bool_t found = false;
  UInt_t id = 0;

  for(UInt_t i = 0; i < yields_.size(); i++)
  {
    if(yields_[i].first == yield)
    {
      found = true;
      id = i;
      break;
    }
  }

  if(!found)
    return NULL;

  TString name = "Yield_";
  name = name + my_toa(yield);
  name = name + "-";
  name = name + yields_[id].second;
  TCanvas* canvas = new TCanvas(name, name, 800, 600);
  Double_t MaxVal = 0.0;
  TH1F* first;

  for(UInt_t i = 0; i < hypothesis_.size(); i++)
  {
    TH1F* temp = GetHist(i, id);
    temp->SetBit(TH1::kNoTitle);
    temp->SetBit(TH1::kNoStats);

    if(temp->GetMaximum() > MaxVal)
      MaxVal = temp->GetMaximum();

    if(i==0)
    {
      temp->Draw();
      first = temp;
    }
    else
    {
      switch((i-1)%6)
      {
        case 0:
          temp->SetLineColor(kPink + (i-1)/6);
          break;
        case 1:
          temp->SetLineColor(kAzure + (i-1)/6);
          break;
        case 2:
          temp->SetLineColor(kSpring + (i-1)/6);
          break;
        case 3:
          temp->SetLineColor(kOrange + (i-1)/6);
          break;
        case 4:
          temp->SetLineColor(kViolet + (i-1)/6);
          break;
        case 5:
          temp->SetLineColor(kTeal + (i-1)/6);
          break;
      }
      temp->Draw("same");
    }
  }

  first->GetYaxis()->SetRangeUser(0, MaxVal*1.05);
  canvas->BuildLegend(0.65, 0.75, 0.88, 0.88);

  return canvas;
}

TH1F* AsymTester::GetHist(UInt_t hyp, UInt_t yield)
{
  Yield_t nEvents = yields_[yield].first;
  Double_t Overshoot = (ASYM_UP-ASYM_DOWN)/(2.*(Double_t)((Int_t)nEvents));
  TH1F* AsymmetryDistrib = new TH1F(hypothesis_[hyp].name+"_AsymmetryDistribution_"+yields_[yield].second, hypothesis_[hyp].name, nEvents+1, ASYM_DOWN-Overshoot, ASYM_UP+Overshoot);

  for(UInt_t i = 0; i < hypothesis_[hyp].asym[yields_[yield].first].size(); i++)
  {
    AsymmetryDistrib->Fill(hypothesis_[hyp].asym[yields_[yield].first][i]);
  }

  AsymmetryDistrib->Scale(1/AsymmetryDistrib->Integral());

  return AsymmetryDistrib;
}

void AsymTester::PrintLLR(Yield_t yield)
{
  Bool_t found = false;
  UInt_t id = 0;

  for(UInt_t i = 0; i < yields_.size(); i++)
  {
    if(yields_[i].first == yield)
    {
      found = true;
      id = i;
      break;
    }
  }

  if(!found)
  {
    std::cout << "There is no " << yield << " yield.";
    return;
  }

  if(hypothesis_.size() < 2)
  {
    std::cout << "At least 2 hypothesis are needed!";
    return;
  }

  std::cout << std::endl << "Printing LLR for a yield of " << yield << " (" << yields_[id].second << ")" << std::endl;
  std::cout << "Only H0 and H1 are used, other hypothesis are ignored." << std::endl;
  std::cout << "\tH0 -> " << hypothesis_[0].name << std::endl << "\tH1 -> " << hypothesis_[1].name << std::endl;

  Double_t LLR = GetLLR(id, 'p');

  std::cout << "\tLLR: " << LLR << std::endl;
}

Double_t AsymTester::GetLLR(Yield_t id, char op)
{
  if(hypothesis_.size() < 2)
  {
    return 0.;
  }

  TCanvas* temp = new TCanvas("temp", "temp", 200, 200);
  TH1F* H0 = GetHist(0, id);
  TH1F* H1 = GetHist(1, id);

  std::vector<Double_t> LLRH0;
  std::vector<Double_t> LLRH1;

  for(UInt_t i = 0; i < hypothesis_[0].asym[yields_[id].first].size(); i++)
  {
    Double_t asym = hypothesis_[0].asym[yields_[id].first][i];
    Double_t LLR = (H1->GetBinContent(H1->FindBin(asym)))/(H0->GetBinContent(H0->FindBin(asym)));
    LLR = TMath::Log(LLR);

    LLRH0.push_back(LLR);

    for(Int_t j = i; j > 0; j--)
    {
      if(LLRH0[j] < LLRH0[j-1])
        std::swap(LLRH0[j], LLRH0[j-1]);
      else
        break;
    }
  }
  for(UInt_t i = 0; i < hypothesis_[1].asym[yields_[id].first].size(); i++)
  {
    Double_t asym = hypothesis_[1].asym[yields_[id].first][i];
    Double_t LLR = (H0->GetBinContent(H0->FindBin(asym)))/(H1->GetBinContent(H1->FindBin(asym)));
    LLR = -TMath::Log(LLR);

    LLRH1.push_back(LLR);

    for(Int_t j = i; j > 0; j--)
    {
      if(LLRH1[j] < LLRH1[j-1])
        std::swap(LLRH1[j], LLRH1[j-1]);
      else
        break;
    }
  }

  delete H0;
  delete H1;
  delete temp;

  if(op == 'p')
  {
    TString name = "LLR";
    name = name+yields_[id].second;
    TCanvas* canvas = new TCanvas(name, name, 800, 600);

    Double_t min = 0.0;
    UInt_t size0 = LLRH0.size();
    UInt_t size1 = LLRH1.size();
    Double_t max = 0.0;
    UInt_t nBins = 10;

    Double_t Inf = -TMath::Log(0);

    for(UInt_t i = 0; i < size0; i++)
    {
      if(LLRH0[i] != Inf)
        if(LLRH0[i] > max)
          max = LLRH0[i];
      if(LLRH0[i] != -Inf)
        if(LLRH0[i] < min)
          min = LLRH0[i];
    }
    for(UInt_t i = 0; i < size1; i++)
    {
      if(LLRH1[i] != Inf)
        if(LLRH1[i] > max)
          max = LLRH1[i];
      if(LLRH1[i] != -Inf)
        if(LLRH1[i] < min)
          min = LLRH1[i];
    }

    min = min - 0.01;
    max = max + 0.01;

    TH1F* LLR0 = new TH1F(hypothesis_[0].name+yields_[id].second, hypothesis_[0].name, nBins, min, max);
    TH1F* LLR1 = new TH1F(hypothesis_[1].name+yields_[id].second, hypothesis_[1].name, nBins, min, max);

    for(UInt_t i = 0; i < size0; i++)
      LLR0->Fill(LLRH0[i]);
    for(UInt_t i = 0; i < size1; i++)
      LLR1->Fill(LLRH1[i]);

    LLR0->SetBit(TH1::kNoTitle);
    LLR0->SetBit(TH1::kNoStats);
    LLR0->Draw();
    LLR1->SetLineColor(kRed);
    LLR1->SetBit(TH1::kNoTitle);
    LLR1->SetBit(TH1::kNoStats);
    LLR1->Draw("same");
  }

  Double_t MedianH0 = 0.0;
  Double_t MedianH1 = 0.0;

  UInt_t size = LLRH0.size();
  if(size%2)
  {
    MedianH0 = LLRH0[size/2];
  }
  else
  {
    MedianH0 = (LLRH0[size/2]+LLRH0[(size/2)-1])/(2.);
  }
  size = LLRH1.size();
  if(size%2)
  {
    MedianH1 = LLRH1[size/2];
  }
  else
  {
    MedianH1 = (LLRH1[size/2]+LLRH1[(size/2)-1])/(2.);
  }

  UInt_t count = 0;
  for(UInt_t i = 0; i < LLRH1.size(); i++)
  {
    if(LLRH1[i] < MedianH0)
    {
      count++;
    }
    if(LLRH1[i] == MedianH0)
    {
      if(count < LLRH1.size()/2)
        count++;
    }
    if(LLRH1[i] > MedianH0)
    {
      break;
    }
  }

  if(count >= LLRH1.size()/2)
    return ((Double_t)(LLRH1.size()-count))/((Double_t)LLRH1.size());
  return ((Double_t)count)/((Double_t)LLRH1.size());
}

Double_t AsymTester::GetMedian(UInt_t hyp, UInt_t yield)
{
  UInt_t size = hypothesis_[hyp].asym[yields_[yield].first].size();
  if(size%2)
  {
    return hypothesis_[hyp].asym[yields_[yield].first][size/2];
  }
  else
  {
    Double_t a = hypothesis_[hyp].asym[yields_[yield].first][size/2-1];
    Double_t b = hypothesis_[hyp].asym[yields_[yield].first][size/2];
    return (a+b)/2;
  }
}

void AsymTester::SetThreshold(Double_t t)
{
  threshold_ = t;
  GenAsym();
}

/*std::vector<std::vector<Double_t> > GenToys(TH1F* Distrib, Int_t nEvents, Int_t nToys = 2000)
{
  TAxis* axis = Distrib->GetXaxis();
  RooRealVar x("x", "x", axis->GetXmin(), axis->GetXmax());
  RooDataHist hist("hist", "hist", x, Distrib);
  RooHistPdf pdf("pdf", "pdf", x, hist);

  std::vector<Double_t> toy;
  std::vector<std::vector<Double_t> > toys;

  for(Int_t i = 0; i < nToys; i++)
  {
    RooDataSet* toyDS = pdf.generate(x, nEvents, AutoBinned(false));

    for(Int_t j = 0; j < nEvents; j++)
    {
      const RooArgSet* ev = toyDS->get(j);
      toy.push_back(ev->getRealValue(x.GetName()));
      //std::cout << toy[j] << std::endl;
      for(Int_t k = j; k > 0; k--)
      {
        if(toy[k] < toy[k-1])
          std::swap(toy[k], toy[k-1]);
        else
          break;
      }
    }
    toys.push_back(toy);
    toy.clear();
    toy.resize(0);

    delete toyDS;
  }
  return toys;
}

std::vector<Double_t> GenAsym(std::vector<std::vector<Double_t> > toys, Double_t threshold)
{
  std::vector<Double_t> asym;

  for(Int_t i = 0; i < toys.size(); i++)
  {
    Int_t Low = 0;
    Int_t High = 0;
    for(Int_t j = 0; j < toys[i].size(); j++)
    {
      if(toys[i][j] < threshold)
        Low++;
      else
        High++;
    }

    asym.push_back(((Double_t)(High-Low))/((Double_t)(Low+High)));

    for(Int_t j = i; j > 0; j--)
    {
      if(asym[j] < asym[j-1])
        std::swap(asym[j], asym[j-1]);
      else
        break;
    }
  }

  return asym;
}

TH1F* GenAsymDistrib(std::vector<std::vector<Double_t> > toys, Double_t threshold)
{
  std::vector<Double_t> asym = GenAsym(toys, threshold);
  Int_t nEvents = toys[0].size();

  Double_t Overshoot = (ASYM_UP-ASYM_DOWN)/(2.*(Double_t)((Int_t)nEvents));
  TH1F* AsymmetryDistrib = new TH1F("Asymmetry Distribution", "Asymmetry Distribution", nEvents+1, ASYM_DOWN-Overshoot, ASYM_UP+Overshoot);

  for(Int_t i = 0; i < asym.size(); i++)
  {
    AsymmetryDistrib->Fill(asym[i]);
  }

  return AsymmetryDistrib;
}

TH1F* GenAsymDistrib(TH1F* angularDistrib, Int_t nEvents, TString name, Int_t nToys = 2000)
{
  RooRealVar x("x", "x", 0, 2*TMath::Pi());
  RooDataHist hist("hist", "hist", x, angularDistrib);

  RooHistPdf pdf("pdf", "pdf", x, hist);

  TDirectory* workDir = gDirectory->mkdir(name);
  TDirectory* currentDir = gDirectory;
  workDir->cd();

  Double_t* Asymmetry = new Double_t[nToys];
  Double_t Overshoot = (ASYM_UP-ASYM_DOWN)/(2.*(Double_t)((Int_t)nEvents));
  TH1F* AsymmetryDistrib = new TH1F("Asymmetry Distribution", "Asymmetry Distribution", nEvents+1, ASYM_DOWN-Overshoot, ASYM_UP+Overshoot);

  for(Int_t i = 0; i < nToys; i++)
  {
    RooDataSet* toy = pdf.generate(x, nEvents);

    TH1F* temp = new TH1F(my_toa(i), "toy distribution", 40, 0, TMath::Pi());
    toy->fillHistogram(temp, x);
    temp->Write();
    temp->Rebin(20);

    Double_t nPlus = temp->GetBinContent(2), nMinus = temp->GetBinContent(1);
    Asymmetry[i] = (nPlus-nMinus)/(nPlus+nMinus);

    AsymmetryDistrib->Fill(Asymmetry[i]);
    //std::cout << Asymmetry[i] << std::endl;

    for(Int_t j = i; j > 0; j--)
    {
      if(Asymmetry[j] < Asymmetry[j-1])
        std::swap(Asymmetry[j], Asymmetry[j-1]);
      else
        break;
    }

    delete temp;
    delete toy;
  }

  AsymmetryDistrib->Write();
  AsymmetryDistrib->Draw();
  AsymmetryDistrib->SetNameTitle(name, name);
  AsymmetryDistrib->Scale(1/AsymmetryDistrib->Integral());
  //Output: Min:Median:Max
  std::cout << name << ": " << Asymmetry[0] << ":" << ((nToys%2)?(Asymmetry[nToys/2]):((Asymmetry[nToys/2]+Asymmetry[nToys/2-1])/2.)) << ":" << Asymmetry[nToys-1] << std::endl;

  delete [] Asymmetry;
  currentDir->cd();
  return AsymmetryDistrib;
}*/

void AsymmetryDistribution()
{
  TString baseDir = "/afs/cern.ch/user/c/cbeiraod/working/Properties/scripts";
  TString toDraw = "abs(deltaPhiJJ)";
  TCut    category = base && tight;

  TFile* VBF_f  = new TFile(baseDir+"/VBFNLO4.root",    "READ");
  TFile* Back_f = new TFile(baseDir+"/Background.root", "READ");
  TFile* Output_f = new TFile(baseDir+"/Asymmetry.root", "RECREATE");
  Output_f->cd();

  TTree* VBF_s0_t     = (TTree*)VBF_f->Get("vbfnlo_m125_8TeV"       );
  TTree* VBF_s2_t     = (TTree*)VBF_f->Get("vbfnlo_spin2_m125_8TeV" );
  //TTree* VBF_cpeven_t = (TTree*)VBF_f->Get("vbfnlo_cpeven_vbf_m125_8TeV");

  std::vector<TTree*> Backgrounds;
  Backgrounds.push_back((TTree*)Back_f->Get("diphojet_8TeV"  ));
  Backgrounds.push_back((TTree*)Back_f->Get("gjet_20_8TeV_pf"));
  Backgrounds.push_back((TTree*)Back_f->Get("gjet_40_8TeV_pf"));

  TH1F** histos = new TH1F*[Backgrounds.size()];
  Float_t* weights = new Float_t[Backgrounds.size()];
  TCanvas* canvas = new TCanvas("Canvas", "Canvas", 800, 600);

  TH1F* Background = new TH1F("Background", "Background"  , 40, 0, TMath::Pi());
  TH1F* VBF_s0     = new TH1F("VBF_s0"    , "VBF_{s0}"    , 40, 0, TMath::Pi());
  TH1F* VBF_s2     = new TH1F("VBF_s2"    , "VBF_{s2}"    , 40, 0, TMath::Pi());
  //TH1F* VBF_cpeven = new TH1F("VBF_cpeven", "VBF_{cpeven}", 40, 0, TMath::Pi());

  for(Int_t i = 0; i < (Int_t)Backgrounds.size(); i++)
  {
    TString histName = Backgrounds[i]->GetName();
    histName += "_h";//my_toa(i);

    histos[i] = new TH1F(histName, histName, 40, 0, TMath::Pi());

    Backgrounds[i]->Draw(toDraw + ">>" + histName, category, "hist");

    if((Int_t)(Backgrounds[i]->GetEntries()) == 0)
    {
      weights[i] = 0;
    }
    else
    {
      Backgrounds[i]->SetBranchAddress("weight", &(weights[i]));
      Backgrounds[i]->GetEntry(0);

      //std::cout << weights[i] << std::endl;
    }

    Background->Add(histos[i], (Double_t)weights[i]);
  }

  VBF_s0_t->Draw(toDraw + ">>VBF_s0", category, "hist");
  VBF_s2_t->Draw(toDraw + ">>VBF_s2", category, "hist");
  //VBF_cpeven_t->Draw(toDraw + ">>VBF_cpeven", category, "hist");

  Background->Write();
  for(Int_t i = 0; i < (Int_t)Backgrounds.size(); i++)
    histos[i]->Write();
  VBF_s0->Write();
  VBF_s2->Write();
  //VBF_cpeven->Write();

  AsymTester asymmetryMeasurement;
  asymmetryMeasurement.AddYield(12, "GammaGamma");
  asymmetryMeasurement.AddYield(24, "GammaGamma+WW");
  asymmetryMeasurement.AddYield(48, "GammaGamma+WW+tautau");
  asymmetryMeasurement.AddHyp("Spin0", VBF_s0);
  asymmetryMeasurement.AddHyp("Spin2", VBF_s2);
  //asymmetryMeasurement.AddHyp("cpeven", VBF_cpeven);

  asymmetryMeasurement.GenIt();
  asymmetryMeasurement.SetThreshold(TMath::Pi()/2);

  asymmetryMeasurement.PlotYield(12);
  asymmetryMeasurement.PlotYield(24);
  asymmetryMeasurement.PlotYield(48);

  asymmetryMeasurement.PrintLLR(12);
  asymmetryMeasurement.PrintLLR(24);
  asymmetryMeasurement.PrintLLR(48);

  /*TH1F* Asym_s0     = GenAsymDistrib(GenToys(VBF_s0, 12, 2000), TMath::Pi()/2);
  TH1F* Asym_s2     = GenAsymDistrib(GenToys(VBF_s2, 12, 2000), TMath::Pi()/2);
  TH1F* Asym_cpeven = GenAsymDistrib(GenToys(VBF_cpeven, 12, 2000), TMath::Pi()/2);
  //TH1F* Asym_s0     = GenAsymDistrib(VBF_s0    , EXPECTED_EVENTS, "Spin0");
  //TH1F* Asym_s2     = GenAsymDistrib(VBF_s2    , EXPECTED_EVENTS, "Spin2");
  //TH1F* Asym_cpeven = GenAsymDistrib(VBF_cpeven, EXPECTED_EVENTS, "CPeven");

  Asym_s0->SetBit(TH1::kNoTitle);
  Asym_s0->SetBit(TH1::kNoStats);
  Asym_s2->SetBit(TH1::kNoTitle);
  Asym_s2->SetBit(TH1::kNoStats);
  Asym_cpeven->SetBit(TH1::kNoTitle);
  Asym_cpeven->SetBit(TH1::kNoStats);
  Asym_s2->SetLineColor(kRed);
  Asym_cpeven->SetLineColor(kGreen);

  Asym_s0->Draw();
  Asym_s2->Draw("same");
  //Asym_cpeven->Draw("same");

  Double_t MaxVal = Asym_s0->GetMaximum();
  if(Asym_s2->GetMaximum() > MaxVal)
    MaxVal = Asym_s2->GetMaximum();
  if(Asym_cpeven->GetMaximum() > MaxVal)
    MaxVal = Asym_cpeven->GetMaximum();
  Asym_s0->GetYaxis()->SetRangeUser(0, MaxVal*1.05);
  canvas->BuildLegend(0.65, 0.75, 0.88, 0.88);

  Asym_s0->Write();
  Asym_s2->Write();
  Asym_cpeven->Write();*/

  //Background->Draw();
  //VBF_cpeven_t->Draw("TMath::Cos(thetaJ1)", base, "norm hist");

  //Output_f->Close();
  Back_f->Close();
  VBF_f->Close();

  return;
}
