#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TMath.h"
#include "Math/DistFunc.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/NumberCountingUtils.h"
#include "Math/QuantFuncMathMore.h"

using namespace std;

namespace separation
{
  #define JK_PARAM 20

  class  NormalQuantile
  {
  public:
	  NormalQuantile() {};
	  double operator() (double *x, double *p) {
		  return ROOT::Math::normal_quantile_c(x[0], p[0]);
	  };
  };

  template <class T>
  std::pair<std::pair<double, double>, std::pair<double, double> > p_value(std::vector<T> H0, std::vector<T> H1)
  {
    for(UInt_t i = 0; i < H0.size()-2; i++)
      for(UInt_t j = 0; j < H0.size()-1-i; j++)
        if(H0[j] > H0[j+1])
          std::swap(H0[j], H0[j+1]);
    for(UInt_t i = 0; i < H1.size()-2; i++)
      for(UInt_t j = 0; j < H1.size()-1-i; j++)
        if(H1[j] > H1[j+1])
          std::swap(H1[j], H1[j+1]);

    UInt_t size = H0.size();
    T H0med = (size%2)?(H0[size/2]):((H0[size/2]+H0[size/2-1])/2);
    size = H1.size();
    T H1med = (size%2)?(H1[size/2]):((H1[size/2]+H1[size/2-1])/2);

    double H0pval = 0;
    double H1pval = 0;
    for(UInt_t i = 0; i < H0.size(); i++)
      if(H0[i] > H1med)
      {
        H0pval = H0.size()-(i-1);
        break;
      }
    for(UInt_t i = 0; i < H1.size(); i++)
      if(H1[i] > H0med)
      {
        H1pval = i;
        break;
      }

    separation::NormalQuantile quant;
    double width = 1.0;

    std::pair<std::pair<double, double>, std::pair<double, double> > retVal;
    retVal.first.first = H0pval/H0.size();
    retVal.first.second = quant(&retVal.first.first, &width);
    retVal.second.first = H1pval/H1.size();
    retVal.second.second = quant(&retVal.second.first, &width);

    return retVal;
  }

  template <class T>
  std::pair<std::pair<double, double>, std::pair<double, double> > p_value_CP(std::vector<T> H0, std::vector<T> H1)
  {
    for(UInt_t i = 0; i < H0.size()-2; i++)
      for(UInt_t j = 0; j < H0.size()-1-i; j++)
        if(H0[j] > H0[j+1])
          std::swap(H0[j], H0[j+1]);
    for(UInt_t i = 0; i < H1.size()-2; i++)
      for(UInt_t j = 0; j < H1.size()-1-i; j++)
        if(H1[j] > H1[j+1])
          std::swap(H1[j], H1[j+1]);

    UInt_t size = H0.size();
    T H0med = (size%2)?(H0[size/2]):((H0[size/2]+H0[size/2-1])/2);
    size = H1.size();
    T H1med = (size%2)?(H1[size/2]):((H1[size/2]+H1[size/2-1])/2);

    Int_t H0pval = 0;
    Int_t H1pval = 0;
    for(UInt_t i = 0; i < H0.size(); i++)
      if(H0[i] > H1med)
      {
        H0pval = H0.size()-(i-1);
        break;
      }
    for(UInt_t i = 0; i < H1.size(); i++)
      if(H1[i] > H0med)
      {
        H1pval = i;
        break;
      }

    separation::NormalQuantile quant;
    double width = 1.0;

    std::pair<std::pair<double, double>, std::pair<double, double> > retVal;
    retVal.first.first = TEfficiency::ClopperPearson(H0.size(), H0pval, 0.68, 1);
    retVal.first.second = TEfficiency::ClopperPearson(H0.size(), H0pval, 0.68, 0);
    retVal.second.first = TEfficiency::ClopperPearson(H1.size(), H1pval, 0.68, 1);
    retVal.second.second = TEfficiency::ClopperPearson(H1.size(), H1pval, 0.68, 0);

    return retVal;
  }

  template <class T>
  std::pair<
            std::pair<
                      std::pair<double, double>,
                      std::pair<double, double>
                      >,
            std::pair<
                      std::pair<double, double>,
                      std::pair<double, double>
                      >
            > JackKnife(std::vector<T> H0, std::vector<T> H1, std::pair<std::pair<double, double>, std::pair<double, double> > (*func)(std::vector<T>, std::vector<T>))
  {
    std::pair<
              std::pair<
                        std::pair<double, double>,
                        std::pair<double, double>
                        >,
              std::pair<
                        std::pair<double, double>,
                        std::pair<double, double>
                        >
              > retVal;

    UInt_t H0_size = H0.size()/JK_PARAM;
    UInt_t H1_size = H1.size()/JK_PARAM;

    std::random_shuffle(H0.begin(),H0.end());
    std::random_shuffle(H1.begin(),H1.end());

    std::vector<double> H0pvals, H1pvals, H0sigmas, H1sigmas;
    double H0pvalue = 0, H0sigma = 0, H1pvalue = 0, H1sigma = 0;

    for(UInt_t i = 0; i < JK_PARAM; i++)
    {
      std::vector<T> H0_small, H1_small;

      for(UInt_t j = 0; j < H0.size(); j++)
      {
        if(j < i*H0_size || j > (i+1)*H0_size-1)
        {
          H0_small.push_back(H0[j]);
        }
      }
      for(UInt_t j = 0; j < H1.size(); j++)
      {
        if(j < i*H1_size || j > (i+1)*H1_size-1)
        {
          H1_small.push_back(H1[j]);
        }
      }

      std::pair<std::pair<double, double>, std::pair<double, double> > temp = func(H0_small, H1_small);
      H0pvals.push_back(temp.first.first);
      H0pvalue += temp.first.first;
      H0sigmas.push_back(temp.first.second);
      H0sigma += temp.first.second;
      H1pvals.push_back(temp.second.first);
      H1pvalue += temp.second.first;
      H1sigmas.push_back(temp.second.second);
      H1sigma += temp.second.second;
    }

    UInt_t size = H0pvals.size();
    H0pvalue /= size;
    H1pvalue /= size;
    H0sigma /= size;
    H1sigma /= size;


    double H0pvalueErr = 0, H0sigmaErr = 0, H1pvalueErr = 0, H1sigmaErr = 0;
    for(UInt_t i = 0; i < H0pvals.size(); i++)
    {
      H0pvals[i] -= H0pvalue;
      H1pvals[i] -= H1pvalue;
      H0sigmas[i] -= H0sigma;
      H1sigmas[i] -= H1sigma;

      H0pvals[i] *= H0pvals[i];
      H1pvals[i] *= H1pvals[i];
      H0sigmas[i] *= H0sigmas[i];
      H1sigmas[i] *= H1sigmas[i];

      H0pvalueErr += H0pvals[i];
      H1pvalueErr += H1pvals[i];
      H0sigmaErr += H0sigmas[i];
      H1sigmaErr += H1sigmas[i];
    }

    H0pvalueErr /= JK_PARAM*(JK_PARAM-1);
    H1pvalueErr /= JK_PARAM*(JK_PARAM-1);
    H0sigmaErr /= JK_PARAM*(JK_PARAM-1);
    H1sigmaErr /= JK_PARAM*(JK_PARAM-1);

    H0pvalueErr = sqrt(H0pvalueErr);
    H1pvalueErr = sqrt(H1pvalueErr);
    H0sigmaErr = sqrt(H0sigmaErr);
    H1sigmaErr = sqrt(H1sigmaErr);

    double onesigma = ROOT::Math::tdistribution_quantile_c(0.16, JK_PARAM-1);


    retVal.first.first.first = H0pvalue;
    retVal.second.first.first = H1pvalue;
    retVal.first.second.first = H0sigma;
    retVal.second.second.first = H1sigma;

    retVal.first.first.second = H0pvalueErr*onesigma;
    retVal.second.first.second = H1pvalueErr*onesigma;
    retVal.first.second.second = H0sigmaErr*onesigma;
    retVal.second.second.second = H1sigmaErr*onesigma;

    return retVal;
  }
}

map<string,double> calcSeparation(TTree *tree)
{
  map<string,double> retVal;

  double q_smtoy;
  double q_gravtoy;
  tree->SetBranchAddress("q_smtoy",&q_smtoy);
  tree->SetBranchAddress("q_gravtoy",&q_gravtoy);
  std::vector<double> SMToys;
  std::vector<double> GravToys;

  for(Int_t i=0;i<tree->GetEntries();i++)
  {
    tree->GetEntry(i);
    SMToys.push_back(q_smtoy);
    GravToys.push_back(q_gravtoy);
  }

  std::pair<std::pair<double, double>, std::pair<double, double> > p_values = separation::p_value<double>(SMToys, GravToys);

  retVal.insert(std::pair<string, double>("SMProb", p_values.first.first));
  retVal.insert(std::pair<string, double>("GravProb", p_values.second.first));
  retVal.insert(std::pair<string, double>("SMSigma", p_values.first.second));
  retVal.insert(std::pair<string, double>("GravSigma", p_values.second.second));

  separation::NormalQuantile quant;
  TF1 *func = new TF1("my_quantile", &quant, 0, 1, 1);
  double width = 1.0;

  retVal.insert(std::pair<string, double>("SMProbErr", sqrt(retVal["SMProb"]*(1-retVal["SMProb"])/SMToys.size())));
  retVal.insert(std::pair<string, double>("GravProbErr", sqrt(retVal["GravProb"]*(1-retVal["GravProb"])/GravToys.size())));
  retVal.insert(std::pair<string, double>("SMSigmaErr", abs(func->Derivative(retVal["SMProb"], &width)*retVal["SMProbErr"])));
  retVal.insert(std::pair<string, double>("GravSigmaErr", abs(func->Derivative(retVal["GravProb"], &width)*retVal["GravProbErr"])));

  std::pair<
            std::pair<
                      std::pair<double, double>,
                      std::pair<double, double>
                      >,
            std::pair<
                      std::pair<double, double>,
                      std::pair<double, double>
                      >
            > p_values_jk = separation::JackKnife<double>(SMToys, GravToys, separation::p_value<double>);

  retVal.insert(std::pair<string, double>("SMProb_jk", p_values_jk.first.first.first));
  retVal.insert(std::pair<string, double>("GravProb_jk", p_values_jk.second.first.first));
  retVal.insert(std::pair<string, double>("SMSigma_jk", p_values_jk.first.second.first));
  retVal.insert(std::pair<string, double>("GravSigma_jk", p_values_jk.second.second.first));

  retVal.insert(std::pair<string, double>("SMProbErr_jk", p_values_jk.first.first.second));
  retVal.insert(std::pair<string, double>("GravProbErr_jk", p_values_jk.second.first.second));
  retVal.insert(std::pair<string, double>("SMSigmaErr_jk", p_values_jk.first.second.second));
  retVal.insert(std::pair<string, double>("GravSigmaErr_jk", p_values_jk.second.second.second));

  p_values = separation::p_value_CP<double>(SMToys, GravToys);
  retVal.insert(std::pair<string, double>("SMProb_Upper", p_values.first.first));
  retVal.insert(std::pair<string, double>("SMProb_Lower", p_values.first.second));
  retVal.insert(std::pair<string, double>("GravProb_Upper", p_values.second.first));
  retVal.insert(std::pair<string, double>("GravProb_Lower", p_values.second.second));


  return retVal;
}

int main(int argc, char* argv[]){

  string filename;
  bool isMassFac=false;
  int nCats = -1;
  if (argc!=2 && argc!=3 && argc!=4 && argc!=5) {
    cout << "usage ./bin/diyPlot <infilename> --isMassFac --nCats" << endl;
    return 1;
  }
  else {
    filename=string(argv[1]);
    for (int i=0; i<argc; i++)
    {
      if (string(argv[i])=="--isMassFac") isMassFac=true;
      if (string(argv[i])=="--nCats" && i < argc-1)
      {
        nCats = atoi(argv[i+1]);
        i++;
      }
    }
  }
  string dir="cutBased";
  if (isMassFac) dir="massFac";
  system(Form("mkdir -p plots/%s",dir.c_str()));

  TFile *inFile = TFile::Open(filename.c_str());
  TTree *tree = (TTree*)inFile->Get("limit");

  TFile *outFile = new TFile(Form("SepStats_%s.root",dir.c_str()),"RECREATE");

  TTree * oldTree;
  bool usingFitStatus = false;
  if(tree->GetBranch("fitStatus"))
  {
    usingFitStatus = true;
    oldTree = tree;
    tree = oldTree->CopyTree("fitStatus[0] == 0 && fitStatus[1] == 0 && fitStatus[2] == 0 && fitStatus[3] == 0");
  }

  //double q_smtoy;
  //double q_gravtoy;
  int nSpinCats;
  double muSM_perCTbin[100];
  double muGRAV_perCTbin[100];

  //tree->SetBranchAddress("q_smtoy",&q_smtoy);
  //tree->SetBranchAddress("q_gravtoy",&q_gravtoy);
  tree->SetBranchAddress("nSpinCats",&nSpinCats);
  tree->SetBranchAddress("muSM_perCTbin",muSM_perCTbin);
  tree->SetBranchAddress("muGRAV_perCTbin",muGRAV_perCTbin);

  // NOTE nSpinCats should not change throughout the tree
  tree->GetEntry(0);

  vector<TH1F*> muSM_percat_hists;
  vector<TH1F*> muGRAV_percat_hists;
  for (int s=0; s<nSpinCats; s++) {
    muSM_percat_hists.push_back(new TH1F(Form("muSM_spin%d",s),Form("muSM_spin%d",s),100,-10,10));
    muGRAV_percat_hists.push_back(new TH1F(Form("muGRAV_spin%d",s),Form("muGRAV_spin%d",s),100,-10,10));
  }

  int ntoys = tree->GetEntries();

  for (int entry=0; entry<tree->GetEntries(); entry++){
    tree->GetEntry(entry);
    for (int s=0; s<nSpinCats; s++){
      muSM_percat_hists[s]->Fill(muSM_perCTbin[s]);
      muGRAV_percat_hists[s]->Fill(muGRAV_perCTbin[s]);
    }
  }

  // -----------------------
  // make model independent plot
  // -----------------------
  TPaveText pt(0.1,0.9,0.45,0.99,"NDC");
  pt.SetTextAlign(12);
  pt.SetTextSize(0.04);
  pt.SetFillColor(0);
  pt.AddText("CMS Expected");
  pt.SetBorderSize(0);
  TPaveText pt2(0.55,0.9,0.9,0.99,"NDC");
  pt2.SetTextAlign(32);
  pt2.SetTextSize(0.04);
  pt2.SetFillColor(0);
  // pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.051 fb^{-1}; #sqrt{s} = 8 TeV, L = 30.0 fb^{-1}");
  pt2.AddText(" #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  pt2.SetBorderSize(0);

  // model indep
  TCanvas *canv = new TCanvas();
  TF1 *f = new TF1("f","0.",0.,1.);
  f->SetLineColor(kBlack);
  f->SetLineWidth(2);
  TGraphErrors *graphSM = new TGraphErrors();
  TGraphErrors *graphSME = new TGraphErrors();
  TGraphErrors *graphGRAV = new TGraphErrors();
  cout << "Model Indep " << endl;
  float cTvals[5]={0.1,0.2875,0.4625,0.65,0.875};
  for (int s=0; s<nSpinCats; s++){
    muSM_percat_hists[s]->Draw();
    canv->Print(Form("plots/%s/muSM_CT%d.pdf",dir.c_str(),s));
    canv->Print(Form("plots/%s/muSM_CT%d.png",dir.c_str(),s));
    muGRAV_percat_hists[s]->Draw();
    canv->Print(Form("plots/%s/muGRAV_CT%d.pdf",dir.c_str(),s));
    canv->Print(Form("plots/%s/muGRAV_CT%d.png",dir.c_str(),s));
    cout << s << " " << muSM_percat_hists[s]->GetMean() << " " << muGRAV_percat_hists[s]->GetMean() << endl;
    //graphSM->SetPoint(s,0.5/nSpinCats+(s*1./nSpinCats),muSM_percat_hists[s]->GetMean());
    //graphSM->SetPointError(s,0.,muSM_percat_hists[s]->GetMeanError());
    //graphSME->SetPoint(s,0.5/nSpinCats+(s*1./nSpinCats),muSM_percat_hists[s]->GetMean());
    //graphSME->SetPointError(s,0.,muSM_percat_hists[s]->GetRMS());
    //graphGRAV->SetPoint(s,0.5/nSpinCats+(s*1./nSpinCats),1./muGRAV_percat_hists[s]->GetMean());
    //graphGRAV->SetPointError(s,0.,muGRAV_percat_hists[s]->GetMeanError());
    graphSM->SetPoint(s,cTvals[s],muSM_percat_hists[s]->GetMean());
    graphSM->SetPointError(s,0.,muSM_percat_hists[s]->GetMeanError());
    graphSME->SetPoint(s,cTvals[s],muSM_percat_hists[s]->GetMean());
    graphSME->SetPointError(s,0.,muSM_percat_hists[s]->GetRMS());
    graphGRAV->SetPoint(s,cTvals[s],1./muGRAV_percat_hists[s]->GetMean());
    graphGRAV->SetPointError(s,0.,muGRAV_percat_hists[s]->GetMeanError());
  }
  TH2F *h = new TH2F("h","",1,0,1.,1,-0.5,3.0);
  h->SetStats(0);
  h->SetLineColor(0);
  h->Draw();
  //graphSM->SetMarkerStyle(kOpenSquare);
  graphSM->SetMarkerColor(kRed);
  graphSM->SetLineColor(kRed);
  graphSM->SetLineWidth(2);
  graphSM->SetFillColor(kRed);
  graphSM->SetFillStyle(3004);
  graphSME->SetLineColor(kRed);
  graphSME->SetLineWidth(2);
  graphSME->SetFillColor(kRed);
  graphSME->SetFillStyle(3004);
  //graphGRAV->SetMarkerStyle(kOpenCircle);
  graphGRAV->SetMarkerColor(kBlue);
  graphGRAV->SetLineColor(kBlue);
  graphGRAV->SetLineWidth(2);
  //graphSME->GetYaxis()->SetRangeUser(-0.5,3.0);
  //graphSME->GetXaxis()->SetRangeUser(0.,1.);
  graphSME->GetXaxis()->SetTitle("|cos(#theta*_{CS}|");
  graphSME->GetYaxis()->SetTitle("#sigma_{X} / #sigma_{SM}");
  TLegend *leg3 = new TLegend(0.11,0.75,0.4,0.89);
  leg3->SetFillColor(0);
  leg3->SetLineColor(0);
  leg3->AddEntry(graphSM,"X#rightarrow#gamma#gamma 0^{+}","LPF");
  leg3->AddEntry(graphGRAV,"X#rightarrow#gamma#gamma 2^{+}_{m}","LEP");
  graphSME->Draw("E3same");
  graphSM->Draw("LEPsame");
  graphGRAV->Draw("LEPsame");
  f->Draw("same");
  leg3->Draw("same");
  pt.Draw("same");
  pt2.Draw("same");
  canv->Print(Form("plots/%s/modIndep.pdf",dir.c_str()));
  canv->Print(Form("plots/%s/modIndep.png",dir.c_str()));
  canv->Write();

  // -----------------------
  // make plots of simultaneous mu
  // -----------------------
  TH1F *muSMSMhist = new TH1F("muSMSMH","muSMSMH",120,-2,4);
  TH1F *muSMGRAVhist = new TH1F("muSMGRAVH","muSMGRAVH",120,-2,4);
  TH1F *muGRAVSMhist = new TH1F("muGRAVSMH","muGRAVSMH",120,-2,4);
  TH1F *muGRAVGRAVhist = new TH1F("muGRAVGRAVH","muGRAVGRAVH",120,-2,4);

  tree->Draw("muSMSM>>muSMSMH","muSMSM>=-10. && muSMSM<=10.","goff");
  tree->Draw("muSMGRAV>>muSMGRAVH","muSMGRAV>=-10. && muSMGRAV<=10.","goff");
  tree->Draw("muGRAVSM>>muGRAVSMH","muGRAVSM>=-10. && muGRAVSM<=10.","goff");
  tree->Draw("muGRAVGRAV>>muGRAVGRAVH","muGRAVGRAV>=-10. && muGRAVGRAV<=10.","goff");

  gStyle->SetOptStat(0);

  muSMSMhist->SetLineColor(kBlue);
  muSMSMhist->SetFillColor(kBlue-9);
  muSMSMhist->SetTitle("SM PDF");
  muSMGRAVhist->SetLineColor(kRed);
  muSMGRAVhist->SetLineWidth(2);
  muSMGRAVhist->SetFillColor(kRed);
  muSMGRAVhist->SetFillStyle(3495);
  muSMGRAVhist->SetTitle("SM PDF");

  muGRAVSMhist->SetLineColor(kBlue);
  muGRAVSMhist->SetFillColor(kBlue-9);
  muGRAVSMhist->SetTitle("GRAV PDF");
  muGRAVGRAVhist->SetLineColor(kRed);
  muGRAVGRAVhist->SetLineWidth(2);
  muGRAVGRAVhist->SetFillColor(kRed);
  muGRAVGRAVhist->SetFillStyle(3495);
  muGRAVGRAVhist->SetTitle("GRAV PDF");

  TLegend *leg = new TLegend(0.2,0.6,0.4,0.89);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(muSMSMhist,"SM toy","F");
  leg->AddEntry(muSMGRAVhist,"GRAV toy","F");

  muSMSMhist->Draw("HIST");
  muSMGRAVhist->Draw("HISTFSAME");
  leg->Draw("SAME");
  canv->Print(Form("plots/%s/mu_SM.pdf",dir.c_str()));
  canv->Print(Form("plots/%s/mu_SM.png",dir.c_str()));
  muGRAVSMhist->Draw("HISTF");
  muGRAVGRAVhist->Draw("HISTSAME");
  leg->Draw("SAME");
  canv->Print(Form("plots/%s/mu_GRAV.pdf",dir.c_str()));
  canv->Print(Form("plots/%s/mu_GRAV.png",dir.c_str()));

  // -----------------------
  // make likelihood plot
  // -----------------------
  double lmin=-15.;
  double lmax=15.;
  TH1F *testStatSMH = new TH1F("testStatSMH","testStatSMH",600,lmin,lmax);
  TH1F *testStatGRAVH = new TH1F("testStatGRAVH","testStatGRAVH",600,lmin,lmax);
  //TH1F *testStatOBSH = new TH1F("testStatOBSH","testStatOBSH",600,lmin,lmax);

  tree->Draw("q_smtoy>>testStatSMH","q_smtoy>=-900. && q_smtoy<=900.","goff");
  tree->Draw("q_gravtoy>>testStatGRAVH","q_gravtoy>=-900. && q_gravtoy<=900.","goff");

  separation::NormalQuantile quant;
  TF1 *func = new TF1("quant", &quant, 0, 1, 1);
  double width = 1.0;

  // calc prob from hist
  double SMprobHist = testStatSMH->Integral(testStatSMH->FindBin(testStatGRAVH->GetMean()),testStatSMH->GetNbinsX())/testStatSMH->Integral();
  double SMprobHistErr = sqrt(SMprobHist*(1-SMprobHist)/testStatSMH->Integral());
  double GRAVprobHist = testStatGRAVH->Integral(0,testStatGRAVH->FindBin(testStatSMH->GetMean()))/testStatGRAVH->Integral();
  double GRAVprobHistErr = sqrt(GRAVprobHist*(1-GRAVprobHist)/testStatGRAVH->Integral());

  // now rebin histograms so they look better on a plot
  testStatSMH->Rebin(10);
  testStatGRAVH->Rebin(10);

  // fit TH1Fs with a gaussian
  TF1 *gausSM = new TF1("gausSM","gaus",lmin,lmax);
  TF1 *gausGRAV = new TF1("gausGRAV","gaus",lmin,lmax);

  testStatSMH->Fit(gausSM,"QN");
  testStatGRAVH->Fit(gausGRAV,"QN");

  // calc prob from fit
  double SMmean = gausSM->GetParameter(1);
  double GRAVmean = gausGRAV->GetParameter(1);

  double SMprob = gausSM->Integral(GRAVmean,15.)/gausSM->Integral(-15.,15.);
  double SMprobErr = sqrt(SMprob*(1-SMprob)/testStatSMH->Integral());
  double GRAVprob = gausGRAV->Integral(-15.,SMmean)/gausGRAV->Integral(-15.,15.);
  double GRAVprobErr = sqrt(GRAVprob*(1-GRAVprob)/testStatGRAVH->Integral());

  double SMsigma = ROOT::Math::normal_quantile_c(SMprob,width);
  double SMsigmaErr = abs(func->Derivative(SMprob, &width)*SMprobErr);
  double GRAVsigma = ROOT::Math::normal_quantile_c(GRAVprob,width);
  double GRAVsigmaErr = abs(func->Derivative(GRAVprob, &width)*GRAVprobErr);

  double SMsigmaHist = ROOT::Math::normal_quantile_c(SMprobHist,width);
  double SMsigmaHistErr = abs(func->Derivative(SMprobHist, &width)*SMprobHistErr);
  double GRAVsigmaHist = ROOT::Math::normal_quantile_c(GRAVprobHist,width);
  double GRAVsigmaHistErr = abs(func->Derivative(GRAVprobHist, &width)*GRAVprobHistErr);

  cout << "SM: " << SMmean << " GRAV: " << GRAVmean << endl;
  cout << "UNBINNED:" << endl;
  cout << "Prob( q > median(2) | 0 ) = "<< SMprob << " = " << SMsigma << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprob << " = " << GRAVsigma << " sigma " << endl;
  //cout << "BINNED:  -- ErfcInverse" << endl;
  //cout << "Prob( q > median(2) | 0 ) = "<< SMprobHist << " = " << SMsigmaHist << " sigma " << endl;
  //cout << "Prob( q < median(0) | 2 ) = "<< GRAVprobHist << " = " << GRAVsigmaHist << " sigma " << endl;
  cout << "BINNED:  -- normal_quantile_c" << endl;
  cout << "Prob( q > median(2) | 0 ) = "<< SMprobHist << " = " << SMsigmaHist << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprobHist << " = " << GRAVsigmaHist << " sigma " << endl;
  cout << "BINNED:  -- normal_quantile_c" << endl;
  cout << "Prob( q > median(2) | 0 ) = "<< SMprobHist << " = " << ROOT::Math::normal_quantile_c(SMprobHist,1.0) << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprobHist << " = " << ROOT::Math::normal_quantile_c(GRAVprobHist,1.0) << " sigma " << endl;
  cout << "BINNED:  -- RooStats::PValueToSignificance" << endl;
  cout << "Prob( q > median(2) | 0 ) = "<< SMprobHist << " = " << RooStats::PValueToSignificance(SMprobHist) << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprobHist << " = " << RooStats::PValueToSignificance(GRAVprobHist) << " sigma " << endl;

  map<string,double> temp = calcSeparation(tree);

  cout << "UNBINNED:  -- normal_quantile_c" << endl;
  cout << "Prob( q > median(2) | 0 ) = " << temp["SMProb"] << "+-" << temp["SMProbErr"] << " = " << temp["SMSigma"] << "+-" << temp["SMSigmaErr"] << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = " << temp["GravProb"] << "+-" << temp["GravProbErr"] << " = " << temp["GravSigma"] << "+-" << temp["GravSigmaErr"] << " sigma " << endl;

  cout << "UNBINNED_JK:  -- normal_quantile_c" << endl;
  cout << "Prob( q > median(2) | 0 ) = " << temp["SMProb_jk"] << "+-" << temp["SMProbErr_jk"] << " = " << temp["SMSigma_jk"] << "+-" << temp["SMSigmaErr_jk"] << " sigma " << endl;
  cout << "Prob( q < median(0) | 2 ) = " << temp["GravProb_jk"] << "+-" << temp["GravProbErr_jk"] << " = " << temp["GravSigma_jk"] << "+-" << temp["GravSigmaErr_jk"] << " sigma " << endl;

  if(nCats >= 0)
  {
    TFile *haddable = new TFile(Form("%dCats_separation.root",nCats),"RECREATE");
    TTree *tree = new TTree("Separation","Separation");



    tree->Branch("nCats", &nCats);


    tree->Branch("FittedSMProb", &SMprob);
    tree->Branch("FittedGravProb", &GRAVprob);
    tree->Branch("FittedSMProbErr", &SMprobErr);
    tree->Branch("FittedGravProbErr", &GRAVprobErr);

    tree->Branch("UnbinnedSMProb", &temp["SMProb"]);
    tree->Branch("UnbinnedGravProb", &temp["GravProb"]);
    tree->Branch("UnbinnedSMProbCPUpper", &temp["SMProb_Upper"]);
    tree->Branch("UnbinnedSMProbCPLower", &temp["SMProb_Lower"]);
    tree->Branch("UnbinnedGravProbCPUpper", &temp["GravProb_Upper"]);
    tree->Branch("UnbinnedGravProbCPLower", &temp["GravProb_Lower"]);
    tree->Branch("UnbinnedSMProbErr", &temp["SMProbErr"]);
    tree->Branch("UnbinnedGravProbErr", &temp["GravProbErr"]);

    tree->Branch("BinnedSMProb", &SMprobHist);
    tree->Branch("BinnedGravProb", &GRAVprobHist);
    tree->Branch("BinnedSMProbErr", &SMprobHistErr);
    tree->Branch("BinnedGravProbErr", &GRAVprobHistErr);


    tree->Branch("FittedSMSigma", &SMsigma);
    tree->Branch("FittedGravSigma", &GRAVsigma);
    tree->Branch("FittedSMSigmaErr", &SMsigmaErr);
    tree->Branch("FittedGravSigmaErr", &GRAVsigmaErr);

    tree->Branch("UnbinnedSMSigma", &temp["SMSigma"]);
    tree->Branch("UnbinnedGravSigma", &temp["GravSigma"]);
    tree->Branch("UnbinnedSMSigmaErr", &temp["SMSigmaErr"]);
    tree->Branch("UnbinnedGravSigmaErr", &temp["GravSigmaErr"]);

    tree->Branch("BinnedSMSigma", &SMsigmaHist);
    tree->Branch("BinnedGravSigma", &GRAVsigmaHist);
    tree->Branch("BinnedSMSigmaErr", &SMsigmaHistErr);
    tree->Branch("BinnedGravSigmaErr", &GRAVsigmaHistErr);



    tree->Fill();
    //tree->Print();
    haddable->Write();

    outFile->cd();
    haddable->Close();
  }

  // set style
  testStatSMH->SetLineColor(kMagenta-3);
  testStatGRAVH->SetLineColor(kBlue+1);
  testStatSMH->SetTitle("");
  testStatGRAVH->SetTitle("");
  testStatGRAVH->GetXaxis()->SetTitle("S = -2 #times ln(L_{SM}/L_{GRAV})");
  testStatSMH->GetXaxis()->SetTitle("S = -2 #times ln(L_{SM}/L_{GRAV})");
  testStatGRAVH->GetYaxis()->SetTitle("Number of toys");
  testStatSMH->GetYaxis()->SetTitle("Number of toys");
  testStatSMH->GetYaxis()->SetTitleOffset(1.2);
  testStatGRAVH->GetYaxis()->SetTitleOffset(1.2);
  testStatSMH->GetXaxis()->SetTitleOffset(1.1);
  testStatGRAVH->GetXaxis()->SetTitleOffset(1.1);
  testStatGRAVH->Draw("HIST");
  testStatSMH->Draw("HISTSAME");
  gausSM->SetLineWidth(3);
  gausSM->SetLineColor(kMagenta-3);
  gausGRAV->SetLineWidth(3);
  gausGRAV->SetLineColor(kBlue+1);
  gausSM->Draw("CSAME");
  gausGRAV->Draw("CSAME");
  // now clone to draw fill areas
  TF1 *gausSMclone = (TF1*)gausSM->Clone();
  TF1 *gausGRAVclone = (TF1*)gausGRAV->Clone();
  gausSMclone->SetFillColor(kMagenta-7);
  gausSMclone->SetFillStyle(3001);
  gausGRAVclone->SetFillColor(kBlue-7);
  gausGRAVclone->SetFillStyle(3004);
  gausSMclone->SetRange(GRAVmean,lmax);
  gausGRAVclone->SetRange(lmin,SMmean);
  gausSMclone->Draw("FCSAME");
  gausGRAVclone->Draw("FCSAME");

  pt.Draw();
  pt2.Draw();
  TPaveText pt3(0.7,0.85,0.89,0.89,"NDC");
  pt3.AddText(Form("N = %d toys",ntoys));
  pt3.SetBorderSize(0);
  pt3.SetFillColor(0);
  pt3.Draw();
  TPaveText pt4(0.6,0.7,0.89,0.85,"NDC");
  //pt4.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",RooStats::PValueToSignificance(GRAVprobHist)));
  //pt4.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",RooStats::PValueToSignificance(SMprobHist)));
  pt4.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",ROOT::Math::normal_quantile_c(GRAVprob,1.0)));
  pt4.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",ROOT::Math::normal_quantile_c(SMprob,1.0)));
  pt4.SetBorderSize(0);
  pt4.SetFillColor(0);
  pt4.Draw();

  TLegend *leg2 = new TLegend(0.11,0.75,0.4,0.89);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(gausSMclone,"X#rightarrow#gamma#gamma 0^{+}","f");
  leg2->AddEntry(gausGRAVclone,"X#rightarrow#gamma#gamma 2^{+}_{m}","f");
  leg2->Draw();
  canv->Print(Form("plots/%s/fitted_sep.pdf",dir.c_str()));
  canv->Print(Form("plots/%s/fitted_sep.png",dir.c_str()));

  testStatSMH->SetFillColor(kMagenta-7);
  testStatSMH->SetFillStyle(3001);
  testStatGRAVH->SetFillColor(kBlue-7);
  testStatGRAVH->SetFillStyle(3004);
  testStatSMH->Draw("HISTF");
  testStatGRAVH->Draw("HISTFSAME");
  leg2->Draw();
  pt.Draw();
  pt2.Draw();
  pt3.Draw();
  TPaveText pt5(0.6,0.7,0.89,0.85,"NDC");
  //pt5.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",RooStats::PValueToSignificance(GRAVprobHist)));
  //pt5.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",RooStats::PValueToSignificance(SMprobHist)));
  pt5.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",ROOT::Math::normal_quantile_c(GRAVprobHist,1.0)));
  pt5.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",ROOT::Math::normal_quantile_c(SMprobHist,1.0)));
  pt5.SetBorderSize(0);
  pt5.SetFillColor(0);
  pt5.Draw();
  canv->Print(Form("plots/%s/binned_sep.pdf",dir.c_str()));
  canv->Print(Form("plots/%s/binned_sep.png",dir.c_str()));

  outFile->cd();
  muSMSMhist->Write();
  muSMGRAVhist->Write();
  muGRAVSMhist->Write();
  muGRAVGRAVhist->Write();
  testStatSMH->Write();
  testStatGRAVH->Write();
  //testStatOBSH->Write();

  outFile->Close();
  inFile->Close();

  return 0;
}

