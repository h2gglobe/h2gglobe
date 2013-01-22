#include <iostream>
#include <string>
#include <vector>
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
#include "TStyle.h"
#include "TMath.h"
#include "Math/DistFunc.h"

using namespace std;

int main(int argc, char* argv[]){
 
  string filename;
  bool isMassFac=false;
  if (argc!=2 && argc!=3) {
    cout << "usage ./bin/diyPlot <infilename> --isMassFac" << endl;
    return 1;
  }
  else {
    filename=string(argv[1]);
    for (int i=0; i<argc; i++) if (string(argv[i])=="--isMassFac") isMassFac=true;
  }
  string dir="cutBased";
  if (isMassFac) dir="massFac";
  system(Form("mkdir -p plots/%s",dir.c_str()));

  TFile *inFile = TFile::Open(filename.c_str());
  TTree *tree = (TTree*)inFile->Get("limit");

  TFile *outFile = new TFile(Form("SepStats_%s.root",dir.c_str()),"RECREATE");

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
  for (int s=0; s<nSpinCats; s++){
    muSM_percat_hists[s]->Draw();
    canv->Print(Form("plots/%s/muSM_CT%d.pdf",dir.c_str(),s));
    canv->Print(Form("plots/%s/muSM_CT%d.png",dir.c_str(),s));
    muGRAV_percat_hists[s]->Draw();
    canv->Print(Form("plots/%s/muGRAV_CT%d.pdf",dir.c_str(),s));
    canv->Print(Form("plots/%s/muGRAV_CT%d.png",dir.c_str(),s));
    cout << s << " " << muSM_percat_hists[s]->GetMean() << " " << muGRAV_percat_hists[s]->GetMean() << endl;
    graphSM->SetPoint(s,0.1+(s*0.2),muSM_percat_hists[s]->GetMean());
    graphSM->SetPointError(s,0.,muSM_percat_hists[s]->GetMeanError());
    graphSME->SetPoint(s,0.1+(s*0.2),muSM_percat_hists[s]->GetMean());
    graphSME->SetPointError(s,0.,muSM_percat_hists[s]->GetRMS());
    graphGRAV->SetPoint(s,0.1+(s*0.2),1./muGRAV_percat_hists[s]->GetMean());
    graphGRAV->SetPointError(s,0.,muGRAV_percat_hists[s]->GetMeanError());
  }
  TH2F *h = new TH2F("h","",1,0,1.,1,-0.5,3.0);
  h->SetStats(0);
  h->SetLineColor(0);
  h->Draw();
  graphSM->SetMarkerStyle(kOpenSquare);
  graphSM->SetMarkerColor(kRed);
  graphSM->SetLineColor(kRed);
  graphSM->SetLineWidth(2);
  graphSM->SetFillColor(kRed);
  graphSM->SetFillStyle(3004);
  graphSME->SetLineColor(kRed);
  graphSME->SetLineWidth(2);
  graphSME->SetFillColor(kRed);
  graphSME->SetFillStyle(3004);
  graphGRAV->SetMarkerStyle(kOpenCircle);
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
  graphSM->Draw("LPsame");
  graphGRAV->Draw("LPsame");
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
 
  // calc prob from hist
  double SMprobHist = testStatSMH->Integral(testStatSMH->FindBin(testStatGRAVH->GetMean()),testStatSMH->GetNbinsX())/testStatSMH->Integral();
  double GRAVprobHist = testStatGRAVH->Integral(0,testStatGRAVH->FindBin(testStatSMH->GetMean()))/testStatGRAVH->Integral();
  double SMsigmaHist = TMath::ErfcInverse(SMprobHist)*TMath::Sqrt(2.0);
  double GRAVsigmaHist = TMath::ErfcInverse(GRAVprobHist)*TMath::Sqrt(2.0);

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
  double GRAVprob = gausGRAV->Integral(-15.,SMmean)/gausGRAV->Integral(-15.,15.);

  //double SMsigma = ROOT::Math::normal_quantile_c(SMprob,1.0);
  //double GRAVsigma = ROOT::Math::normal_quantile_c(GRAVprob,1.0);
  double SMsigma = TMath::ErfcInverse(SMprob)*TMath::Sqrt(2.0);
  double GRAVsigma = TMath::ErfcInverse(GRAVprob)*TMath::Sqrt(2.0);

  cout << "SM: " << SMmean << " GRAV: " << GRAVmean << endl;
  cout << "UNBINNED:" << endl; 
  cout << "Prob( q > median(2) | 0 ) = "<< SMprob << " = " << SMsigma << " sigma " << endl; 
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprob << " = " << GRAVsigma << " sigma " << endl; 
  cout << "BINNED:" << endl; 
  cout << "Prob( q > median(2) | 0 ) = "<< SMprobHist << " = " << SMsigmaHist << " sigma " << endl; 
  cout << "Prob( q < median(0) | 2 ) = "<< GRAVprobHist << " = " << GRAVsigmaHist << " sigma " << endl; 

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
  pt4.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",GRAVsigma));
  pt4.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",SMsigma));
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
  pt5.AddText(Form("p (q<med(2^{+}_{m}) | 0^{+}) = %4.2f#sigma",GRAVsigmaHist));
  pt5.AddText(Form("p (q>med(0^{+}) | 2^{+}_{m}) = %4.2f#sigma",SMsigmaHist));
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

