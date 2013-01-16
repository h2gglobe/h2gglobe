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
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Math/DistFunc.h"

using namespace std;

void extractSignif(TH1F* hSM, vector<float> &v_SM, TH1F* hPS, vector<float> &v_PS, TH1F *hObs, vector<float> &v_Obs){
  
  hSM->SetName("hSM;S = -2 #times ln(L_{1}/L_{2});Number of Toys");
  hPS->SetName("hPS;S = -2 #times ln(L_{1}/L_{2});Number of Toys");
  hObs->SetName("hObserved");
  hSM->SetTitle("");
  hPS->SetTitle("");
  hObs->SetTitle("");

  sort(v_SM.begin(),v_SM.end());//sort in ascending order
  sort(v_PS.begin(),v_PS.end()); 
  sort(v_Obs.begin(),v_Obs.end());
  int ntoysSM= hSM->GetEntries();
  int ntoysPS= hPS->GetEntries();

  //we assume that SM is on the right and PS on the left of zero
  if(v_PS.at(0)>v_SM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model shouldstay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return;
  }
  float medianSM=v_SM.at(int(ntoysSM/2));
  float medianPS=v_PS.at(int(ntoysPS/2));
  cout<<"Toys generated "<<ntoysSM<<"\t"<<ntoysPS<<endl;
  cout<<"Mean of SM/PS hypothesis: "<<hSM->GetMean()<<"   /   "<<hPS->GetMean()<<endl;
  cout<<"RMS  of SM/PS hypothesis: "<<hSM->GetRMS()<<"   /   "<<hPS->GetRMS()<<endl;
  cout<<"Median of SM/PS hypothesis: "<<medianSM<<"   /   "<<medianPS<<endl;

  float coverage=0.0;
  float diff=10.0;

  float integralSM=hSM->Integral();
  float integralPS=hPS->Integral();
  cout << "Integral SM/PS: " << integralSM << " / " << integralPS << endl;

  float tailSM=hSM->Integral(1,hSM->FindBin(medianPS))/integralSM;
  float tailPS=hPS->Integral(hPS->FindBin(medianSM),hPS->GetNbinsX())/integralPS;

  cout << "Tails SM/PS: " << tailSM << " / " << tailPS << endl;
  cout<<"Prob( q < median(P) | S ) = "<<tailSM<<"  ("<<ROOT::Math::normal_quantile_c(tailSM,1.0) <<" sigma)"<<endl;
  cout<<"Prob( q > median(S) | P ) = "<<tailPS<<"  ("<<ROOT::Math::normal_quantile_c(tailPS,1.0) <<" sigma)"<<endl;

  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracPS=hPS->Integral(i,hPS->GetNbinsX()) / integralPS;
    if(fabs(fracSM-fracPS)<diff){
      diff=fabs(fracSM-fracPS);
      coverage=(fracSM+fracPS)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;

  //Fancy plot
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","c1",800,800);
  c1->cd();
  hSM->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hSM->SetYTitle("Generated experiments");
  hSM->GetYaxis()->SetTitleOffset(1.15);
  hPS->SetXTitle("S = -2 #times ln(L_{1}/L_{2})");
  hPS->SetYTitle("Generated experiments");
  hSM->SetLineColor(kMagenta-3);
  hSM->SetFillColor(kMagenta-3);
  hSM->SetLineWidth(2);
  hSM->SetFillStyle(3605);
  hPS->SetLineColor(kBlue+1);
  hPS->SetFillColor(kBlue+1);
  hPS->SetLineWidth(2);
  hPS->SetFillStyle(3695);
  hObs->SetLineColor(kGreen+3);
  hObs->SetLineWidth(2);
  hSM->Draw();
  hPS->Draw("sames");
  //  hObs->Draw("sames");

  TLegend *leg = new TLegend(0.15,0.65,0.4,0.89);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSM,"  SM, Spin-0","f");
  leg->AddEntry(hPS,"  Graviton, Spin-2","f");
  leg->Draw();


  TPaveText pt(0.16,0.95,0.45,0.99,"NDC");
  pt.SetFillColor(0);
  pt.AddText("CMS Expected");
  pt.SetBorderSize(0);
  TPaveText pt2(0.55,0.95,0.99,0.99,"NDC");
  pt2.SetFillColor(0);
  // pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.051 fb^{-1}; #sqrt{s} = 8 TeV, L = 30.0 fb^{-1}");
  pt2.AddText(" #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  pt2.SetBorderSize(0);
  pt.Draw();
  pt2.Draw();
  c1->SaveAs("sigsep_combine.pdf");
}

int main(int argc, char* argv[]){
 
  string filename;
  if (argc!=2) {
    cout << "usage ./bin/diyPlot <infilename>" << endl;
    return 1;
  }
  else {
    filename=string(argv[1]);
  }
  TFile *inFile = TFile::Open(filename.c_str());
  TTree *tree = (TTree*)inFile->Get("limit");

  TFile *outFile = new TFile("SepStats.root","RECREATE");

  double q_smtoy;
  double q_gravtoy;

  tree->SetBranchAddress("q_smtoy",&q_smtoy);
  tree->SetBranchAddress("q_gravtoy",&q_gravtoy);

  TH1F *muSMSMhist = new TH1F("muSMSMH","muSMSMH",100,-10,10);
  TH1F *muSMGRAVhist = new TH1F("muSMGRAVH","muSMGRAVH",100,-10,10);
  TH1F *muGRAVSMhist = new TH1F("muGRAVSMH","muGRAVSMH",100,-10,10);
  TH1F *muGRAVGRAVhist = new TH1F("muGRAVGRAVH","muGRAVGRAVH",100,-10,10);

  TH1F *testStatSMH = new TH1F("testStatSMH","testStatSMH",50,-10,10);
  TH1F *testStatGRAVH = new TH1F("testStatGRAVH","testStatGRAVH",50,-10,10);
  TH1F *testStatOBSH = new TH1F("testStatOBSH","testStatOBSH",50,-10,10);
  
  vector<float> v_SM;
  vector<float> v_GRAV;
  vector<float> v_OBS;

  int ntoys = tree->GetEntries();

  tree->Draw("muSMSM>>muSMSMH","muSMSM>=-10. && muSMSM<=10.","goff");
  tree->Draw("muSMGRAV>>muSMGRAVH","muSMGRAV>=-10. && muSMGRAV<=10.","goff");
  tree->Draw("muGRAVSM>>muGRAVSMH","muGRAVSM>=-10. && muGRAVSM<=10.","goff");
  tree->Draw("muGRAVGRAV>>muGRAVGRAVH","muGRAVGRAV>=-10. && muGRAVGRAV<=10.","goff");

  tree->Draw("q_smtoy>>testStatSMH","q_smtoy>=-1000. && q_smtoy<=1000.","goff");
  tree->Draw("q_gravtoy>>testStatGRAVH","q_gravtoy>=-1000. && q_gravtoy<=1000.","goff");

  gStyle->SetOptStat(0);

  muSMSMhist->SetLineColor(kBlue);
  muSMSMhist->SetFillColor(kBlue-9);
  muSMSMhist->SetTitle("SM PDF");
  muSMSMhist->GetXaxis()->SetRangeUser(-2.,4.);
  muSMGRAVhist->SetLineColor(kRed);
  muSMGRAVhist->SetLineWidth(2);
  muSMGRAVhist->SetFillColor(kRed);
  muSMGRAVhist->SetFillStyle(3495);
  muSMGRAVhist->SetTitle("SM PDF");
  muSMGRAVhist->GetXaxis()->SetRangeUser(-2.,4.);

  muGRAVSMhist->SetLineColor(kBlue);
  muGRAVSMhist->SetFillColor(kBlue-9);
  muGRAVSMhist->SetTitle("GRAV PDF");
  muGRAVSMhist->GetXaxis()->SetRangeUser(-2.,4.);
  muGRAVGRAVhist->SetLineColor(kRed);
  muGRAVGRAVhist->SetLineWidth(2);
  muGRAVGRAVhist->SetFillColor(kRed);
  muGRAVGRAVhist->SetFillStyle(3495);
  muGRAVGRAVhist->SetTitle("GRAV PDF");
  muGRAVGRAVhist->GetXaxis()->SetRangeUser(-2.,4.);

  TLegend *leg = new TLegend(0.2,0.6,0.4,0.89);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(muSMSMhist,"SM toy","F");
  leg->AddEntry(muSMGRAVhist,"GRAV toy","F");

  TCanvas *canv = new TCanvas();
  muSMSMhist->Draw("HIST");
  muSMGRAVhist->Draw("HISTFSAME");
  leg->Draw("SAME");
  canv->Print("mu_SM.pdf");
  muGRAVSMhist->Draw("HISTF");
  muGRAVGRAVhist->Draw("HISTSAME");
  leg->Draw("SAME");
  canv->Print("mu_GRAV.pdf");

  for (int entry=0; entry<tree->GetEntries(); entry++){
    tree->GetEntry(entry);
    v_SM.push_back(float(q_smtoy));
    v_GRAV.push_back(float(q_gravtoy));
  }

  TCanvas *canv2 = new TCanvas("c2","c2",600,600);
  
  TF1 *gausSM = new TF1("gausSM","gaus",-10,10);
  TF1 *gausGRAV = new TF1("gausGRAV","gaus",-10,10);
  testStatSMH->Fit(gausSM,"N");
  testStatGRAVH->Fit(gausGRAV,"N");

  testStatSMH->SetLineColor(kMagenta-3);
  testStatGRAVH->SetLineColor(kBlue+1);
  //testStatSMH->SetFillColor(kMagenta-3);
  //testStatGRAVH->SetFillColor(kBlue+1);
  //testStatSMH->SetFillStyle(3605);
  //testStatGRAVH->SetFillStyle(3695);
  testStatGRAVH->SetTitle("");
  testStatSMH->SetTitle("");
  testStatGRAVH->GetXaxis()->SetTitle("S = -2 #times ln(L_{SM}/L_{GRAV})");
  testStatSMH->GetXaxis()->SetTitle("S = -2 #times ln(L_{SM}/L_{GRAV})");
  testStatGRAVH->GetYaxis()->SetTitle("Number of toys");
  testStatSMH->GetYaxis()->SetTitle("Number of toys");
  testStatSMH->GetYaxis()->SetTitleOffset(1.2);
  testStatGRAVH->GetYaxis()->SetTitleOffset(1.2);
  testStatSMH->GetXaxis()->SetTitleOffset(1.1);
  testStatGRAVH->GetXaxis()->SetTitleOffset(1.1);
  testStatGRAVH->Draw("HISTF");
  testStatSMH->Draw("HISTFSAME");
  gausSM->SetLineWidth(3);
  gausSM->SetLineColor(kMagenta-3);
  gausGRAV->SetLineWidth(3);
  gausGRAV->SetLineColor(kBlue+1);
  //gausSM->Draw("SAME");
  //gausGRAV->Draw("SAME");

  TH1F *testStatSMHFit   = new TH1F("testStatSMHFit","testStatSMHFit",1000,-10.,10.); 
  TH1F *testStatGRAVHFit   = new TH1F("testStatGRAVHFit","testStatGRAVHFit",1000,-10.,10.); 
  
  for (int b=1; b<=testStatSMHFit->GetNbinsX(); b++){
    
    testStatSMHFit->SetBinContent(b,gausSM->Eval(testStatSMHFit->GetBinCenter(b)));
    testStatGRAVHFit->SetBinContent(b,gausGRAV->Eval(testStatGRAVHFit->GetBinCenter(b)));
  }
 
  testStatSMHFit->SetLineColor(kMagenta-3);
  testStatGRAVHFit->SetLineColor(kBlue+1);
  testStatSMHFit->SetLineWidth(3);
  testStatGRAVHFit->SetLineWidth(3);
  testStatSMHFit->Draw("HISTFSAME");
  testStatGRAVHFit->Draw("HISTFSAME");

  double SMmean = gausSM->GetParameter(1);
  double GRAVmean = gausGRAV->GetParameter(1);

  double SMprob = gausSM->Integral(GRAVmean,10.)/gausSM->Integral(-10.,10.);
  double GRAVprob = gausGRAV->Integral(-10.,SMmean)/gausGRAV->Integral(-10.,10.);

  double SMsigma = ROOT::Math::normal_quantile_c(SMprob,1.0);
  double GRAVsigma = ROOT::Math::normal_quantile_c(GRAVprob,1.0);

  TH1F *testStatSMHFitClone = (TH1F*)testStatSMHFit->Clone(Form("%s_clone",testStatSMHFit->GetName()));
  TH1F *testStatGRAVHFitClone = (TH1F*)testStatGRAVHFit->Clone(Form("%s_clone",testStatGRAVHFit->GetName()));

  cout << "SM: " << SMmean << " GRAV: " << GRAVmean << endl;
  cout << "SM: " << testStatSMHFitClone->FindBin(SMmean) << " GRAV: " << testStatGRAVHFitClone->FindBin(GRAVmean) << endl;
  // get frac of GRAV < SM median
  for (int b=testStatSMHFitClone->FindBin(SMmean); b<=testStatSMHFitClone->GetNbinsX(); b++) testStatGRAVHFitClone->SetBinContent(b,0.);
  // get frac of SM > GRAV median
  for (int b=1; b<=testStatGRAVHFitClone->FindBin(GRAVmean); b++) testStatSMHFitClone->SetBinContent(b,0.);
  testStatSMHFitClone->SetFillColor(kMagenta-7);
  testStatSMHFitClone->SetFillStyle(4050);
  testStatGRAVHFitClone->SetFillColor(kBlue-7);
  testStatGRAVHFitClone->SetFillStyle(4050);
  testStatSMHFitClone->Draw("HISTFSAME");
  testStatGRAVHFitClone->Draw("HISTFSAME");

  
  TPaveText pt(0.16,0.95,0.45,0.99,"NDC");
  pt.SetFillColor(0);
  pt.AddText("CMS Expected");
  pt.SetBorderSize(0);
  TPaveText pt2(0.55,0.95,0.99,0.99,"NDC");
  pt2.SetFillColor(0);
  // pt2.AddText(" #sqrt{s} = 7 TeV, L = 5.051 fb^{-1}; #sqrt{s} = 8 TeV, L = 30.0 fb^{-1}");
  pt2.AddText(" #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  pt2.SetBorderSize(0);
  pt.Draw();
  pt2.Draw();
  TPaveText pt3(0.7,0.85,0.89,0.89,"NDC");
  pt3.AddText(Form("N = %d toys",ntoys));
  pt3.SetBorderSize(0);
  pt3.SetFillColor(0);
  pt3.Draw();
  TPaveText pt4(0.65,0.75,0.89,0.85,"NDC");
  pt4.AddText(Form("Separation ~%1.1f#sigma",GRAVsigma));
  pt4.SetBorderSize(0);
  pt4.SetFillColor(0);
  pt4.Draw();

  TLegend *leg2 = new TLegend(0.15,0.65,0.4,0.89);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(testStatSMHFitClone,"  SM, Spin-0","f");
  leg2->AddEntry(testStatGRAVHFitClone,"  Graviton, Spin-2","f");
  leg2->Draw();
  canv2->Print("fitted_sep.pdf");
  
  cout << "UNBINNED:" << endl; 
  cout << "Prob( q < median(P) | S ) = "<< SMprob << " = " << SMsigma << " sigma " << endl; 
  cout << "Prob( q > median(S) | P ) = "<< GRAVprob << " = " << GRAVsigma << " sigma " << endl; 

  extractSignif(testStatSMH,v_SM,testStatGRAVH,v_GRAV,testStatOBSH,v_OBS);
  
  outFile->cd();
  muSMSMhist->Write();
  muSMGRAVhist->Write();
  muGRAVSMhist->Write();
  muGRAVGRAVhist->Write();
  testStatSMH->Write();
  testStatGRAVH->Write();
  testStatOBSH->Write();

  outFile->Close();
  inFile->Close();

  return 0;
}

