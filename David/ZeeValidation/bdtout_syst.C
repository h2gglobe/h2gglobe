#include <algorithm>

void bdtout_syst(bool passMVAcut=false, bool equalArea=true) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString txt_passcut = "";
  if (passMVAcut) txt_passcut = "_passcut";

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  //TFile *file = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_zee_moriond_preapproval_phoPD_MCtriggers_ptreweight_noEcalIso.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  bdtout_cat0_DYJetsToLL->Rebin(2);
  bdtout_cat1_DYJetsToLL->Rebin(2);
  bdtout_cat2_DYJetsToLL->Rebin(2);
  bdtout_cat3_DYJetsToLL->Rebin(2);
  bdtout_cat4_DYJetsToLL->Rebin(2);
  bdtoutEB_cat0_DYJetsToLL->Rebin(2);
  bdtoutEBEE_cat0_DYJetsToLL->Rebin(2);
  bdtoutEE_cat0_DYJetsToLL->Rebin(2);
  bdtout_cat0_Data->Rebin(2);
  bdtout_cat1_Data->Rebin(2);
  bdtout_cat2_Data->Rebin(2);
  bdtout_cat3_Data->Rebin(2);
  bdtout_cat4_Data->Rebin(2);
  bdtoutEB_cat0_Data->Rebin(2);
  bdtoutEBEE_cat0_Data->Rebin(2);
  bdtoutEE_cat0_Data->Rebin(2);

  bdtout_up_cat0_DYJetsToLL->Rebin(2);
  bdtout_up_cat1_DYJetsToLL->Rebin(2);
  bdtout_up_cat2_DYJetsToLL->Rebin(2);
  bdtout_up_cat3_DYJetsToLL->Rebin(2);
  bdtout_up_cat4_DYJetsToLL->Rebin(2);
  bdtoutEB_up_cat0_DYJetsToLL->Rebin(2);
  bdtoutEBEE_up_cat0_DYJetsToLL->Rebin(2);
  bdtoutEE_up_cat0_DYJetsToLL->Rebin(2);

  bdtout_down_cat0_DYJetsToLL->Rebin(2);
  bdtout_down_cat1_DYJetsToLL->Rebin(2);
  bdtout_down_cat2_DYJetsToLL->Rebin(2);
  bdtout_down_cat3_DYJetsToLL->Rebin(2);
  bdtout_down_cat4_DYJetsToLL->Rebin(2);
  bdtoutEB_down_cat0_DYJetsToLL->Rebin(2);
  bdtoutEBEE_down_cat0_DYJetsToLL->Rebin(2);
  bdtoutEE_down_cat0_DYJetsToLL->Rebin(2);

  bdtout_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output");
  bdtout_cat1_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat0");
  bdtout_cat2_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat1");
  bdtout_cat3_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat2");
  bdtout_cat4_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat3");
  bdtoutEB_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EB)");
  bdtoutEBEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (one EB, one EE)");
  bdtoutEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EE)");

  bdtout_cat0_Data->GetYaxis()->SetTitle("");
  bdtout_cat1_Data->GetYaxis()->SetTitle("");
  bdtout_cat2_Data->GetYaxis()->SetTitle("");
  bdtout_cat3_Data->GetYaxis()->SetTitle("");
  bdtout_cat4_Data->GetYaxis()->SetTitle("");
  bdtoutEB_cat0_Data->GetYaxis()->SetTitle("");
  bdtoutEBEE_cat0_Data->GetYaxis()->SetTitle("");
  bdtoutEE_cat0_Data->GetYaxis()->SetTitle("");

  bdtout_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEB_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEBEE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEE_cat0_Data->GetXaxis()->SetTitleSize(0.05);

  bdtout_cat0_DYJetsToLL->SetFillColor(38);
  bdtout_cat1_DYJetsToLL->SetFillColor(38);
  bdtout_cat2_DYJetsToLL->SetFillColor(38);
  bdtout_cat3_DYJetsToLL->SetFillColor(38);
  bdtout_cat4_DYJetsToLL->SetFillColor(38);
  bdtoutEB_cat0_DYJetsToLL->SetFillColor(38);
  bdtoutEBEE_cat0_DYJetsToLL->SetFillColor(38);
  bdtoutEE_cat0_DYJetsToLL->SetFillColor(38);

  bdtout_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEB_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEBEE_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEE_cat0_DYJetsToLL->SetLineColor(1);

  bdtout_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat1_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat2_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat3_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat4_DYJetsToLL->SetLineColor(2);
  bdtoutEB_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEBEE_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEE_up_cat0_DYJetsToLL->SetLineColor(2);

  bdtout_down_cat0_DYJetsToLL->SetLineColor(2);
  bdtout_down_cat1_DYJetsToLL->SetLineColor(2);
  bdtout_down_cat2_DYJetsToLL->SetLineColor(2);
  bdtout_down_cat3_DYJetsToLL->SetLineColor(2);
  bdtout_down_cat4_DYJetsToLL->SetLineColor(2);
  bdtoutEB_down_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEBEE_down_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEE_down_cat0_DYJetsToLL->SetLineColor(2);

  bdtout_cat0_Data->SetMarkerStyle(20);
  bdtout_cat1_Data->SetMarkerStyle(20);
  bdtout_cat2_Data->SetMarkerStyle(20);
  bdtout_cat3_Data->SetMarkerStyle(20);
  bdtout_cat4_Data->SetMarkerStyle(20);
  bdtoutEB_cat0_Data->SetMarkerStyle(20);
  bdtoutEBEE_cat0_Data->SetMarkerStyle(20);
  bdtoutEE_cat0_Data->SetMarkerStyle(20);

  bdtout_cat0_Data->SetMarkerSize(0.4);
  bdtout_cat1_Data->SetMarkerSize(0.4);
  bdtout_cat2_Data->SetMarkerSize(0.4);
  bdtout_cat3_Data->SetMarkerSize(0.4);
  bdtout_cat4_Data->SetMarkerSize(0.4);
  bdtoutEB_cat0_Data->SetMarkerSize(0.4);
  bdtoutEBEE_cat0_Data->SetMarkerSize(0.4);
  bdtoutEE_cat0_Data->SetMarkerSize(0.4);

  bdtout_cat0_Data->SetLineColor(1);
  bdtout_cat1_Data->SetLineColor(1);
  bdtout_cat2_Data->SetLineColor(1);
  bdtout_cat3_Data->SetLineColor(1);
  bdtout_cat4_Data->SetLineColor(1);
  bdtoutEB_cat0_Data->SetLineColor(1);
  bdtoutEBEE_cat0_Data->SetLineColor(1);
  bdtoutEE_cat0_Data->SetLineColor(1);

  if (passMVAcut) {
    bdtout_cat0_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtout_cat1_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtout_cat2_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtout_cat3_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtout_cat4_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtoutEB_cat0_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtoutEBEE_cat0_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    bdtoutEE_cat0_Data->GetXaxis()->SetRangeUser(-0.05,1.);
  }

  bdtout_cat0_syst = (TH1*)bdtout_up_cat0_DYJetsToLL->Clone();
  bdtout_cat0_syst->Add(bdtout_down_cat0_DYJetsToLL);
  bdtout_cat0_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat0_syst->GetNbinsX()+1; i++) {
    float up = bdtout_up_cat0_DYJetsToLL->GetBinContent(i);
    float down = bdtout_down_cat0_DYJetsToLL->GetBinContent(i);
    bdtout_cat0_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtout_cat1_syst = (TH1*)bdtout_up_cat1_DYJetsToLL->Clone();
  bdtout_cat1_syst->Add(bdtout_down_cat1_DYJetsToLL);
  bdtout_cat1_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat1_syst->GetNbinsX()+1; i++) {
    float up = bdtout_up_cat1_DYJetsToLL->GetBinContent(i);
    float down = bdtout_down_cat1_DYJetsToLL->GetBinContent(i);
    bdtout_cat1_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtout_cat2_syst = (TH1*)bdtout_up_cat2_DYJetsToLL->Clone();
  bdtout_cat2_syst->Add(bdtout_down_cat2_DYJetsToLL);
  bdtout_cat2_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat2_syst->GetNbinsX()+1; i++) {
    float up = bdtout_up_cat2_DYJetsToLL->GetBinContent(i);
    float down = bdtout_down_cat2_DYJetsToLL->GetBinContent(i);
    bdtout_cat2_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtout_cat3_syst = (TH1*)bdtout_up_cat3_DYJetsToLL->Clone();
  bdtout_cat3_syst->Add(bdtout_down_cat3_DYJetsToLL);
  bdtout_cat3_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat3_syst->GetNbinsX()+1; i++) {
    float up = bdtout_up_cat3_DYJetsToLL->GetBinContent(i);
    float down = bdtout_down_cat3_DYJetsToLL->GetBinContent(i);
    bdtout_cat3_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtout_cat4_syst = (TH1*)bdtout_up_cat4_DYJetsToLL->Clone();
  bdtout_cat4_syst->Add(bdtout_down_cat4_DYJetsToLL);
  bdtout_cat4_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat4_syst->GetNbinsX()+1; i++) {
    float up = bdtout_up_cat4_DYJetsToLL->GetBinContent(i);
    float down = bdtout_down_cat4_DYJetsToLL->GetBinContent(i);
    bdtout_cat4_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtoutEB_syst = (TH1*)bdtoutEB_up_cat0_DYJetsToLL->Clone();
  bdtoutEB_syst->Add(bdtoutEB_down_cat0_DYJetsToLL);
  bdtoutEB_syst->Scale(0.5);
  for (int i=1; i<bdtoutEB_syst->GetNbinsX()+1; i++) {
    float up = bdtoutEB_up_cat0_DYJetsToLL->GetBinContent(i);
    float down = bdtoutEB_down_cat0_DYJetsToLL->GetBinContent(i);
    bdtoutEB_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtoutEBEE_syst = (TH1*)bdtoutEBEE_up_cat0_DYJetsToLL->Clone();
  bdtoutEBEE_syst->Add(bdtoutEBEE_down_cat0_DYJetsToLL);
  bdtoutEBEE_syst->Scale(0.5);
  for (int i=1; i<bdtoutEBEE_syst->GetNbinsX()+1; i++) {
    float up = bdtoutEBEE_up_cat0_DYJetsToLL->GetBinContent(i);
    float down = bdtoutEBEE_down_cat0_DYJetsToLL->GetBinContent(i);
    bdtoutEBEE_syst->SetBinError(i,fabs(up-down)/2.);
  }

  bdtoutEE_syst = (TH1*)bdtoutEE_up_cat0_DYJetsToLL->Clone();
  bdtoutEE_syst->Add(bdtoutEE_down_cat0_DYJetsToLL);
  bdtoutEE_syst->Scale(0.5);
  for (int i=1; i<bdtoutEE_syst->GetNbinsX()+1; i++) {
    float up = bdtoutEE_up_cat0_DYJetsToLL->GetBinContent(i);
    float down = bdtoutEE_down_cat0_DYJetsToLL->GetBinContent(i);
    bdtoutEE_syst->SetBinError(i,fabs(up-down)/2.);
  }

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg->AddEntry(bdtout_cat0_syst,"MC with idmva±0.01","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg2->AddEntry(bdtout_cat0_syst,"MC with idmva±0.01","F");
  //leg2->AddEntry(bdtout_Hgg90,"H#rightarrow#gamma#gamma, m_{H}=90 GeV");

  TLegend *leg3;
  leg3 = new TLegend(.45,.65,.87,.87);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(10);
  leg3->SetTextSize(.03);
  leg3->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg3->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg3->AddEntry(bdtout_cat0_syst,"MC with idmva±0.01","F");

  TLegend *leg4;
  leg4 = new TLegend(.2,.65,.52,.87);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(10);
  leg4->SetTextSize(.025);
  leg4->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg4->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg4->AddEntry(bdtout_cat0_syst,"MC with idmva±0.01","F");

  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();

  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output",2200,800);
  c_bdtout->Divide(4,2);

  c_bdtout->cd(1);
  float sf;
  if (passMVAcut) {
    sf = bdtout_cat0_Data->Integral(48,100)/bdtout_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  }
  bdtout_cat0_DYJetsToLL->Scale(sf);
  bdtout_up_cat0_DYJetsToLL->Scale(sf);
  bdtout_down_cat0_DYJetsToLL->Scale(sf);
  bdtout_cat0_syst->Scale(sf);
  bdtout_cat0_Data->Draw("e");
  leg2->Draw();
  bdtout_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_cat0_syst->SetFillStyle(3013);
  bdtout_cat0_syst->SetFillColor(2);
  bdtout_cat0_syst->Draw("e2,same");
  bdtout_up_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_down_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_cat0_Data->Draw("e,same");
  //bdtout_Hgg90->Scale(bdtout_cat0_Data->Integral()/bdtout_Hgg90->Integral());
  //bdtout_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout->cd(2);
  if (passMVAcut) {
    sf = bdtoutEB_cat0_Data->Integral(48,100)/bdtoutEB_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEB_cat0_Data->Integral()/bdtoutEB_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEB_cat0_DYJetsToLL->Scale(sf);
  bdtoutEB_up_cat0_DYJetsToLL->Scale(sf);
  bdtoutEB_down_cat0_DYJetsToLL->Scale(sf);
  bdtoutEB_syst->Scale(sf);
  bdtoutEB_cat0_Data->Draw("e");
  leg2->Draw();
  bdtoutEB_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEB_syst->SetFillStyle(3013);
  bdtoutEB_syst->SetFillColor(2);
  bdtoutEB_syst->Draw("e2,same");
  bdtoutEB_up_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEB_down_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEB_cat0_Data->Draw("e,same");
  //bdtoutEB_Hgg90->Scale(bdtoutEB_cat0_Data->Integral()/bdtoutEB_Hgg90->Integral());
  //bdtoutEB_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout->cd(3);
  if (passMVAcut) {
    sf = bdtoutEBEE_cat0_Data->Integral(48,100)/bdtoutEBEE_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEBEE_cat0_Data->Integral()/bdtoutEBEE_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEBEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEBEE_up_cat0_DYJetsToLL->Scale(sf);
  bdtoutEBEE_down_cat0_DYJetsToLL->Scale(sf);
  bdtoutEBEE_syst->Scale(sf);
  bdtoutEBEE_cat0_Data->Draw("e");
  leg->Draw();
  bdtoutEBEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEBEE_syst->SetFillStyle(3013);
  bdtoutEBEE_syst->SetFillColor(2);
  bdtoutEBEE_syst->Draw("e2,same");
  bdtoutEBEE_up_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEBEE_down_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEBEE_cat0_Data->Draw("e,same");
  //bdtoutEBEE_Hgg90->Scale(bdtoutEBEE_cat0_Data->Integral()/bdtoutEBEE_Hgg90->Integral());
  //bdtoutEBEE_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout->cd(4);
  if (passMVAcut) {
    sf = bdtoutEE_cat0_Data->Integral(48,100)/bdtoutEE_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEE_cat0_Data->Integral()/bdtoutEE_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEE_up_cat0_DYJetsToLL->Scale(sf);
  bdtoutEE_down_cat0_DYJetsToLL->Scale(sf);
  bdtoutEE_syst->Scale(sf);
  bdtoutEE_cat0_Data->Draw("e");
  leg->Draw();
  bdtoutEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEE_syst->SetFillStyle(3013);
  bdtoutEE_syst->SetFillColor(2);
  bdtoutEE_syst->Draw("e2,same");
  bdtoutEE_up_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEE_down_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEE_cat0_Data->Draw("e,same");
  //bdtoutEE_Hgg90->Scale(bdtoutEE_cat0_Data->Integral()/bdtoutEE_Hgg90->Integral());
  //bdtoutEE_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  TLine *line;
  if (passMVAcut) {
    line = new TLine(-0.06,1.,1.,1.);
  } else {
    line = new TLine(-1.,1.,1.,1.);
  }
  line->SetLineColor(4);
  line->SetLineWidth(2);

  c_bdtout->cd(5);
  gPad->SetGrid();
  plotRatio(bdtout_cat0_Data,bdtout_cat0_DYJetsToLL,bdtout_up_cat0_DYJetsToLL,bdtout_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  /*
  TFile *fout = new TFile("zee_bdtout_ratio.root","RECREATE");
  fout->cd();
  ratio_bdtout->Write("zee_bdtout_ratio");
  fout->Close();
  file->cd();
  */

  c_bdtout->cd(6);
  gPad->SetGrid();
  plotRatio(bdtoutEB_cat0_Data,bdtoutEB_cat0_DYJetsToLL,bdtoutEB_up_cat0_DYJetsToLL,bdtoutEB_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();


  c_bdtout->cd(7);
  gPad->SetGrid();
  plotRatio(bdtoutEBEE_cat0_Data,bdtoutEBEE_cat0_DYJetsToLL,bdtoutEBEE_up_cat0_DYJetsToLL,bdtoutEBEE_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->cd(8);
  gPad->SetGrid();
  plotRatio(bdtoutEE_cat0_Data,bdtoutEE_cat0_DYJetsToLL,bdtoutEE_up_cat0_DYJetsToLL,bdtoutEE_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->SaveAs("bdtout"+txt_passcut+"_syst.png");


  TCanvas *c_bdtout_basecat = new TCanvas("c_bdtout_basecat","BDT output in CiC categories",2200,800);
  c_bdtout_basecat->Divide(4,2);

  c_bdtout_basecat->cd(1);
  if (passMVAcut) {
    sf = bdtout_cat1_Data->Integral(48,100)/bdtout_cat1_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat1_Data->Integral()/bdtout_cat1_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat1_DYJetsToLL->Scale(sf);
  bdtout_up_cat1_DYJetsToLL->Scale(sf);
  bdtout_down_cat1_DYJetsToLL->Scale(sf);
  bdtout_cat1_syst->Scale(sf);
  bdtout_cat1_Data->Draw("e");
  leg2->Draw();
  bdtout_cat1_DYJetsToLL->Draw("hist,same");
  bdtout_cat1_syst->SetFillStyle(3013);
  bdtout_cat1_syst->SetFillColor(2);
  bdtout_cat1_syst->Draw("e2,same");
  bdtout_up_cat1_DYJetsToLL->Draw("hist,same");
  bdtout_down_cat1_DYJetsToLL->Draw("hist,same");
  bdtout_cat1_Data->Draw("e,same");
  //bdtout_cat1_Hgg90->Scale(bdtout_cat1_Data->Integral()/bdtout_cat1_Hgg90->Integral());
  //bdtout_cat1_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(2);
  if (passMVAcut) {
    sf = bdtout_cat2_Data->Integral(48,100)/bdtout_cat2_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat2_Data->Integral()/bdtout_cat2_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat2_DYJetsToLL->Scale(sf);
  bdtout_up_cat2_DYJetsToLL->Scale(sf);
  bdtout_down_cat2_DYJetsToLL->Scale(sf);
  bdtout_cat2_syst->Scale(sf);
  bdtout_cat2_Data->Draw("e");
  leg2->Draw();
  bdtout_cat2_DYJetsToLL->Draw("hist,same");
  bdtout_cat2_syst->SetFillStyle(3013);
  bdtout_cat2_syst->SetFillColor(2);
  bdtout_cat2_syst->Draw("e2,same");
  bdtout_up_cat2_DYJetsToLL->Draw("hist,same");
  bdtout_down_cat2_DYJetsToLL->Draw("hist,same");
  bdtout_cat2_Data->Draw("e,same");
  //bdtout_cat2_Hgg90->Scale(bdtout_cat2_Data->Integral()/bdtout_cat2_Hgg90->Integral());
  //bdtout_cat2_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(3);
  if (passMVAcut) {
    sf = bdtout_cat3_Data->Integral(48,100)/bdtout_cat3_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat3_Data->Integral()/bdtout_cat3_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat3_DYJetsToLL->Scale(sf);
  bdtout_up_cat3_DYJetsToLL->Scale(sf);
  bdtout_down_cat3_DYJetsToLL->Scale(sf);
  bdtout_cat3_syst->Scale(sf);
  bdtout_cat3_Data->Draw("e");
  leg->Draw();
  bdtout_cat3_DYJetsToLL->Draw("hist,same");
  bdtout_cat3_syst->SetFillStyle(3013);
  bdtout_cat3_syst->SetFillColor(2);
  bdtout_cat3_syst->Draw("e2,same");
  bdtout_up_cat3_DYJetsToLL->Draw("hist,same");
  bdtout_down_cat3_DYJetsToLL->Draw("hist,same");
  bdtout_cat3_Data->Draw("e,same");
  //bdtout_cat3_Hgg90->Scale(bdtout_cat3_Data->Integral()/bdtout_cat3_Hgg90->Integral());
  //bdtout_cat3_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(4);
  if (passMVAcut) {
    sf = bdtout_cat4_Data->Integral(48,100)/bdtout_cat4_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat4_Data->Integral()/bdtout_cat4_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat4_DYJetsToLL->Scale(sf);
  bdtout_up_cat4_DYJetsToLL->Scale(sf);
  bdtout_down_cat4_DYJetsToLL->Scale(sf);
  bdtout_cat4_syst->Scale(sf);
  bdtout_cat4_Data->Draw("e");
  leg->Draw();
  bdtout_cat4_DYJetsToLL->Draw("hist,same");
  bdtout_cat4_syst->SetFillStyle(3013);
  bdtout_cat4_syst->SetFillColor(2);
  bdtout_cat4_syst->Draw("e2,same");
  bdtout_up_cat4_DYJetsToLL->Draw("hist,same");
  bdtout_down_cat4_DYJetsToLL->Draw("hist,same");
  bdtout_cat4_Data->Draw("e,same");
  //bdtout_cat4_Hgg90->Scale(bdtout_cat4_Data->Integral()/bdtout_cat4_Hgg90->Integral());
  //bdtout_cat4_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(5);
  gPad->SetGrid();
  plotRatio(bdtout_cat1_Data,bdtout_cat1_DYJetsToLL,bdtout_up_cat1_DYJetsToLL,bdtout_down_cat1_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(6);
  gPad->SetGrid();
  plotRatio(bdtout_cat2_Data,bdtout_cat2_DYJetsToLL,bdtout_up_cat2_DYJetsToLL,bdtout_down_cat2_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(7);
  gPad->SetGrid();
  plotRatio(bdtout_cat3_Data,bdtout_cat3_DYJetsToLL,bdtout_up_cat3_DYJetsToLL,bdtout_down_cat3_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(8);
  gPad->SetGrid();
  plotRatio(bdtout_cat4_Data,bdtout_cat4_DYJetsToLL,bdtout_up_cat4_DYJetsToLL,bdtout_down_cat4_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->SaveAs("bdtout_basecat"+txt_passcut+"_syst.png");


}

void plotRatio(TH1* bdtout_Data, TH1* bdtout_DYJetsToLL, TH1* bdtout_up_DYJetsToLL, TH1* bdtout_down_DYJetsToLL, bool passMVAcut)
{

  ratio_bdtout = (TH1*)bdtout_Data->Clone();
  ratio_bdtout_up = (TH1*)bdtout_up_DYJetsToLL->Clone();
  ratio_bdtout_down = (TH1*)bdtout_down_DYJetsToLL->Clone();
  ratio_bdtout->Divide(bdtout_DYJetsToLL);
  ratio_bdtout_up->Divide(bdtout_DYJetsToLL);
  ratio_bdtout_down->Divide(bdtout_DYJetsToLL);
  ratio_bdtout_up->SetLineColor(2);
  ratio_bdtout_down->SetLineColor(2);

  ratio_bdtout_syst = (TH1*)ratio_bdtout_up->Clone();
  ratio_bdtout_syst->Add(ratio_bdtout_down);
  ratio_bdtout_syst->Scale(0.5);
  TBox *ratio_bdtout_syst_box[400];
  for (int i=1; i<ratio_bdtout_syst->GetNbinsX()+1; i++) {
    float up = max(ratio_bdtout_up->GetBinContent(i),ratio_bdtout_down->GetBinContent(i));
    float down = min(ratio_bdtout_up->GetBinContent(i),ratio_bdtout_down->GetBinContent(i));
    ratio_bdtout_syst->SetBinError(i,fabs(up-down)/2.);
    ratio_bdtout_syst_box[i-1] = new TBox(ratio_bdtout_up->GetXaxis()->GetBinLowEdge(i),max(0.5,down),ratio_bdtout_up->GetXaxis()->GetBinLowEdge(i+1),min(2.,up));
    ratio_bdtout_syst_box[i-1]->SetFillStyle(3013);
    ratio_bdtout_syst_box[i-1]->SetFillColor(2);
  }

  ratio_bdtout_syst->SetFillStyle(3013);
  ratio_bdtout_syst->SetFillColor(2);
  ratio_bdtout_syst->SetMarkerStyle(0);
  ratio_bdtout_syst->SetMaximum(1.8);
  ratio_bdtout_syst->SetMinimum(0.2);
  ratio_bdtout_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_bdtout_syst->GetYaxis()->SetTitleSize(0.05);
  int istart=1;
  if (passMVAcut) {
    istart=48;
    ratio_bdtout_syst->GetXaxis()->SetRangeUser(-0.05,1.);
  }
  ratio_bdtout_syst->Draw("e2");
  for (int i=istart; i<ratio_bdtout_syst->GetNbinsX()+1; i++) {
    if (ratio_bdtout_syst_box[i-1]->GetY2()>0.5&&ratio_bdtout_syst_box[i-1]->GetY1()<2.) ratio_bdtout_syst_box[i-1]->Draw("hist,same");
  }
  ratio_bdtout_up->Draw("hist,same");
  ratio_bdtout_down->Draw("hist,same");
  ratio_bdtout->Draw("e,same");

}
