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
  file->cd();

  bdtout_cat0_Data->Rebin(2);
  bdtout_cat1_Data->Rebin(2);
  bdtout_cat2_Data->Rebin(2);
  bdtout_cat3_Data->Rebin(2);
  bdtout_cat4_Data->Rebin(2);
  bdtoutEB_cat0_Data->Rebin(2);
  bdtoutEBEE_cat0_Data->Rebin(2);
  bdtoutEE_cat0_Data->Rebin(2);

  bdtout_cat0_DYJetsToLL->Rebin(2);
  bdtout_cat1_DYJetsToLL->Rebin(2);
  bdtout_cat2_DYJetsToLL->Rebin(2);
  bdtout_cat3_DYJetsToLL->Rebin(2);
  bdtout_cat4_DYJetsToLL->Rebin(2);
  bdtoutEB_cat0_DYJetsToLL->Rebin(2);
  bdtoutEBEE_cat0_DYJetsToLL->Rebin(2);
  bdtoutEE_cat0_DYJetsToLL->Rebin(2);

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

  bdtout_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat1_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat2_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat3_DYJetsToLL->SetLineColor(2);
  bdtout_up_cat4_DYJetsToLL->SetLineColor(2);
  bdtoutEB_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEBEE_up_cat0_DYJetsToLL->SetLineColor(2);
  bdtoutEE_up_cat0_DYJetsToLL->SetLineColor(2);

  bdtout_up_cat0_DYJetsToLL->SetMarkerStyle(0);
  bdtout_up_cat1_DYJetsToLL->SetMarkerStyle(0);
  bdtout_up_cat2_DYJetsToLL->SetMarkerStyle(0);
  bdtout_up_cat3_DYJetsToLL->SetMarkerStyle(0);
  bdtout_up_cat4_DYJetsToLL->SetMarkerStyle(0);
  bdtoutEB_up_cat0_DYJetsToLL->SetMarkerStyle(0);
  bdtoutEBEE_up_cat0_DYJetsToLL->SetMarkerStyle(0);
  bdtoutEE_up_cat0_DYJetsToLL->SetMarkerStyle(0);

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

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  hist_syst = (TH1*)bdtout_up_cat0_DYJetsToLL->Clone();
  hist_syst->SetFillStyle(3013);
  hist_syst->SetFillColor(2);

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg->AddEntry(hist_syst,"MC with idmva±0.01","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(bdtout_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(bdtout_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  leg2->AddEntry(hist_syst,"MC with idmva±0.01","F");

  TLine *line;
  if (passMVAcut) {
    line = new TLine(-0.06,1.,1.,1.);
  } else {
    line = new TLine(-1.,1.,1.,1.);
  }
  line->SetLineColor(4);
  line->SetLineWidth(2);

  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();


  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output",2200,800);
  c_bdtout->Divide(4,2);

  c_bdtout->cd(1);
  plotDataMC(bdtout_cat0_Data, bdtout_cat0_DYJetsToLL, bdtout_up_cat0_DYJetsToLL, bdtout_down_cat0_DYJetsToLL, passMVAcut, equalArea, leg2);

  c_bdtout->cd(2);
  plotDataMC(bdtoutEB_cat0_Data, bdtoutEB_cat0_DYJetsToLL, bdtoutEB_up_cat0_DYJetsToLL, bdtoutEB_down_cat0_DYJetsToLL, passMVAcut, equalArea, leg2);

  c_bdtout->cd(3);
  plotDataMC(bdtoutEBEE_cat0_Data, bdtoutEBEE_cat0_DYJetsToLL, bdtoutEBEE_up_cat0_DYJetsToLL, bdtoutEBEE_down_cat0_DYJetsToLL, passMVAcut, equalArea, leg);

  c_bdtout->cd(4);
  plotDataMC(bdtoutEE_cat0_Data, bdtoutEE_cat0_DYJetsToLL, bdtoutEE_up_cat0_DYJetsToLL, bdtoutEE_down_cat0_DYJetsToLL, passMVAcut, equalArea, leg);

  c_bdtout->cd(5);
  plotRatio(bdtout_cat0_Data,bdtout_cat0_DYJetsToLL,bdtout_up_cat0_DYJetsToLL,bdtout_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->cd(6);
  plotRatio(bdtoutEB_cat0_Data,bdtoutEB_cat0_DYJetsToLL,bdtoutEB_up_cat0_DYJetsToLL,bdtoutEB_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->cd(7);
  plotRatio(bdtoutEBEE_cat0_Data,bdtoutEBEE_cat0_DYJetsToLL,bdtoutEBEE_up_cat0_DYJetsToLL,bdtoutEBEE_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->cd(8);
  plotRatio(bdtoutEE_cat0_Data,bdtoutEE_cat0_DYJetsToLL,bdtoutEE_up_cat0_DYJetsToLL,bdtoutEE_down_cat0_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout->SaveAs("bdtout"+txt_passcut+"_syst.png");


  TCanvas *c_bdtout_basecat = new TCanvas("c_bdtout_basecat","BDT output in CiC categories",2200,800);
  c_bdtout_basecat->Divide(4,2);

  c_bdtout_basecat->cd(1);
  plotDataMC(bdtout_cat1_Data, bdtout_cat1_DYJetsToLL, bdtout_up_cat1_DYJetsToLL, bdtout_down_cat1_DYJetsToLL, passMVAcut, equalArea, leg2);

  c_bdtout_basecat->cd(2);
  plotDataMC(bdtout_cat2_Data, bdtout_cat2_DYJetsToLL, bdtout_up_cat2_DYJetsToLL, bdtout_down_cat2_DYJetsToLL, passMVAcut, equalArea, leg2);

  c_bdtout_basecat->cd(3);
  plotDataMC(bdtout_cat3_Data, bdtout_cat3_DYJetsToLL, bdtout_up_cat3_DYJetsToLL, bdtout_down_cat3_DYJetsToLL, passMVAcut, equalArea, leg);

  c_bdtout_basecat->cd(4);
  plotDataMC(bdtout_cat4_Data, bdtout_cat4_DYJetsToLL, bdtout_up_cat4_DYJetsToLL, bdtout_down_cat4_DYJetsToLL, passMVAcut, equalArea, leg);

  c_bdtout_basecat->cd(5);
  plotRatio(bdtout_cat1_Data,bdtout_cat1_DYJetsToLL,bdtout_up_cat1_DYJetsToLL,bdtout_down_cat1_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(6);
  plotRatio(bdtout_cat2_Data,bdtout_cat2_DYJetsToLL,bdtout_up_cat2_DYJetsToLL,bdtout_down_cat2_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(7);
  plotRatio(bdtout_cat3_Data,bdtout_cat3_DYJetsToLL,bdtout_up_cat3_DYJetsToLL,bdtout_down_cat3_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->cd(8);
  plotRatio(bdtout_cat4_Data,bdtout_cat4_DYJetsToLL,bdtout_up_cat4_DYJetsToLL,bdtout_down_cat4_DYJetsToLL, passMVAcut);
  line->Draw();

  c_bdtout_basecat->SaveAs("bdtout_basecat"+txt_passcut+"_syst.png");


}


void plotDataMC(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, bool passMVAcut, bool equalArea, TLegend* leg)
{

  hist_syst = (TH1*)hist_MC_up->Clone();
  hist_syst->Add(hist_MC_down);
  hist_syst->Scale(0.5);
  for (int i=1; i<hist_syst->GetNbinsX()+1; i++) {
    float up = hist_MC_up->GetBinContent(i);
    float down = hist_MC_down->GetBinContent(i);
    hist_syst->SetBinError(i,fabs(up-down)/2.);
  }

  float sf;
  if (passMVAcut) {
    sf = hist_Data->Integral(48,100)/hist_MC->Integral(48,100);
  } else {
    sf = hist_Data->Integral()/hist_MC->Integral();
    if (!equalArea) sf = sf_presel;
  }
  hist_MC->Scale(sf);
  hist_MC_up->Scale(sf);
  hist_MC_down->Scale(sf);
  hist_syst->Scale(sf);
  hist_Data->Draw("e");
  leg->Draw();
  hist_MC->Draw("hist,same");
  hist_syst->SetFillStyle(3013);
  hist_syst->SetFillColor(2);
  hist_syst->Draw("e2,same");
  hist_MC_up->Draw("hist,same");
  hist_MC_down->Draw("hist,same");
  hist_Data->Draw("e,same");
  gPad->RedrawAxis();

}


void plotRatio(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, bool passMVAcut)
{

  gPad->SetGrid();

  ratio = (TH1*)hist_Data->Clone();
  ratio_up = (TH1*)hist_MC_up->Clone();
  ratio_down = (TH1*)hist_MC_down->Clone();
  ratio->Divide(hist_MC);
  ratio_up->Divide(hist_MC);
  ratio_down->Divide(hist_MC);
  ratio_up->SetLineColor(2);
  ratio_down->SetLineColor(2);

  ratio_syst = (TH1*)ratio_up->Clone();
  ratio_syst->Add(ratio_down);
  ratio_syst->Scale(0.5);
  TBox *ratio_syst_box[400];
  for (int i=1; i<ratio_syst->GetNbinsX()+1; i++) {
    float up = max(ratio_up->GetBinContent(i),ratio_down->GetBinContent(i));
    float down = min(ratio_up->GetBinContent(i),ratio_down->GetBinContent(i));
    ratio_syst->SetBinError(i,fabs(up-down)/2.);
    ratio_syst_box[i-1] = new TBox(ratio_up->GetXaxis()->GetBinLowEdge(i),max(0.5,down),ratio_up->GetXaxis()->GetBinLowEdge(i+1),min(2.,up));
    ratio_syst_box[i-1]->SetFillStyle(3013);
    ratio_syst_box[i-1]->SetFillColor(2);
  }

  ratio_syst->SetFillStyle(3013);
  ratio_syst->SetFillColor(2);
  ratio_syst->SetMarkerStyle(0);
  ratio_syst->SetMaximum(1.8);
  ratio_syst->SetMinimum(0.2);
  ratio_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_syst->GetYaxis()->SetTitleSize(0.05);
  int istart=1;
  if (passMVAcut) {
    istart=48;
    ratio_syst->GetXaxis()->SetRangeUser(-0.05,1.);
  }
  ratio_syst->Draw("e2");
  for (int i=istart; i<ratio_syst->GetNbinsX()+1; i++) {
    if (ratio_syst_box[i-1]->GetY2()>0.5&&ratio_syst_box[i-1]->GetY1()<2.) ratio_syst_box[i-1]->Draw("hist,same");
  }
  ratio_up->Draw("hist,same");
  ratio_down->Draw("hist,same");
  ratio->Draw("e,same");

}
