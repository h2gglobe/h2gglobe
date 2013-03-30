#include <algorithm>

void diphomva_public(bool passMVAcut=false, bool equalArea=true) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString txt_passcut = "";
  if (passMVAcut) txt_passcut = "_passcut";

  TFile *file_data = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_zee_moriond_preapproval_v2_runABCD_idmvaSigmaECorr.root");
  file_data->cd();

  bdtout_cat0_DataClone = (TH1*)bdtout_cat0_Data->Clone();

  TFile *file_MC = TFile::Open("zee_envelopes_runABCD_idmvaSigmaECorr.root");
  file_MC->cd();

  bdtout_cat0_DYJetsToLL->Rebin(2);
  bdtout_cat0_DataClone->Rebin(2);
  bdtout_cat0_DYJetsToLL_top->Rebin(2);
  bdtout_cat0_DYJetsToLL_bottom->Rebin(2);
  bdtout_cat0_DataClone->GetXaxis()->SetTitle("diphoton BDT output");
  bdtout_cat0_DataClone->GetYaxis()->SetTitle("");
  bdtout_cat0_DataClone->GetXaxis()->SetTitleSize(0.045);
  bdtout_cat0_DYJetsToLL->SetFillColor(38);
  bdtout_cat0_DYJetsToLL->SetLineColor(1);
  bdtout_cat0_DYJetsToLL_top->SetLineColor(2);
  bdtout_cat0_DYJetsToLL_top->SetMarkerStyle(0);
  bdtout_cat0_DYJetsToLL_bottom->SetLineColor(2);
  bdtout_cat0_DataClone->SetMarkerStyle(20);
  bdtout_cat0_DataClone->SetMarkerSize(.7);
  bdtout_cat0_DataClone->SetLineColor(1);

  if (passMVAcut) bdtout_cat0_DataClone->GetXaxis()->SetRangeUser(-0.05,1.);

  bdtout_cat0_syst = (TH1*)bdtout_cat0_DYJetsToLL_top->Clone();
  bdtout_cat0_syst->Add(bdtout_cat0_DYJetsToLL_bottom);
  bdtout_cat0_syst->Scale(0.5);
  for (int i=1; i<bdtout_cat0_syst->GetNbinsX()+1; i++) {
    float up = bdtout_cat0_DYJetsToLL_top->GetBinContent(i);
    float down = bdtout_cat0_DYJetsToLL_bottom->GetBinContent(i);
    bdtout_cat0_syst->SetBinError(i,fabs(up-down)/2.);
  }

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.04);

  TLegend *leg;
  leg = new TLegend(.15,.65,.4,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.03);
  leg->AddEntry(bdtout_cat0_DataClone,"Data");
  leg->AddEntry(bdtout_cat0_DYJetsToLL,"Z#rightarrow e^{+}e^{-} MC","F");
  leg->AddEntry(bdtout_cat0_syst,"MC systematic uncertainty","F");

  float sf_presel = bdtout_cat0_DataClone->Integral()/bdtout_cat0_DYJetsToLL->Integral();

  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output",1200,800);

  float sf;
  if (passMVAcut) {
    sf = bdtout_cat0_DataClone->Integral(48,100)/bdtout_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat0_DataClone->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  }
  bdtout_cat0_DYJetsToLL->Scale(sf);
  bdtout_cat0_DYJetsToLL_top->Scale(sf);
  bdtout_cat0_DYJetsToLL_bottom->Scale(sf);
  bdtout_cat0_syst->Scale(sf);
  bdtout_cat0_DataClone->Draw("e");
  leg->Draw();
  txt->DrawLatex(0.45,0.8,"#scale[0.8]{#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

  bdtout_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_cat0_syst->SetFillStyle(3013);
  bdtout_cat0_syst->SetFillColor(2);
  bdtout_cat0_syst->Draw("e2,same");
  bdtout_cat0_DYJetsToLL_top->Draw("hist,same");
  bdtout_cat0_DYJetsToLL_bottom->Draw("hist,same");
  bdtout_cat0_DataClone->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout->SaveAs("diphomva.png");
  c_bdtout->SaveAs("diphomva.pdf");

}
