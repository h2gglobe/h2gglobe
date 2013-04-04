#include <algorithm>

void idmva() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  phoidMva_EB_Data = (TH1*)pho1_phoidMva_EB_cat0_Data->Clone();
  phoidMva_EB_Data->Add(pho2_phoidMva_EB_cat0_Data);
  phoidMva_EE_Data = (TH1*)pho1_phoidMva_EE_cat0_Data->Clone();
  phoidMva_EE_Data->Add(pho2_phoidMva_EE_cat0_Data);

  phoidMva_EB_DYJetsToLL = (TH1*)pho1_phoidMva_EB_cat0_DYJetsToLL->Clone();
  phoidMva_EB_DYJetsToLL->Add(pho2_phoidMva_EB_cat0_DYJetsToLL);
  phoidMva_EE_DYJetsToLL = (TH1*)pho1_phoidMva_EE_cat0_DYJetsToLL->Clone();
  phoidMva_EE_DYJetsToLL->Add(pho2_phoidMva_EE_cat0_DYJetsToLL);

  phoidMva_EB_up_DYJetsToLL = (TH1*)pho1_phoidMva_EB_up_cat0_DYJetsToLL->Clone();
  phoidMva_EB_up_DYJetsToLL->Add(pho2_phoidMva_EB_up_cat0_DYJetsToLL);
  phoidMva_EE_up_DYJetsToLL = (TH1*)pho1_phoidMva_EE_up_cat0_DYJetsToLL->Clone();
  phoidMva_EE_up_DYJetsToLL->Add(pho2_phoidMva_EE_up_cat0_DYJetsToLL);

  phoidMva_EB_down_DYJetsToLL = (TH1*)pho1_phoidMva_EB_down_cat0_DYJetsToLL->Clone();
  phoidMva_EB_down_DYJetsToLL->Add(pho2_phoidMva_EB_down_cat0_DYJetsToLL);
  phoidMva_EE_down_DYJetsToLL = (TH1*)pho1_phoidMva_EE_down_cat0_DYJetsToLL->Clone();
  phoidMva_EE_down_DYJetsToLL->Add(pho2_phoidMva_EE_down_cat0_DYJetsToLL);

  phoidMva_EB_up_DYJetsToLL->SetLineColor(2);
  phoidMva_EB_down_DYJetsToLL->SetLineColor(2);
  phoidMva_EE_up_DYJetsToLL->SetLineColor(2);
  phoidMva_EE_down_DYJetsToLL->SetLineColor(2);

  phoidMva_EB_Data->GetYaxis()->SetTitle("");
  phoidMva_EE_Data->GetYaxis()->SetTitle("");
  phoidMva_EB_Data->GetXaxis()->SetTitle("photon ID MVA output (EB)");
  phoidMva_EE_Data->GetXaxis()->SetTitle("photon ID MVA output (EE)");
  phoidMva_EB_Data->GetXaxis()->SetTitleSize(0.05);
  phoidMva_EE_Data->GetXaxis()->SetTitleSize(0.05);
  phoidMva_EB_DYJetsToLL->SetFillColor(38);
  phoidMva_EE_DYJetsToLL->SetFillColor(38);
  phoidMva_EB_Data->SetMarkerStyle(20);
  phoidMva_EE_Data->SetMarkerStyle(20);
  phoidMva_EB_Data->SetLineColor(1);
  phoidMva_EE_Data->SetLineColor(1);
  phoidMva_EB_Data->SetMarkerSize(.5);
  phoidMva_EE_Data->SetMarkerSize(.5);

  phoidMva_EB_syst = (TH1*)phoidMva_EB_up_DYJetsToLL->Clone();
  phoidMva_EB_syst->Add(phoidMva_EB_down_DYJetsToLL);
  phoidMva_EB_syst->Scale(0.5);
  for (int i=1; i<phoidMva_EB_syst->GetNbinsX()+1; i++) {
    float up = phoidMva_EB_up_DYJetsToLL->GetBinContent(i);
    float down = phoidMva_EB_down_DYJetsToLL->GetBinContent(i);
    phoidMva_EB_syst->SetBinError(i,fabs(up-down)/2.);
  }

  phoidMva_EE_syst = (TH1*)phoidMva_EE_up_DYJetsToLL->Clone();
  phoidMva_EE_syst->Add(phoidMva_EE_down_DYJetsToLL);
  phoidMva_EE_syst->Scale(0.5);
  for (int i=1; i<phoidMva_EE_syst->GetNbinsX()+1; i++) {
    float up = phoidMva_EE_up_DYJetsToLL->GetBinContent(i);
    float down = phoidMva_EE_down_DYJetsToLL->GetBinContent(i);
    phoidMva_EE_syst->SetBinError(i,fabs(up-down)/2.);
  }

  TLegend *leg2;
  leg2 = new TLegend(.16,.55,.45,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(phoidMva_EB_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(phoidMva_EB_DYJetsToLL,"DYJetsToLL MC","F");
  leg2->AddEntry(phoidMva_EB_syst,"MC ±0.01","F");

  TCanvas *c_idmva = new TCanvas("c_idmva","Photon ID MVA output",1600,900);
  c_idmva->Divide(2,2);

  c_idmva->cd(1);
  float sf = phoidMva_EB_Data->Integral()/phoidMva_EB_DYJetsToLL->Integral();
  phoidMva_EB_DYJetsToLL->Scale(sf);
  phoidMva_EB_down_DYJetsToLL->Scale(sf);
  phoidMva_EB_up_DYJetsToLL->Scale(sf);
  phoidMva_EB_syst->Scale(sf);
  phoidMva_EB_DYJetsToLL_line = (TH1*)phoidMva_EB_DYJetsToLL->Clone();
  phoidMva_EB_DYJetsToLL_line->SetFillColor(0);
  phoidMva_EB_Data->GetXaxis()->SetRangeUser(-0.5,0.5);
  float max = phoidMva_EB_Data->GetMaximum();
  phoidMva_EB_Data->SetMaximum(max*1.1);
  phoidMva_EB_Data->Draw("e");
  phoidMva_EB_DYJetsToLL->Draw("hist,same");
  phoidMva_EB_syst->SetFillStyle(3013);
  phoidMva_EB_syst->SetFillColor(2);
  phoidMva_EB_syst->Draw("e2,same");
  phoidMva_EB_up_DYJetsToLL->Draw("hist,same");
  phoidMva_EB_down_DYJetsToLL->Draw("hist,same");
  phoidMva_EB_DYJetsToLL_line->Draw("hist,same");
  phoidMva_EB_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_idmva->cd(2);
  sf = phoidMva_EE_Data->Integral()/phoidMva_EE_DYJetsToLL->Integral();
  phoidMva_EE_DYJetsToLL->Scale(sf);
  phoidMva_EE_down_DYJetsToLL->Scale(sf);
  phoidMva_EE_up_DYJetsToLL->Scale(sf);
  phoidMva_EE_syst->Scale(sf);
  phoidMva_EE_DYJetsToLL_line = (TH1*)phoidMva_EE_DYJetsToLL->Clone();
  phoidMva_EE_DYJetsToLL_line->SetFillColor(0);
  phoidMva_EE_Data->GetXaxis()->SetRangeUser(-0.5,0.5);
  float max = phoidMva_EE_Data->GetMaximum();
  phoidMva_EE_Data->SetMaximum(max*1.1);
  phoidMva_EE_Data->Draw("e");
  phoidMva_EE_DYJetsToLL->Draw("hist,same");
  phoidMva_EE_syst->SetFillStyle(3013);
  phoidMva_EE_syst->SetFillColor(2);
  phoidMva_EE_syst->Draw("e2,same");
  phoidMva_EE_up_DYJetsToLL->Draw("hist,same");
  phoidMva_EE_down_DYJetsToLL->Draw("hist,same");
  phoidMva_EE_DYJetsToLL_line->Draw("hist,same");
  phoidMva_EE_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_idmva->cd(3);
  gPad->SetGrid();
  ratio_EB = (TH1*)phoidMva_EB_Data->Clone();
  ratio_EB_up = (TH1*)phoidMva_EB_up_DYJetsToLL->Clone();
  ratio_EB_down = (TH1*)phoidMva_EB_down_DYJetsToLL->Clone();
  ratio_EB->Divide(phoidMva_EB_DYJetsToLL);
  ratio_EB_up->Divide(phoidMva_EB_DYJetsToLL);
  ratio_EB_down->Divide(phoidMva_EB_DYJetsToLL);
  ratio_EB_up->SetLineColor(2);
  ratio_EB_down->SetLineColor(2);

  ratio_EB_syst = (TH1*)ratio_EB_up->Clone();
  ratio_EB_syst->Add(ratio_EB_down);
  ratio_EB_syst->Scale(0.5);
  TBox *ratio_EB_syst_box[400];
  for (int i=1; i<ratio_EB_syst->GetNbinsX()+1; i++) {
    float up = max(ratio_EB_up->GetBinContent(i),ratio_EB_down->GetBinContent(i));
    float down = min(ratio_EB_up->GetBinContent(i),ratio_EB_down->GetBinContent(i));
    ratio_EB_syst->SetBinError(i,fabs(up-down)/2.);
    ratio_EB_syst_box[i-1] = new TBox(ratio_EB_up->GetXaxis()->GetBinLowEdge(i),max(0.5,down),ratio_EB_up->GetXaxis()->GetBinLowEdge(i+1),min(2.,up));
    ratio_EB_syst_box[i-1]->SetFillStyle(3013);
    ratio_EB_syst_box[i-1]->SetFillColor(2);
  }

  ratio_EB_syst->SetFillStyle(3013);
  ratio_EB_syst->SetFillColor(2);
  ratio_EB_syst->SetMarkerStyle(0);
  ratio_EB_syst->SetMaximum(2.);
  ratio_EB_syst->SetMinimum(0.);
  ratio_EB_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_EB_syst->GetYaxis()->SetTitleSize(0.05);
  ratio_EB_syst->GetXaxis()->SetRangeUser(-0.5,0.5);
  ratio_EB_syst->Draw("e2");
  for (int i=1; i<ratio_EB_syst->GetNbinsX()+1; i++) {
    if (ratio_EB_syst_box[i-1]->GetY2()>0.5&&ratio_EB_syst_box[i-1]->GetY1()<2.) ratio_EB_syst_box[i-1]->Draw("hist,same");
    //ratio_EB_syst_box[i-1]->Draw("hist,same");
  }
  ratio_EB_up->Draw("hist,same");
  ratio_EB_down->Draw("hist,same");
  ratio_EB->Draw("e,same");
  TLine *line = new TLine(-0.5,1.,0.5,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  gPad->RedrawAxis();
  line->Draw();

  c_idmva->cd(4);
  gPad->SetGrid();
  ratio_EE = (TH1*)phoidMva_EE_Data->Clone();
  ratio_EE_up = (TH1*)phoidMva_EE_up_DYJetsToLL->Clone();
  ratio_EE_down = (TH1*)phoidMva_EE_down_DYJetsToLL->Clone();
  ratio_EE->Divide(phoidMva_EE_DYJetsToLL);
  ratio_EE_up->Divide(phoidMva_EE_DYJetsToLL);
  ratio_EE_down->Divide(phoidMva_EE_DYJetsToLL);
  ratio_EE_up->SetLineColor(2);
  ratio_EE_down->SetLineColor(2);

  ratio_EE_syst = (TH1*)ratio_EE_up->Clone();
  ratio_EE_syst->Add(ratio_EE_down);
  ratio_EE_syst->Scale(0.5);
  TBox *ratio_EE_syst_box[400];
  for (int i=1; i<ratio_EE_syst->GetNbinsX()+1; i++) {
    float up = max(ratio_EE_up->GetBinContent(i),ratio_EE_down->GetBinContent(i));
    float down = min(ratio_EE_up->GetBinContent(i),ratio_EE_down->GetBinContent(i));
    ratio_EE_syst->SetBinError(i,fabs(up-down)/2.);
    ratio_EE_syst_box[i-1] = new TBox(ratio_EE_up->GetXaxis()->GetBinLowEdge(i),max(0.5,down),ratio_EE_up->GetXaxis()->GetBinLowEdge(i+1),min(2.,up));
    ratio_EE_syst_box[i-1]->SetFillStyle(3013);
    ratio_EE_syst_box[i-1]->SetFillColor(2);
  }

  ratio_EE_syst->SetFillStyle(3013);
  ratio_EE_syst->SetFillColor(2);
  ratio_EE_syst->SetMarkerStyle(0);
  ratio_EE_syst->SetMaximum(2.);
  ratio_EE_syst->SetMinimum(0.);
  ratio_EE_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_EE_syst->GetYaxis()->SetTitleSize(0.05);
  ratio_EE_syst->GetXaxis()->SetRangeUser(-0.5,0.5);
  ratio_EE_syst->Draw("e2");
  for (int i=1; i<ratio_EE_syst->GetNbinsX()+1; i++) {
    if (ratio_EE_syst_box[i-1]->GetY2()>0.5&&ratio_EE_syst_box[i-1]->GetY1()<2.) ratio_EE_syst_box[i-1]->Draw("hist,same");
    //ratio_EE_syst_box[i-1]->Draw("hist,same");
  }
  ratio_EE_up->Draw("hist,same");
  ratio_EE_down->Draw("hist,same");
  ratio_EE->Draw("e,same");
  gPad->RedrawAxis();
  line->Draw();

  c_idmva->SaveAs("idmva.png");

}
