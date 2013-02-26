void sigmaE() {

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

  sigmaEOverE_EB_Data = (TH1*)pho1_sigmaEOverE_EB_cat0_Data->Clone();
  sigmaEOverE_EB_Data->Add(pho2_sigmaEOverE_EB_cat0_Data);
  sigmaEOverE_EE_Data = (TH1*)pho1_sigmaEOverE_EE_cat0_Data->Clone();
  sigmaEOverE_EE_Data->Add(pho2_sigmaEOverE_EE_cat0_Data);

  sigmaEOverE_EB_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EB_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EB_DYJetsToLL->Add(pho2_sigmaEOverE_EB_cat0_DYJetsToLL);
  sigmaEOverE_EE_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EE_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EE_DYJetsToLL->Add(pho2_sigmaEOverE_EE_cat0_DYJetsToLL);

  sigmaEOverE_EB_up_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EB_up_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EB_up_DYJetsToLL->Add(pho2_sigmaEOverE_EB_up_cat0_DYJetsToLL);
  sigmaEOverE_EE_up_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EE_up_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EE_up_DYJetsToLL->Add(pho2_sigmaEOverE_EE_up_cat0_DYJetsToLL);

  sigmaEOverE_EB_down_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EB_down_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EB_down_DYJetsToLL->Add(pho2_sigmaEOverE_EB_down_cat0_DYJetsToLL);
  sigmaEOverE_EE_down_DYJetsToLL = (TH1*)pho1_sigmaEOverE_EE_down_cat0_DYJetsToLL->Clone();
  sigmaEOverE_EE_down_DYJetsToLL->Add(pho2_sigmaEOverE_EE_down_cat0_DYJetsToLL);

  sigmaEOverE_EB_Data->Rebin(2);
  sigmaEOverE_EB_DYJetsToLL->Rebin(2);
  sigmaEOverE_EB_up_DYJetsToLL->Rebin(2);
  sigmaEOverE_EB_down_DYJetsToLL->Rebin(2);
  sigmaEOverE_EE_Data->Rebin(2);
  sigmaEOverE_EE_DYJetsToLL->Rebin(2);
  sigmaEOverE_EE_up_DYJetsToLL->Rebin(2);
  sigmaEOverE_EE_down_DYJetsToLL->Rebin(2);

  sigmaEOverE_EB_up_DYJetsToLL->SetLineColor(2);
  sigmaEOverE_EB_down_DYJetsToLL->SetLineColor(2);
  sigmaEOverE_EE_up_DYJetsToLL->SetLineColor(2);
  sigmaEOverE_EE_down_DYJetsToLL->SetLineColor(2);

  sigmaEOverE_EB_Data->GetYaxis()->SetTitle("");
  sigmaEOverE_EE_Data->GetYaxis()->SetTitle("");
  sigmaEOverE_EB_Data->GetXaxis()->SetTitle("#sigma_{E}/E (EB)");
  sigmaEOverE_EE_Data->GetXaxis()->SetTitle("#sigma_{E}/E (EE)");
  sigmaEOverE_EB_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaEOverE_EE_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaEOverE_EB_DYJetsToLL->SetFillColor(38);
  sigmaEOverE_EE_DYJetsToLL->SetFillColor(38);
  sigmaEOverE_EB_Data->SetMarkerStyle(20);
  sigmaEOverE_EE_Data->SetMarkerStyle(20);
  sigmaEOverE_EB_Data->SetMarkerSize(.8);
  sigmaEOverE_EE_Data->SetMarkerSize(.8);

  sigmaEOverE_EB_syst = (TH1*)sigmaEOverE_EB_up_DYJetsToLL->Clone();
  sigmaEOverE_EB_syst->Add(sigmaEOverE_EB_down_DYJetsToLL);
  sigmaEOverE_EB_syst->Scale(0.5);
  for (int i=1; i<sigmaEOverE_EB_syst->GetNbinsX()+1; i++) {
    float up = sigmaEOverE_EB_up_DYJetsToLL->GetBinContent(i);
    float down = sigmaEOverE_EB_down_DYJetsToLL->GetBinContent(i);
    sigmaEOverE_EB_syst->SetBinError(i,fabs(up-down)/2.);
  }

  sigmaEOverE_EE_syst = (TH1*)sigmaEOverE_EE_up_DYJetsToLL->Clone();
  sigmaEOverE_EE_syst->Add(sigmaEOverE_EE_down_DYJetsToLL);
  sigmaEOverE_EE_syst->Scale(0.5);
  for (int i=1; i<sigmaEOverE_EE_syst->GetNbinsX()+1; i++) {
    float up = sigmaEOverE_EE_up_DYJetsToLL->GetBinContent(i);
    float down = sigmaEOverE_EE_down_DYJetsToLL->GetBinContent(i);
    sigmaEOverE_EE_syst->SetBinError(i,fabs(up-down)/2.);
  }

  TLegend *leg2;
  leg2 = new TLegend(.5,.55,.87,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(sigmaEOverE_EB_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(sigmaEOverE_EB_DYJetsToLL,"DYJetsToLL MC","F");
  leg2->AddEntry(sigmaEOverE_EB_syst,"MC, #sigma_{E} scale factor 10%","F");
  //leg2->AddEntry(sigmaEOverE_EB_down_DYJetsToLL,"MC, #sigma_{E} scale factor -0.1");
  //leg2->AddEntry(sigmaEOverE_EB_up_DYJetsToLL,"MC #sigma_{E} scale factor +0.1");

  TCanvas *c_sigmae = new TCanvas("c_sigmae","sigma_E/E",1200,900);
  c_sigmae->Divide(2,2);

  c_sigmae->cd(1);
  float sf = sigmaEOverE_EB_Data->Integral()/sigmaEOverE_EB_DYJetsToLL->Integral();
  sigmaEOverE_EB_DYJetsToLL->Scale(sf);
  sigmaEOverE_EB_down_DYJetsToLL->Scale(sf);
  sigmaEOverE_EB_up_DYJetsToLL->Scale(sf);
  sigmaEOverE_EB_syst->Scale(sf);
  sigmaEOverE_EB_DYJetsToLL_line = (TH1*)sigmaEOverE_EB_DYJetsToLL->Clone();
  sigmaEOverE_EB_DYJetsToLL_line->SetFillColor(0);
  sigmaEOverE_EB_Data->GetXaxis()->SetRangeUser(0.,0.06);
  sigmaEOverE_EB_up_DYJetsToLL->GetXaxis()->SetRangeUser(0.,0.06);
  float max = sigmaEOverE_EB_Data->GetMaximum();
  sigmaEOverE_EB_Data->SetMaximum(max*1.1);
  sigmaEOverE_EB_Data->Draw("e");
  sigmaEOverE_EB_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EB_syst->SetFillStyle(3013);
  sigmaEOverE_EB_syst->SetFillColor(2);
  sigmaEOverE_EB_syst->Draw("e2,same");
  sigmaEOverE_EB_up_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EB_down_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EB_DYJetsToLL_line->Draw("hist,same");
  sigmaEOverE_EB_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_sigmae->cd(2);
  sf = sigmaEOverE_EE_Data->Integral()/sigmaEOverE_EE_DYJetsToLL->Integral();
  sigmaEOverE_EE_DYJetsToLL->Scale(sf);
  sigmaEOverE_EE_down_DYJetsToLL->Scale(sf);
  sigmaEOverE_EE_up_DYJetsToLL->Scale(sf);
  sigmaEOverE_EE_syst->Scale(sf);
  sigmaEOverE_EE_DYJetsToLL_line = (TH1*)sigmaEOverE_EE_DYJetsToLL->Clone();
  sigmaEOverE_EE_DYJetsToLL_line->SetFillColor(0);
  sigmaEOverE_EE_Data->GetXaxis()->SetRangeUser(0.,0.06);
  sigmaEOverE_EE_up_DYJetsToLL->GetXaxis()->SetRangeUser(0.,0.06);
  float max = sigmaEOverE_EE_Data->GetMaximum();
  sigmaEOverE_EE_Data->SetMaximum(max*1.1);
  sigmaEOverE_EE_Data->Draw("e");
  sigmaEOverE_EE_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EE_syst->SetFillStyle(3013);
  sigmaEOverE_EE_syst->SetFillColor(2);
  sigmaEOverE_EE_syst->Draw("e2,same");
  sigmaEOverE_EE_up_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EE_down_DYJetsToLL->Draw("hist,same");
  sigmaEOverE_EE_DYJetsToLL_line->Draw("hist,same");
  sigmaEOverE_EE_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_sigmae->cd(3);
  gPad->SetGrid();
  ratio_EB = (TH1*)sigmaEOverE_EB_Data->Clone();
  ratio_EB_up = (TH1*)sigmaEOverE_EB_up_DYJetsToLL->Clone();
  ratio_EB_down = (TH1*)sigmaEOverE_EB_down_DYJetsToLL->Clone();
  ratio_EB->Divide(sigmaEOverE_EB_DYJetsToLL);
  ratio_EB_up->Divide(sigmaEOverE_EB_DYJetsToLL);
  ratio_EB_down->Divide(sigmaEOverE_EB_DYJetsToLL);
  ratio_EB_up->SetLineColor(2);
  ratio_EB_down->SetLineColor(2);

  ratio_EB_syst = (TH1*)ratio_EB_up->Clone();
  ratio_EB_syst->Add(ratio_EB_down);
  ratio_EB_syst->Scale(0.5);
  for (int i=1; i<ratio_EB_syst->GetNbinsX()+1; i++) {
    float up = ratio_EB_up->GetBinContent(i);
    float down = ratio_EB_down->GetBinContent(i);
    ratio_EB_syst->SetBinError(i,fabs(up-down)/2.);
  }

  ratio_EB_syst->SetFillStyle(3013);
  ratio_EB_syst->SetFillColor(2);
  ratio_EB_syst->SetMarkerStyle(0);
  ratio_EB_syst->SetMaximum(2.);
  ratio_EB_syst->SetMinimum(0.5);
  ratio_EB_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_EB_syst->GetYaxis()->SetTitleSize(0.05);
  ratio_EB_syst->Draw("e2");
  ratio_EB_up->Draw("hist,same");
  ratio_EB_down->Draw("hist,same");
  ratio_EB->Draw("e,same");
  TLine *line = new TLine(0.,1.,0.06,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();

  c_sigmae->cd(4);
  gPad->SetGrid();
  ratio_EE = (TH1*)sigmaEOverE_EE_Data->Clone();
  ratio_EE_up = (TH1*)sigmaEOverE_EE_up_DYJetsToLL->Clone();
  ratio_EE_down = (TH1*)sigmaEOverE_EE_down_DYJetsToLL->Clone();
  ratio_EE->Divide(sigmaEOverE_EE_DYJetsToLL);
  ratio_EE_up->Divide(sigmaEOverE_EE_DYJetsToLL);
  ratio_EE_down->Divide(sigmaEOverE_EE_DYJetsToLL);
  ratio_EE_up->SetLineColor(2);
  ratio_EE_down->SetLineColor(2);

  ratio_EE_syst = (TH1*)ratio_EE_up->Clone();
  ratio_EE_syst->Add(ratio_EE_down);
  ratio_EE_syst->Scale(0.5);
  for (int i=1; i<ratio_EE_syst->GetNbinsX()+1; i++) {
    float up = ratio_EE_up->GetBinContent(i);
    float down = ratio_EE_down->GetBinContent(i);
    ratio_EE_syst->SetBinError(i,fabs(up-down)/2.);
  }

  ratio_EE_syst->SetFillStyle(3013);
  ratio_EE_syst->SetFillColor(2);
  ratio_EE_syst->SetMarkerStyle(0);
  ratio_EE_syst->SetMaximum(2.);
  ratio_EE_syst->SetMinimum(0.5);
  ratio_EE_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_EE_syst->GetYaxis()->SetTitleSize(0.05);
  ratio_EE_syst->Draw("e2");
  ratio_EE_up->Draw("hist,same");
  ratio_EE_down->Draw("hist,same");
  ratio_EE->Draw("e,same");
  line->Draw();

  c_sigmae->SaveAs("sigmaE.png");

}
