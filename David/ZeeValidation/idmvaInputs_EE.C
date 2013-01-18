void idmvaInputs_EE() {

  bool isEE=true;

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  TString det = isEE ? "EE" : "EB";

  pfchargedisogood03_EE_cat0_Data->GetXaxis()->SetTitle("pfchargedisogood03 (GeV) ("+det+")");
  pfchargedisobad03_EE_cat0_Data->GetXaxis()->SetTitle("pfchargedisobad03 (GeV) ("+det+")");
  pfphotoniso03_EE_cat0_Data->GetXaxis()->SetTitle("pfphotoniso03 (GeV) ("+det+")");
  s4ratio_EE_cat0_Data->GetXaxis()->SetTitle("s4ratio ("+det+")");
  sieie_EE_cat0_Data->GetXaxis()->SetTitle("sieie ("+det+")");
  sieip_EE_cat0_Data->GetXaxis()->SetTitle("covieip ("+det+")");
  etawidth_EE_cat0_Data->GetXaxis()->SetTitle("etawidth ("+det+")");
  phiwidth_EE_cat0_Data->GetXaxis()->SetTitle("phiwidth ("+det+")");
  if (isEE) ESEffSigmaRR_EE_cat0_Data->GetXaxis()->SetTitle("ESEffSigmaRR ("+det+")");
  r9_EE_cat0_Data->GetXaxis()->SetTitle("r9 ("+det+")");
  rho_EE_cat0_Data->GetXaxis()->SetTitle("rho ("+det+")");

  pfchargedisogood03_EE_cat0_Data->GetYaxis()->SetTitle("");
  pfchargedisobad03_EE_cat0_Data->GetYaxis()->SetTitle("");
  pfphotoniso03_EE_cat0_Data->GetYaxis()->SetTitle("");
  s4ratio_EE_cat0_Data->GetYaxis()->SetTitle("");
  sieie_EE_cat0_Data->GetYaxis()->SetTitle("");
  sieip_EE_cat0_Data->GetYaxis()->SetTitle("");
  etawidth_EE_cat0_Data->GetYaxis()->SetTitle("");
  phiwidth_EE_cat0_Data->GetYaxis()->SetTitle("");
  if (isEE) ESEffSigmaRR_EE_cat0_Data->GetYaxis()->SetTitle("");
  r9_EE_cat0_Data->GetYaxis()->SetTitle("");
  rho_EE_cat0_Data->GetYaxis()->SetTitle("");

  pfchargedisogood03_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pfchargedisobad03_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pfphotoniso03_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  s4ratio_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sieie_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sieip_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  etawidth_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  phiwidth_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  if (isEE) ESEffSigmaRR_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  r9_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  rho_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);

  pfchargedisogood03_EE_cat0_DYJetsToLL->SetFillColor(38);
  pfchargedisobad03_EE_cat0_DYJetsToLL->SetFillColor(38);
  pfphotoniso03_EE_cat0_DYJetsToLL->SetFillColor(38);
  s4ratio_EE_cat0_DYJetsToLL->SetFillColor(38);
  sieie_EE_cat0_DYJetsToLL->SetFillColor(38);
  sieip_EE_cat0_DYJetsToLL->SetFillColor(38);
  etawidth_EE_cat0_DYJetsToLL->SetFillColor(38);
  phiwidth_EE_cat0_DYJetsToLL->SetFillColor(38);
  if (isEE) ESEffSigmaRR_EE_cat0_DYJetsToLL->SetFillColor(38);
  r9_EE_cat0_DYJetsToLL->SetFillColor(38);
  rho_EE_cat0_DYJetsToLL->SetFillColor(38);

  pfchargedisogood03_EE_cat0_DYJetsToLL->SetLineColor(1);
  pfchargedisobad03_EE_cat0_DYJetsToLL->SetLineColor(1);
  pfphotoniso03_EE_cat0_DYJetsToLL->SetLineColor(1);
  s4ratio_EE_cat0_DYJetsToLL->SetLineColor(1);
  sieie_EE_cat0_DYJetsToLL->SetLineColor(1);
  sieip_EE_cat0_DYJetsToLL->SetLineColor(1);
  etawidth_EE_cat0_DYJetsToLL->SetLineColor(1);
  phiwidth_EE_cat0_DYJetsToLL->SetLineColor(1);
  if (isEE) ESEffSigmaRR_EE_cat0_DYJetsToLL->SetLineColor(1);
  r9_EE_cat0_DYJetsToLL->SetLineColor(1);
  rho_EE_cat0_DYJetsToLL->SetLineColor(1);

  pfchargedisogood03_EE_cat0_Data->SetMarkerStyle(20);
  pfchargedisobad03_EE_cat0_Data->SetMarkerStyle(20);
  pfphotoniso03_EE_cat0_Data->SetMarkerStyle(20);
  s4ratio_EE_cat0_Data->SetMarkerStyle(20);
  sieie_EE_cat0_Data->SetMarkerStyle(20);
  sieip_EE_cat0_Data->SetMarkerStyle(20);
  etawidth_EE_cat0_Data->SetMarkerStyle(20);
  phiwidth_EE_cat0_Data->SetMarkerStyle(20);
  if (isEE) ESEffSigmaRR_EE_cat0_Data->SetMarkerStyle(20);
  r9_EE_cat0_Data->SetMarkerStyle(20);
  rho_EE_cat0_Data->SetMarkerStyle(20);

  pfchargedisogood03_EE_cat0_Data->SetMarkerSize(0.4);
  pfchargedisobad03_EE_cat0_Data->SetMarkerSize(0.4);
  pfphotoniso03_EE_cat0_Data->SetMarkerSize(0.4);
  s4ratio_EE_cat0_Data->SetMarkerSize(0.4);
  sieie_EE_cat0_Data->SetMarkerSize(0.4);
  sieip_EE_cat0_Data->SetMarkerSize(0.4);
  etawidth_EE_cat0_Data->SetMarkerSize(0.4);
  phiwidth_EE_cat0_Data->SetMarkerSize(0.4);
  if (isEE) ESEffSigmaRR_EE_cat0_Data->SetMarkerSize(0.4);
  r9_EE_cat0_Data->SetMarkerSize(0.4);
  rho_EE_cat0_Data->SetMarkerSize(0.4);

  pfchargedisogood03_EE_cat0_Data->SetLineColor(1);
  pfchargedisobad03_EE_cat0_Data->SetLineColor(1);
  pfphotoniso03_EE_cat0_Data->SetLineColor(1);
  s4ratio_EE_cat0_Data->SetLineColor(1);
  sieie_EE_cat0_Data->SetLineColor(1);
  sieip_EE_cat0_Data->SetLineColor(1);
  etawidth_EE_cat0_Data->SetLineColor(1);
  phiwidth_EE_cat0_Data->SetLineColor(1);
  if (isEE) ESEffSigmaRR_EE_cat0_Data->SetLineColor(1);
  r9_EE_cat0_Data->SetLineColor(1);
  rho_EE_cat0_Data->SetLineColor(1);


  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(pfchargedisogood03_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(pfchargedisogood03_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(pfchargedisogood03_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(pfchargedisogood03_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg3;
  leg3 = new TLegend(.45,.65,.87,.87);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(10);
  leg3->SetTextSize(.03);
  leg3->AddEntry(pfchargedisogood03_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg3->AddEntry(pfchargedisogood03_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg4;
  leg4 = new TLegend(.2,.65,.52,.87);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(10);
  leg4->SetTextSize(.025);
  leg4->AddEntry(pfchargedisogood03_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg4->AddEntry(pfchargedisogood03_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);


  TCanvas *c_idmvain_1 = new TCanvas("c_idmvain_1","ID MVA input variables",1500,900);
  c_idmvain_1->Divide(4,3);

  c_idmvain_1->cd(1);
  float sf = pfchargedisogood03_EE_cat0_Data->Integral()/pfchargedisogood03_EE_cat0_DYJetsToLL->Integral();
  pfchargedisogood03_EE_cat0_DYJetsToLL->Scale(sf);
  pfchargedisogood03_EE_cat0_Data->Draw("e");
  pfchargedisogood03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfchargedisogood03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(2);
  sf = pfchargedisobad03_EE_cat0_Data->Integral()/pfchargedisobad03_EE_cat0_DYJetsToLL->Integral();
  pfchargedisobad03_EE_cat0_DYJetsToLL->Scale(sf);
  pfchargedisobad03_EE_cat0_Data->Draw("e");
  pfchargedisobad03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfchargedisobad03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(3);
  sf = pfphotoniso03_EE_cat0_Data->Integral()/pfphotoniso03_EE_cat0_DYJetsToLL->Integral();
  pfphotoniso03_EE_cat0_DYJetsToLL->Scale(sf);
  pfphotoniso03_EE_cat0_Data->Draw("e");
  pfphotoniso03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfphotoniso03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(4);
  sf = rho_EE_cat0_Data->Integral()/rho_EE_cat0_DYJetsToLL->Integral();
  rho_EE_cat0_DYJetsToLL->Scale(sf);
  rho_EE_cat0_Data->Draw("e");
  rho_EE_cat0_DYJetsToLL->Draw("hist,same");
  rho_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(5);
  gPad->SetLogy();
  pfchargedisogood03_EE_cat0_Data->Draw("e");
  pfchargedisogood03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfchargedisogood03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(6);
  gPad->SetLogy();
  pfchargedisobad03_EE_cat0_Data->Draw("e");
  pfchargedisobad03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfchargedisobad03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(7);
  gPad->SetLogy();
  pfphotoniso03_EE_cat0_Data->Draw("e");
  pfphotoniso03_EE_cat0_DYJetsToLL->Draw("hist,same");
  pfphotoniso03_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(8);
  gPad->SetLogy();
  rho_EE_cat0_Data->Draw("e");
  rho_EE_cat0_DYJetsToLL->Draw("hist,same");
  rho_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_1->cd(9);
  gPad->SetGrid();
  ratio_pfchargedisogood03_EE = (TH1*)pfchargedisogood03_EE_cat0_Data->Clone();
  ratio_pfchargedisogood03_EE->Sumw2();
  ratio_pfchargedisogood03_EE->Divide(pfchargedisogood03_EE_cat0_DYJetsToLL);
  ratio_pfchargedisogood03_EE->SetMaximum(1.8);
  ratio_pfchargedisogood03_EE->SetMinimum(0.2);
  ratio_pfchargedisogood03_EE->Draw("e");
  line1 = (TLine*)line->Clone();
  line1->SetX1(0.);
  line1->SetX2(8.);
  line1->Draw();

  c_idmvain_1->cd(10);
  gPad->SetGrid();
  ratio_pfchargedisobad03_EE = (TH1*)pfchargedisobad03_EE_cat0_Data->Clone();
  ratio_pfchargedisobad03_EE->Sumw2();
  ratio_pfchargedisobad03_EE->Divide(pfchargedisobad03_EE_cat0_DYJetsToLL);
  ratio_pfchargedisobad03_EE->SetMaximum(1.8);
  ratio_pfchargedisobad03_EE->SetMinimum(0.2);
  ratio_pfchargedisobad03_EE->Draw("e");
  line2 = (TLine*)line->Clone();
  line2->SetX1(0.);
  line2->SetX2(16.);
  line2->Draw();

  c_idmvain_1->cd(11);
  gPad->SetGrid();
  ratio_pfphotoniso03_EE = (TH1*)pfphotoniso03_EE_cat0_Data->Clone();
  ratio_pfphotoniso03_EE->Sumw2();
  ratio_pfphotoniso03_EE->Divide(pfphotoniso03_EE_cat0_DYJetsToLL);
  ratio_pfphotoniso03_EE->SetMaximum(1.8);
  ratio_pfphotoniso03_EE->SetMinimum(0.2);
  ratio_pfphotoniso03_EE->Draw("e");
  line3 = (TLine*)line->Clone();
  if (!isEE) {
    line3->SetX1(0.);
    line3->SetX2(12.);
  } else {
    line3->SetX1(0.);
    line3->SetX2(20.);
  }
  line3->Draw();

  c_idmvain_1->cd(12);
  gPad->SetGrid();
  ratio_rho_EE = (TH1*)rho_EE_cat0_Data->Clone();
  ratio_rho_EE->Sumw2();
  ratio_rho_EE->Divide(rho_EE_cat0_DYJetsToLL);
  ratio_rho_EE->SetMaximum(1.8);
  ratio_rho_EE->SetMinimum(0.2);
  ratio_rho_EE->Draw("e");
  line6 = (TLine*)line->Clone();
  if (!isEE) {
    line6->SetX1(-0.001);
    line6->SetX2(0.001);
  } else {
    line6->SetX1(-0.001);
    line6->SetX2(0.001);
  }
  line6->Draw();

  c_idmvain_1->SaveAs("idmvain_1_EE.png");


  TCanvas *c_idmvain_2 = new TCanvas("c_idmvain_2","ID MVA input variables",1500,900);
  c_idmvain_2->Divide(4,3);
  c_idmvain_2->cd(1);

  sf = sieie_EE_cat0_Data->Integral()/sieie_EE_cat0_DYJetsToLL->Integral();
  sieie_EE_cat0_DYJetsToLL->Scale(sf);
  sieie_EE_cat0_Data->Draw("e");
  sieie_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieie_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(2);
  sf = sieip_EE_cat0_Data->Integral()/sieip_EE_cat0_DYJetsToLL->Integral();
  sieip_EE_cat0_DYJetsToLL->Scale(sf);
  sieip_EE_cat0_Data->Draw("e");
  sieip_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieip_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(3);
  sf = etawidth_EE_cat0_Data->Integral()/etawidth_EE_cat0_DYJetsToLL->Integral();
  etawidth_EE_cat0_DYJetsToLL->Scale(sf);
  etawidth_EE_cat0_Data->Draw("e");
  etawidth_EE_cat0_DYJetsToLL->Draw("hist,same");
  etawidth_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(4);
  sf = phiwidth_EE_cat0_Data->Integral()/phiwidth_EE_cat0_DYJetsToLL->Integral();
  phiwidth_EE_cat0_DYJetsToLL->Scale(sf);
  phiwidth_EE_cat0_Data->Draw("e");
  phiwidth_EE_cat0_DYJetsToLL->Draw("hist,same");
  phiwidth_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(5);
  gPad->SetLogy();
  sieie_EE_cat0_Data->Draw("e");
  sieie_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieie_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(6);
  gPad->SetLogy();
  sieip_EE_cat0_Data->Draw("e");
  sieip_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieip_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(7);
  gPad->SetLogy();
  etawidth_EE_cat0_Data->Draw("e");
  etawidth_EE_cat0_DYJetsToLL->Draw("hist,same");
  etawidth_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(8);
  gPad->SetLogy();
  phiwidth_EE_cat0_Data->Draw("e");
  phiwidth_EE_cat0_DYJetsToLL->Draw("hist,same");
  phiwidth_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_idmvain_2->cd(9);
  gPad->SetGrid();
  ratio_sieie_EE = (TH1*)sieie_EE_cat0_Data->Clone();
  ratio_sieie_EE->Sumw2();
  ratio_sieie_EE->Divide(sieie_EE_cat0_DYJetsToLL);
  ratio_sieie_EE->SetMaximum(1.8);
  ratio_sieie_EE->SetMinimum(0.2);
  ratio_sieie_EE->Draw("e");
  line5 = (TLine*)line->Clone();
  if (!isEE) {
    line5->SetX1(0.004);
    line5->SetX2(0.014);
  } else {
    line5->SetX1(0.01);
    line5->SetX2(0.04);
  }
  line5->Draw();

  c_idmvain_2->cd(10);
  gPad->SetGrid();
  ratio_sieip_EE = (TH1*)sieip_EE_cat0_Data->Clone();
  ratio_sieip_EE->Sumw2();
  ratio_sieip_EE->Divide(sieip_EE_cat0_DYJetsToLL);
  ratio_sieip_EE->SetMaximum(1.8);
  ratio_sieip_EE->SetMinimum(0.2);
  ratio_sieip_EE->Draw("e");
  line6 = (TLine*)line->Clone();
  if (!isEE) {
    line6->SetX1(-0.00015);
    line6->SetX2(0.00015);
  } else {
    line6->SetX1(-0.00015);
    line6->SetX2(0.00015);
  }
  line6->Draw();

  c_idmvain_2->cd(11);
  gPad->SetGrid();
  ratio_etawidth_EE = (TH1*)etawidth_EE_cat0_Data->Clone();
  ratio_etawidth_EE->Sumw2();
  ratio_etawidth_EE->Divide(etawidth_EE_cat0_DYJetsToLL);
  ratio_etawidth_EE->SetMaximum(1.8);
  ratio_etawidth_EE->SetMinimum(0.2);
  ratio_etawidth_EE->Draw("e");
  line7 = (TLine*)line->Clone();
  if (!isEE) {
    line7->SetX1(0.);
    line7->SetX2(0.02);
  } else {
    line7->SetX1(0.);
    line7->SetX2(0.05);
  }
  line7->Draw();

  c_idmvain_2->cd(12);
  gPad->SetGrid();
  ratio_phiwidth_EE = (TH1*)phiwidth_EE_cat0_Data->Clone();
  ratio_phiwidth_EE->Sumw2();
  ratio_phiwidth_EE->Divide(phiwidth_EE_cat0_DYJetsToLL);
  ratio_phiwidth_EE->SetMaximum(1.8);
  ratio_phiwidth_EE->SetMinimum(0.2);
  ratio_phiwidth_EE->Draw("e");
  line8 = (TLine*)line->Clone();
  if (!isEE) {
    line8->SetX1(0.);
    line8->SetX2(0.12);
  } else {
    line8->SetX1(0.);
    line8->SetX2(0.14);
  }
  line8->Draw();

  c_idmvain_2->SaveAs("idmvain_2_EE.png");


  TCanvas *c_idmvain_3 = new TCanvas("c_idmvain_3","ID MVA input variables",1200,900);
  c_idmvain_3->Divide(3,3);

  c_idmvain_3->cd(1);
  sf = s4ratio_EE_cat0_Data->Integral()/s4ratio_EE_cat0_DYJetsToLL->Integral();
  s4ratio_EE_cat0_DYJetsToLL->Scale(sf);
  s4ratio_EE_cat0_Data->Draw("e");
  s4ratio_EE_cat0_DYJetsToLL->Draw("hist,same");
  s4ratio_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_idmvain_3->cd(2);
  sf = r9_EE_cat0_Data->Integral()/r9_EE_cat0_DYJetsToLL->Integral();
  r9_EE_cat0_DYJetsToLL->Scale(sf);
  r9_EE_cat0_Data->Draw("e");
  r9_EE_cat0_DYJetsToLL->Draw("hist,same");
  r9_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  if (isEE) {
    c_idmvain_3->cd(3);
    sf = ESEffSigmaRR_EE_cat0_Data->Integral()/ESEffSigmaRR_EE_cat0_DYJetsToLL->Integral();
    ESEffSigmaRR_EE_cat0_DYJetsToLL->Scale(sf);
    ESEffSigmaRR_EE_cat0_Data->Draw("e");
    ESEffSigmaRR_EE_cat0_DYJetsToLL->Draw("hist,same");
    ESEffSigmaRR_EE_cat0_Data->Draw("e,same");
    gPad->RedrawAxis();
    leg2->Draw();
  }

  c_idmvain_3->cd(4);
  gPad->SetLogy();
  s4ratio_EE_cat0_Data->Draw("e");
  s4ratio_EE_cat0_DYJetsToLL->Draw("hist,same");
  s4ratio_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_idmvain_3->cd(5);
  gPad->SetLogy();
  r9_EE_cat0_Data->Draw("e");
  r9_EE_cat0_DYJetsToLL->Draw("hist,same");
  r9_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  if (isEE) {
    c_idmvain_3->cd(6);
    sf = ESEffSigmaRR_EE_cat0_Data->Integral()/ESEffSigmaRR_EE_cat0_DYJetsToLL->Integral();
    ESEffSigmaRR_EE_cat0_DYJetsToLL->Scale(sf);
    ESEffSigmaRR_EE_cat0_Data->Draw("e");
    ESEffSigmaRR_EE_cat0_DYJetsToLL->Draw("hist,same");
    ESEffSigmaRR_EE_cat0_Data->Draw("e,same");
    gPad->RedrawAxis();
    leg2->Draw();
  }

  c_idmvain_3->cd(7);
  gPad->SetGrid();
  ratio_s4ratio_EE = (TH1*)s4ratio_EE_cat0_Data->Clone();
  ratio_s4ratio_EE->Sumw2();
  ratio_s4ratio_EE->Divide(s4ratio_EE_cat0_DYJetsToLL);
  ratio_s4ratio_EE->SetMaximum(1.8);
  ratio_s4ratio_EE->SetMinimum(0.2);
  ratio_s4ratio_EE->Draw("e");
  line8 = (TLine*)line->Clone();
  if (!isEE) {
    line8->SetX1(0.);
    line8->SetX2(0.12);
  } else {
    line8->SetX1(0.);
    line8->SetX2(0.14);
  }
  line8->Draw();

  c_idmvain_3->cd(8);
  gPad->SetGrid();
  ratio_r9_EE = (TH1*)r9_EE_cat0_Data->Clone();
  ratio_r9_EE->Sumw2();
  ratio_r9_EE->Divide(r9_EE_cat0_DYJetsToLL);
  ratio_r9_EE->SetMaximum(1.8);
  ratio_r9_EE->SetMinimum(0.2);
  ratio_r9_EE->Draw("e");
  line5 = (TLine*)line->Clone();
  if (!isEE) {
    line5->SetX1(0.004);
    line5->SetX2(0.014);
  } else {
    line5->SetX1(0.01);
    line5->SetX2(0.04);
  }
  line5->Draw();

  if (isEE) {
    c_idmvain_3->cd(9);
    gPad->SetGrid();
    ratio_ESEffSigmaRR_EE = (TH1*)ESEffSigmaRR_EE_cat0_Data->Clone();
    ratio_ESEffSigmaRR_EE->Sumw2();
    ratio_ESEffSigmaRR_EE->Divide(ESEffSigmaRR_EE_cat0_DYJetsToLL);
    ratio_ESEffSigmaRR_EE->SetMaximum(1.8);
    ratio_ESEffSigmaRR_EE->SetMinimum(0.2);
    ratio_ESEffSigmaRR_EE->Draw("e");
    line8 = (TLine*)line->Clone();
    if (!isEE) {
      line8->SetX1(0.);
      line8->SetX2(0.12);
    } else {
      line8->SetX1(0.);
      line8->SetX2(0.14);
    }
    line8->Draw();
  }

  c_idmvain_3->SaveAs("idmvain_3_EE.png");

}
