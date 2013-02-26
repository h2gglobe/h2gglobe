void CiCInputs_EE() {

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

  cic_isosumoet_EE_cat0_Data->GetXaxis()->SetTitle("cic_isosumoet (GeV) ("+det+")");
  cic_isosumoetbad_EE_cat0_Data->GetXaxis()->SetTitle("cic_isosumoetbad (GeV) ("+det+")");
  cic_trkisooet_EE_cat0_Data->GetXaxis()->SetTitle("cic_trkisooet (GeV) ("+det+")");
  cic_hoe_EE_cat0_Data->GetXaxis()->SetTitle("cic_hoe ("+det+")");
  sieie_EE_cat0_Data->GetXaxis()->SetTitle("sieie ("+det+")");
  r9_EE_cat0_Data->GetXaxis()->SetTitle("r9 ("+det+")");

  cic_isosumoet_EE_cat0_Data->GetYaxis()->SetTitle("");
  cic_isosumoetbad_EE_cat0_Data->GetYaxis()->SetTitle("");
  cic_trkisooet_EE_cat0_Data->GetYaxis()->SetTitle("");
  cic_hoe_EE_cat0_Data->GetYaxis()->SetTitle("");
  sieie_EE_cat0_Data->GetYaxis()->SetTitle("");
  r9_EE_cat0_Data->GetYaxis()->SetTitle("");

  cic_isosumoet_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  cic_isosumoetbad_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  cic_trkisooet_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  cic_hoe_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sieie_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  r9_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);

  cic_isosumoet_EE_cat0_DYJetsToLL->SetFillColor(38);
  cic_isosumoetbad_EE_cat0_DYJetsToLL->SetFillColor(38);
  cic_trkisooet_EE_cat0_DYJetsToLL->SetFillColor(38);
  cic_hoe_EE_cat0_DYJetsToLL->SetFillColor(38);
  sieie_EE_cat0_DYJetsToLL->SetFillColor(38);
  r9_EE_cat0_DYJetsToLL->SetFillColor(38);

  cic_isosumoet_EE_cat0_DYJetsToLL->SetLineColor(1);
  cic_isosumoetbad_EE_cat0_DYJetsToLL->SetLineColor(1);
  cic_trkisooet_EE_cat0_DYJetsToLL->SetLineColor(1);
  cic_hoe_EE_cat0_DYJetsToLL->SetLineColor(1);
  sieie_EE_cat0_DYJetsToLL->SetLineColor(1);
  r9_EE_cat0_DYJetsToLL->SetLineColor(1);

  cic_isosumoet_EE_cat0_Data->SetMarkerStyle(20);
  cic_isosumoetbad_EE_cat0_Data->SetMarkerStyle(20);
  cic_trkisooet_EE_cat0_Data->SetMarkerStyle(20);
  cic_hoe_EE_cat0_Data->SetMarkerStyle(20);
  sieie_EE_cat0_Data->SetMarkerStyle(20);
  r9_EE_cat0_Data->SetMarkerStyle(20);

  cic_isosumoet_EE_cat0_Data->SetMarkerSize(0.4);
  cic_isosumoetbad_EE_cat0_Data->SetMarkerSize(0.4);
  cic_trkisooet_EE_cat0_Data->SetMarkerSize(0.4);
  cic_hoe_EE_cat0_Data->SetMarkerSize(0.4);
  sieie_EE_cat0_Data->SetMarkerSize(0.4);
  r9_EE_cat0_Data->SetMarkerSize(0.4);

  cic_isosumoet_EE_cat0_Data->SetLineColor(1);
  cic_isosumoetbad_EE_cat0_Data->SetLineColor(1);
  cic_trkisooet_EE_cat0_Data->SetLineColor(1);
  cic_hoe_EE_cat0_Data->SetLineColor(1);
  sieie_EE_cat0_Data->SetLineColor(1);
  r9_EE_cat0_Data->SetLineColor(1);


  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(cic_isosumoet_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(cic_isosumoet_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(cic_isosumoet_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(cic_isosumoet_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg3;
  leg3 = new TLegend(.45,.65,.87,.87);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(10);
  leg3->SetTextSize(.03);
  leg3->AddEntry(cic_isosumoet_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg3->AddEntry(cic_isosumoet_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg4;
  leg4 = new TLegend(.2,.65,.52,.87);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(10);
  leg4->SetTextSize(.025);
  leg4->AddEntry(cic_isosumoet_EE_cat0_Data,"Data (19.6fb^{-1})");
  leg4->AddEntry(cic_isosumoet_EE_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);


  TCanvas *c_cicin_1 = new TCanvas("c_cicin_1","ID MVA input variables",1500,900);
  c_cicin_1->Divide(3,3);

  c_cicin_1->cd(1);
  float sf = cic_isosumoet_EE_cat0_Data->Integral()/cic_isosumoet_EE_cat0_DYJetsToLL->Integral();
  cic_isosumoet_EE_cat0_DYJetsToLL->Scale(sf);
  cic_isosumoet_EE_cat0_Data->Draw("e");
  cic_isosumoet_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_isosumoet_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(2);
  sf = cic_isosumoetbad_EE_cat0_Data->Integral()/cic_isosumoetbad_EE_cat0_DYJetsToLL->Integral();
  cic_isosumoetbad_EE_cat0_DYJetsToLL->Scale(sf);
  cic_isosumoetbad_EE_cat0_Data->Draw("e");
  cic_isosumoetbad_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_isosumoetbad_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(3);
  sf = cic_trkisooet_EE_cat0_Data->Integral()/cic_trkisooet_EE_cat0_DYJetsToLL->Integral();
  cic_trkisooet_EE_cat0_DYJetsToLL->Scale(sf);
  cic_trkisooet_EE_cat0_Data->Draw("e");
  cic_trkisooet_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_trkisooet_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(4);
  gPad->SetLogy();
  cic_isosumoet_EE_cat0_Data->Draw("e");
  cic_isosumoet_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_isosumoet_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(5);
  gPad->SetLogy();
  cic_isosumoetbad_EE_cat0_Data->Draw("e");
  cic_isosumoetbad_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_isosumoetbad_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(6);
  gPad->SetLogy();
  cic_trkisooet_EE_cat0_Data->Draw("e");
  cic_trkisooet_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_trkisooet_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_1->cd(7);
  gPad->SetGrid();
  ratio_cic_isosumoet_EE = (TH1*)cic_isosumoet_EE_cat0_Data->Clone();
  ratio_cic_isosumoet_EE->Divide(cic_isosumoet_EE_cat0_DYJetsToLL);
  ratio_cic_isosumoet_EE->SetMaximum(1.8);
  ratio_cic_isosumoet_EE->SetMinimum(0.2);
  ratio_cic_isosumoet_EE->Draw("e");
  line1 = (TLine*)line->Clone();
  line1->SetX1(-4.);
  line1->SetX2(20.);
  line1->Draw();

  c_cicin_1->cd(8);
  gPad->SetGrid();
  ratio_cic_isosumoetbad_EE = (TH1*)cic_isosumoetbad_EE_cat0_Data->Clone();
  ratio_cic_isosumoetbad_EE->Divide(cic_isosumoetbad_EE_cat0_DYJetsToLL);
  ratio_cic_isosumoetbad_EE->SetMaximum(1.8);
  ratio_cic_isosumoetbad_EE->SetMinimum(0.2);
  ratio_cic_isosumoetbad_EE->Draw("e");
  line2 = (TLine*)line->Clone();
  line2->SetX1(-4.);
  line2->SetX2(20.);
  line2->Draw();

  c_cicin_1->cd(9);
  gPad->SetGrid();
  ratio_cic_trkisooet_EE = (TH1*)cic_trkisooet_EE_cat0_Data->Clone();
  ratio_cic_trkisooet_EE->Divide(cic_trkisooet_EE_cat0_DYJetsToLL);
  ratio_cic_trkisooet_EE->SetMaximum(1.8);
  ratio_cic_trkisooet_EE->SetMinimum(0.2);
  ratio_cic_trkisooet_EE->Draw("e");
  line3 = (TLine*)line->Clone();
  line3->SetX1(0.);
  line3->SetX2(12.);
  line3->Draw();

  c_cicin_1->SaveAs("CiCinput1_EE.png");


  TCanvas *c_cicin_2 = new TCanvas("c_cicin_2","ID MVA input variables",1500,900);
  c_cicin_2->Divide(3,3);
  c_cicin_2->cd(3);

  sf = sieie_EE_cat0_Data->Integral()/sieie_EE_cat0_DYJetsToLL->Integral();
  sieie_EE_cat0_DYJetsToLL->Scale(sf);
  sieie_EE_cat0_Data->Draw("e");
  sieie_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieie_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_2->cd(6);
  gPad->SetLogy();
  sieie_EE_cat0_Data->Draw("e");
  sieie_EE_cat0_DYJetsToLL->Draw("hist,same");
  sieie_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_cicin_2->cd(9);
  gPad->SetGrid();
  ratio_sieie_EE = (TH1*)sieie_EE_cat0_Data->Clone();
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

  c_cicin_2->cd(1);
  sf = cic_hoe_EE_cat0_Data->Integral()/cic_hoe_EE_cat0_DYJetsToLL->Integral();
  cic_hoe_EE_cat0_DYJetsToLL->Scale(sf);
  cic_hoe_EE_cat0_Data->Draw("e");
  cic_hoe_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_hoe_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_cicin_2->cd(2);
  sf = r9_EE_cat0_Data->Integral()/r9_EE_cat0_DYJetsToLL->Integral();
  r9_EE_cat0_DYJetsToLL->Scale(sf);
  r9_EE_cat0_Data->Draw("e");
  r9_EE_cat0_DYJetsToLL->Draw("hist,same");
  r9_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_cicin_2->cd(4);
  gPad->SetLogy();
  cic_hoe_EE_cat0_Data->Draw("e");
  cic_hoe_EE_cat0_DYJetsToLL->Draw("hist,same");
  cic_hoe_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_cicin_2->cd(5);
  gPad->SetLogy();
  r9_EE_cat0_Data->Draw("e");
  r9_EE_cat0_DYJetsToLL->Draw("hist,same");
  r9_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_cicin_2->cd(7);
  gPad->SetGrid();
  ratio_cic_hoe_EE = (TH1*)cic_hoe_EE_cat0_Data->Clone();
  ratio_cic_hoe_EE->Divide(cic_hoe_EE_cat0_DYJetsToLL);
  ratio_cic_hoe_EE->SetMaximum(1.8);
  ratio_cic_hoe_EE->SetMinimum(0.2);
  ratio_cic_hoe_EE->Draw("e");
  line8 = (TLine*)line->Clone();
  if (!isEE) {
    line8->SetX1(0.);
    line8->SetX2(0.1);
  } else {
    line8->SetX1(0.);
    line8->SetX2(0.1);
  }
  line8->Draw();

  c_cicin_2->cd(8);
  gPad->SetGrid();
  ratio_r9_EE = (TH1*)r9_EE_cat0_Data->Clone();
  ratio_r9_EE->Divide(r9_EE_cat0_DYJetsToLL);
  ratio_r9_EE->SetMaximum(1.8);
  ratio_r9_EE->SetMinimum(0.2);
  ratio_r9_EE->Draw("e");
  line5 = (TLine*)line->Clone();
  line5->SetX1(0.);
  line5->SetX2(1.);
  line5->Draw();

  c_cicin_2->SaveAs("CiCinput2_EE.png");

}
