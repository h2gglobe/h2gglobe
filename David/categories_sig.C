void categories_sig() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TText *text = new TText();
  text->SetNDC();
  text->SetTextSize(0.05);

  TFile *file = TFile::Open("histograms_CMS-HGG_categories_sig.root");
  file->cd();

  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output");

  bdtout_all_tot = (TH1*)bdtout_cat0_tot->Clone();
  bdtout_all_tot->Add(bdtout_cat1_tot);
  bdtout_all_tot->Add(bdtout_cat2_tot);
  bdtout_all_tot->Add(bdtout_cat3_tot);

  cout << bdtout_all_tot->Integral() << endl;
  cout << bdtout_all_tot->Integral(49,100)+bdtout_all_tot->GetBinContent(48)/2. << endl;
  cout << bdtout_all_tot->GetBinLowEdge(49) << endl;
  cout << bdtout_all_tot->GetBinLowEdge(48) << endl;
  cout << (bdtout_all_tot->Integral(49,100)+bdtout_all_tot->GetBinContent(48)/2.)/bdtout_all_tot->Integral() << endl;

  bdtout_passCiC_all_tot = (TH1*)bdtout_passCiC_cat0_tot->Clone();
  bdtout_passCiC_all_tot->Add(bdtout_passCiC_cat1_tot);
  bdtout_passCiC_all_tot->Add(bdtout_passCiC_cat2_tot);
  bdtout_passCiC_all_tot->Add(bdtout_passCiC_cat3_tot);

  bdtout_failpresel_cat0_tot->Rebin(2);
  bdtout_failpresel_cat1_tot->Rebin(2);
  bdtout_failpresel_cat2_tot->Rebin(2);
  bdtout_failpresel_cat3_tot->Rebin(2);

  bdtout_failpresel_all_tot = (TH1*)bdtout_failpresel_cat0_tot->Clone();
  bdtout_failpresel_all_tot->Add(bdtout_failpresel_cat1_tot);
  bdtout_failpresel_all_tot->Add(bdtout_failpresel_cat2_tot);
  bdtout_failpresel_all_tot->Add(bdtout_failpresel_cat3_tot);

  float nFail_m100180_bdtout005 = bdtout_all_tot->Integral(1,52) + bdtout_all_tot->GetBinContent(53)/2.;
  float nFail_m100180 = bdtout_all_tot->Integral(1,100);
  float frac_bdtout005 = nFail_m100180_bdtout005/nFail_m100180;
  float frac_bdtout005_err = frac_err(nFail_m100180_bdtout005,nFail_m100180);
  cout << bdtout_all_tot->Integral(1,52) << " " <<  bdtout_all_tot->GetBinContent(53)/2. << " " << bdtout_all_tot->Integral(1,100) << " " << bdtout_all_tot->GetBinLowEdge(53) << " " << frac_bdtout005 << " ± " << frac_bdtout005_err << endl;

  bdtout_all_tot->SetLineColor(1);
  bdtout_cat0_tot->SetLineColor(2);
  bdtout_cat1_tot->SetLineColor(3);
  bdtout_cat2_tot->SetLineColor(4);
  bdtout_cat3_tot->SetLineColor(6);

  bdtout_all_tot->SetLineWidth(2);
  bdtout_cat0_tot->SetLineWidth(2);
  bdtout_cat1_tot->SetLineWidth(2);
  bdtout_cat2_tot->SetLineWidth(2);
  bdtout_cat3_tot->SetLineWidth(2);

  bdtout_cat0_tot->SetMarkerColor(2);
  bdtout_cat1_tot->SetMarkerColor(3);
  bdtout_cat2_tot->SetMarkerColor(4);
  bdtout_cat3_tot->SetMarkerColor(6);

  bdtout_all_tot->SetMarkerStyle(20);
  bdtout_cat0_tot->SetMarkerStyle(20);
  bdtout_cat1_tot->SetMarkerStyle(20);
  bdtout_cat2_tot->SetMarkerStyle(20);
  bdtout_cat3_tot->SetMarkerStyle(20);

  bdtout_all_tot->SetMarkerSize(0.7);
  bdtout_cat0_tot->SetMarkerSize(0.7);
  bdtout_cat1_tot->SetMarkerSize(0.7);
  bdtout_cat2_tot->SetMarkerSize(0.7);
  bdtout_cat3_tot->SetMarkerSize(0.7);

  bdtout_passCiC_all_tot->SetLineColor(1);
  bdtout_passCiC_cat0_tot->SetLineColor(2);
  bdtout_passCiC_cat1_tot->SetLineColor(3);
  bdtout_passCiC_cat2_tot->SetLineColor(4);
  bdtout_passCiC_cat3_tot->SetLineColor(6);

  bdtout_passCiC_lowPt_cat0_tot->SetLineColor(2);
  bdtout_passCiC_lowPt_cat1_tot->SetLineColor(3);
  bdtout_passCiC_lowPt_cat2_tot->SetLineColor(4);
  bdtout_passCiC_lowPt_cat3_tot->SetLineColor(6);

  bdtout_passCiC_highPt_cat0_tot->SetLineColor(2);
  bdtout_passCiC_highPt_cat1_tot->SetLineColor(3);
  bdtout_passCiC_highPt_cat2_tot->SetLineColor(4);
  bdtout_passCiC_highPt_cat3_tot->SetLineColor(6);

  bdtout_failpresel_all_tot->SetLineColor(1);
  bdtout_failpresel_cat0_tot->SetLineColor(2);
  bdtout_failpresel_cat1_tot->SetLineColor(3);
  bdtout_failpresel_cat2_tot->SetLineColor(4);
  bdtout_failpresel_cat3_tot->SetLineColor(6);


  bdtout_all_tot->GetXaxis()->SetTitle("di-photon MVA output");

  float boundaries[4] = {-0.05,0.49,0.79,.91};
  float max = bdtout_all_tot->GetMaximum();

  TBox* box = new TBox(-1.,0.,boundaries[0],max*1.05);
  box->SetFillColor(38);
  box->SetFillStyle(3002);

  bdtout_all_tot->Draw("hist");
  box->Draw("hist,same");
  bdtout_all_tot->Draw("hist,same");
  bdtout_cat0_tot->Draw("hist,same");
  bdtout_cat1_tot->Draw("hist,same");
  bdtout_cat2_tot->Draw("hist,same");
  bdtout_cat3_tot->Draw("hist,same");

  TLegend *leg;
  leg = new TLegend(.14,.6,.46,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(bdtout_all_tot,"All");
  leg->AddEntry(bdtout_cat0_tot,"both EB, both R9>0.94");
  leg->AddEntry(bdtout_cat1_tot,"both EB, !both R9>0.94");
  leg->AddEntry(bdtout_cat2_tot,"!both EB, both R9>0.94");
  leg->AddEntry(bdtout_cat3_tot,"!both EB, !both R9>0.94");
  leg->Draw("hist");

  TLine* line[4];
  for (int i=0; i<4; i++) {
    line[i] = new TLine(boundaries[i],0.,boundaries[i],max*1.05);
    line[i]->SetLineColor(4);
    line[i]->SetLineWidth(2);
    line[i]->SetLineStyle(9);
    line[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout->SaveAs("categories_sig.png");


  TCanvas *c_bdtout_passCiC = new TCanvas("c_bdtout_passCiC","BDT output, pass CiC");

  bdtout_passCiC_cat0_tot_clone = (TH1*)bdtout_passCiC_cat0_tot->Clone();
  bdtout_passCiC_cat1_tot_clone = (TH1*)bdtout_passCiC_cat1_tot->Clone();
  bdtout_passCiC_cat2_tot_clone = (TH1*)bdtout_passCiC_cat2_tot->Clone();
  bdtout_passCiC_cat3_tot_clone = (TH1*)bdtout_passCiC_cat3_tot->Clone();

  bdtout_passCiC_all_tot->SetLineWidth(2);
  bdtout_passCiC_cat0_tot_clone->SetLineWidth(2);
  bdtout_passCiC_cat1_tot_clone->SetLineWidth(2);
  bdtout_passCiC_cat2_tot_clone->SetLineWidth(2);
  bdtout_passCiC_cat3_tot_clone->SetLineWidth(2);

  bdtout_passCiC_all_tot->GetXaxis()->SetTitle("di-photon MVA output");

  max = bdtout_passCiC_all_tot->GetMaximum();
  TBox* box_passCiC = new TBox(-1.,0.,boundaries[0],max*1.05);
  box_passCiC->SetFillColor(38);
  box_passCiC->SetFillStyle(3002);

  bdtout_passCiC_all_tot->Draw("hist");
  box_passCiC->Draw("hist,same");
  bdtout_passCiC_all_tot->Draw("hist,same");
  bdtout_passCiC_cat0_tot_clone->Draw("hist,same");
  bdtout_passCiC_cat1_tot_clone->Draw("hist,same");
  bdtout_passCiC_cat2_tot_clone->Draw("hist,same");
  bdtout_passCiC_cat3_tot_clone->Draw("hist,same");

  leg->Draw("hist");

  TLine* line_passCiC[4];
  for (int i=0; i<4; i++) {
    line_passCiC[i] = new TLine(boundaries[i],0.,boundaries[i],max*1.05);
    line_passCiC[i]->SetLineColor(4);
    line_passCiC[i]->SetLineWidth(2);
    line_passCiC[i]->SetLineStyle(9);
    line_passCiC[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout_passCiC->SaveAs("categories_passCiC_sig.png");

  /*
  TCanvas *c_bdtout_failpresel = new TCanvas("c_bdtout_failpresel","BDT output, pass CiC, fail presel");

  bdtout_failpresel_cat0_tot_clone = (TH1*)bdtout_failpresel_cat0_tot->Clone();
  bdtout_failpresel_cat1_tot_clone = (TH1*)bdtout_failpresel_cat1_tot->Clone();
  bdtout_failpresel_cat2_tot_clone = (TH1*)bdtout_failpresel_cat2_tot->Clone();
  bdtout_failpresel_cat3_tot_clone = (TH1*)bdtout_failpresel_cat3_tot->Clone();

  bdtout_failpresel_all_tot->SetLineWidth(2);
  bdtout_failpresel_cat0_tot_clone->SetLineWidth(2);
  bdtout_failpresel_cat1_tot_clone->SetLineWidth(2);
  bdtout_failpresel_cat2_tot_clone->SetLineWidth(2);
  bdtout_failpresel_cat3_tot_clone->SetLineWidth(2);

  bdtout_failpresel_all_tot->GetXaxis()->SetTitle("di-photon MVA output");

  max = bdtout_failpresel_all_tot->GetMaximum();
  TBox* box_failpresel = new TBox(-1.,0.,boundaries[0],max*1.05);
  box_failpresel->SetFillColor(38);
  box_failpresel->SetFillStyle(3002);

  bdtout_failpresel_all_tot->Draw("hist");
  box_failpresel->Draw("hist,same");
  bdtout_failpresel_all_tot->Draw("hist,same");
  bdtout_failpresel_cat0_tot_clone->Draw("hist,same");
  bdtout_failpresel_cat1_tot_clone->Draw("hist,same");
  bdtout_failpresel_cat2_tot_clone->Draw("hist,same");
  bdtout_failpresel_cat3_tot_clone->Draw("hist,same");

  leg->Draw("hist");

  TLine* line_failpresel[4];
  for (int i=0; i<4; i++) {
    line_failpresel[i] = new TLine(boundaries[i],0.,boundaries[i],max*1.05);
    line_failpresel[i]->SetLineColor(4);
    line_failpresel[i]->SetLineWidth(2);
    line_failpresel[i]->SetLineStyle(9);
    line_failpresel[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout_failpresel->SaveAs("categories_failpresel_sig.png");
  */

  TCanvas *c_bdtout_compareCiC = new TCanvas("c_bdtout_compareCiC","BDT output: pass CiC supertight",1000,650);
  c_bdtout_compareCiC->Divide(2,2);

  bdtout_passCiC_cat0_tot->SetFillColor(2);
  bdtout_passCiC_cat1_tot->SetFillColor(3);
  bdtout_passCiC_cat2_tot->SetFillColor(4);
  bdtout_passCiC_cat3_tot->SetFillColor(6);

  bdtout_passCiC_cat0_tot->SetFillStyle(3002);
  bdtout_passCiC_cat1_tot->SetFillStyle(3002);
  bdtout_passCiC_cat2_tot->SetFillStyle(3002);
  bdtout_passCiC_cat3_tot->SetFillStyle(3002);

  bdtout_cat0_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat1_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat2_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat3_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat0_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat1_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat2_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat3_tot->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC->cd(1);

  box_passCiC_cat0 = (TBox*)box->Clone();
  box_passCiC_cat0->SetY2(bdtout_cat0_tot->GetMaximum()*1.05);

  bdtout_cat0_tot->Draw("hist");
  box_passCiC_cat0->Draw("hist,same");
  bdtout_cat0_tot->Draw("hist,same");
  bdtout_passCiC_cat0_tot->Draw("hist,same");

  TLine* line_passCiC_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat0[i]->SetY2(bdtout_cat0_tot->GetMaximum()*1.05);
    line_passCiC_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(2);

  box_passCiC_cat1 = (TBox*)box->Clone();
  box_passCiC_cat1->SetY2(bdtout_cat1_tot->GetMaximum()*1.05);

  bdtout_cat1_tot->Draw("hist");
  box_passCiC_cat1->Draw("hist,same");
  bdtout_cat1_tot->Draw("hist,same");
  bdtout_passCiC_cat1_tot->Draw("hist,same");

  TLine* line_passCiC_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat1[i]->SetY2(bdtout_cat1_tot->GetMaximum()*1.05);
    line_passCiC_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(3);

  box_passCiC_cat2 = (TBox*)box->Clone();
  box_passCiC_cat2->SetY2(bdtout_cat2_tot->GetMaximum()*1.05);

  bdtout_cat2_tot->Draw("hist");
  box_passCiC_cat2->Draw("hist,same");
  bdtout_cat2_tot->Draw("hist,same");
  bdtout_passCiC_cat2_tot->Draw("hist,same");

  TLine* line_passCiC_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat2[i]->SetY2(bdtout_cat2_tot->GetMaximum()*1.05);
    line_passCiC_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(4);

  box_passCiC_cat3 = (TBox*)box->Clone();
  box_passCiC_cat3->SetY2(bdtout_cat3_tot->GetMaximum()*1.05);

  bdtout_cat3_tot->Draw("hist");
  box_passCiC_cat3->Draw("hist,same");
  bdtout_cat3_tot->Draw("hist,same");
  bdtout_passCiC_cat3_tot->Draw("hist,same");

  TLine* line_passCiC_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat3[i]->SetY2(bdtout_cat3_tot->GetMaximum()*1.05);
    line_passCiC_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->SaveAs("categories_compareCiC_sig.png");



  TCanvas *c_bdtout_lowPt = new TCanvas("c_bdtout_lowPt","BDT output: di-photon pT<40 GeV");

  bdtout_lowPt_all_tot = (TH1*)bdtout_lowPt_cat0_tot->Clone();
  bdtout_lowPt_all_tot->Add(bdtout_lowPt_cat1_tot);
  bdtout_lowPt_all_tot->Add(bdtout_lowPt_cat2_tot);
  bdtout_lowPt_all_tot->Add(bdtout_lowPt_cat3_tot);

  bdtout_lowPt_all_tot->SetLineColor(1);
  bdtout_lowPt_cat0_tot->SetLineColor(2);
  bdtout_lowPt_cat1_tot->SetLineColor(3);
  bdtout_lowPt_cat2_tot->SetLineColor(4);
  bdtout_lowPt_cat3_tot->SetLineColor(6);

  bdtout_lowPt_all_tot->SetLineWidth(2);
  bdtout_lowPt_cat0_tot->SetLineWidth(2);
  bdtout_lowPt_cat1_tot->SetLineWidth(2);
  bdtout_lowPt_cat2_tot->SetLineWidth(2);
  bdtout_lowPt_cat3_tot->SetLineWidth(2);

  bdtout_lowPt_all_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_all_tot->GetYaxis()->SetTitle("Events/0.02");
  bdtout_lowPt_all_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_all_tot->GetXaxis()->SetLabelSize(0.05);
  bdtout_lowPt_all_tot->GetYaxis()->SetTitleSize(0.05);
  bdtout_lowPt_all_tot->GetYaxis()->SetLabelSize(0.05);

  max = bdtout_lowPt_all_tot->GetMaximum();
  TBox* box_lowPt = new TBox(-1.,0.,boundaries[0],max*1.05);
  box_lowPt->SetFillColor(38);
  box_lowPt->SetFillStyle(3002);

  bdtout_lowPt_all_tot->Draw("hist");
  box_lowPt->Draw("hist,same");
  bdtout_lowPt_all_tot->Draw("hist,same");
  bdtout_lowPt_cat0_tot->Draw("hist,same");
  bdtout_lowPt_cat1_tot->Draw("hist,same");
  bdtout_lowPt_cat2_tot->Draw("hist,same");
  bdtout_lowPt_cat3_tot->Draw("hist,same");

  TLegend *leg_lowPt = (TLegend*)leg->Clone();
  leg_lowPt->Clear();
  leg_lowPt->AddEntry(bdtout_all_tot,"All (p_{T}(#gamma#gamma) < 40 GeV)");
  leg_lowPt->AddEntry(bdtout_cat0_tot,"both EB, both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat1_tot,"both EB, !both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat2_tot,"!both EB, both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat3_tot,"!both EB, !both R9>0.94");
  leg_lowPt->Draw("hist");

  TLine* line_lowPt[4];
  for (int i=0; i<4; i++) {
    line_lowPt[i] = new TLine(boundaries[i],0.,boundaries[i],max*1.05);
    line_lowPt[i]->SetLineColor(4);
    line_lowPt[i]->SetLineWidth(2);
    line_lowPt[i]->SetLineStyle(9);
    line_lowPt[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout_lowPt->SaveAs("categories_lowPt_sig.png");


  TCanvas *c_bdtout_compareCiC_lowPt = new TCanvas("c_bdtout_compareCiC_lowPt","BDT output: di-photon pT<40 GeV, pass CiC supertight",1000,650);
  c_bdtout_compareCiC_lowPt->Divide(2,2);

  bdtout_passCiC_lowPt_cat0_tot->SetFillColor(2);
  bdtout_passCiC_lowPt_cat1_tot->SetFillColor(3);
  bdtout_passCiC_lowPt_cat2_tot->SetFillColor(4);
  bdtout_passCiC_lowPt_cat3_tot->SetFillColor(6);

  bdtout_passCiC_lowPt_cat0_tot->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat1_tot->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat2_tot->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat3_tot->SetFillStyle(3002);

  bdtout_lowPt_cat0_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat1_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat2_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat3_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat0_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat1_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat2_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat3_tot->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC_lowPt->cd(1);

  box_passCiC_lowPt_cat0 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat0->SetY2(bdtout_lowPt_cat0_tot->GetMaximum()*1.05);

  bdtout_lowPt_cat0_tot->Draw("hist");
  box_passCiC_lowPt_cat0->Draw("hist,same");
  bdtout_lowPt_cat0_tot->Draw("hist,same");
  bdtout_passCiC_lowPt_cat0_tot->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat0[i]->SetY2(bdtout_lowPt_cat0_tot->GetMaximum()*1.05);
    line_passCiC_lowPt_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(2);

  box_passCiC_lowPt_cat1 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat1->SetY2(bdtout_lowPt_cat1_tot->GetMaximum()*1.05);

  bdtout_lowPt_cat1_tot->Draw("hist");
  box_passCiC_lowPt_cat1->Draw("hist,same");
  bdtout_lowPt_cat1_tot->Draw("hist,same");
  bdtout_passCiC_lowPt_cat1_tot->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat1[i]->SetY2(bdtout_lowPt_cat1_tot->GetMaximum()*1.05);
    line_passCiC_lowPt_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(3);

  box_passCiC_lowPt_cat2 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat2->SetY2(bdtout_lowPt_cat2_tot->GetMaximum()*1.05);

  bdtout_lowPt_cat2_tot->Draw("hist");
  box_passCiC_lowPt_cat2->Draw("hist,same");
  bdtout_lowPt_cat2_tot->Draw("hist,same");
  bdtout_passCiC_lowPt_cat2_tot->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat2[i]->SetY2(bdtout_lowPt_cat2_tot->GetMaximum()*1.05);
    line_passCiC_lowPt_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(4);

  box_passCiC_lowPt_cat3 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat3->SetY2(bdtout_lowPt_cat3_tot->GetMaximum()*1.05);

  bdtout_lowPt_cat3_tot->Draw("hist");
  box_passCiC_lowPt_cat3->Draw("hist,same");
  bdtout_lowPt_cat3_tot->Draw("hist,same");
  bdtout_passCiC_lowPt_cat3_tot->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat3[i]->SetY2(bdtout_lowPt_cat3_tot->GetMaximum()*1.05);
    line_passCiC_lowPt_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->SaveAs("categories_compareCiC_lowPt_sig.png");


  TCanvas *c_bdtout_highPt = new TCanvas("c_bdtout_highPt","BDT output: di-photon pT>40 GeV");

  bdtout_highPt_all_tot = (TH1*)bdtout_highPt_cat0_tot->Clone();
  bdtout_highPt_all_tot->Add(bdtout_highPt_cat1_tot);
  bdtout_highPt_all_tot->Add(bdtout_highPt_cat2_tot);
  bdtout_highPt_all_tot->Add(bdtout_highPt_cat3_tot);

  bdtout_highPt_all_tot->SetLineColor(1);
  bdtout_highPt_cat0_tot->SetLineColor(2);
  bdtout_highPt_cat1_tot->SetLineColor(3);
  bdtout_highPt_cat2_tot->SetLineColor(4);
  bdtout_highPt_cat3_tot->SetLineColor(6);

  bdtout_highPt_all_tot->SetLineWidth(2);
  bdtout_highPt_cat0_tot->SetLineWidth(2);
  bdtout_highPt_cat1_tot->SetLineWidth(2);
  bdtout_highPt_cat2_tot->SetLineWidth(2);
  bdtout_highPt_cat3_tot->SetLineWidth(2);

  bdtout_highPt_all_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_all_tot->GetYaxis()->SetTitle("Events/0.02");
  bdtout_highPt_all_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_all_tot->GetXaxis()->SetLabelSize(0.05);
  bdtout_highPt_all_tot->GetYaxis()->SetTitleSize(0.05);
  bdtout_highPt_all_tot->GetYaxis()->SetLabelSize(0.05);

  max = 1.2*bdtout_highPt_all_tot->GetMaximum();
  bdtout_highPt_all_tot->SetMaximum(max);

  TBox* box_highPt = new TBox(-1.,0.,boundaries[0],max);
  box_highPt->SetFillColor(38);
  box_highPt->SetFillStyle(3002);

  bdtout_highPt_all_tot->Draw("hist");
  box_highPt->Draw("hist,same");
  bdtout_highPt_all_tot->Draw("hist,same");
  bdtout_highPt_cat0_tot->Draw("hist,same");
  bdtout_highPt_cat1_tot->Draw("hist,same");
  bdtout_highPt_cat2_tot->Draw("hist,same");
  bdtout_highPt_cat3_tot->Draw("hist,same");

  TLegend *leg_highPt = (TLegend*)leg->Clone();
  leg_highPt->Clear();
  leg_highPt->AddEntry(bdtout_all_tot,"All (p_{T}(#gamma#gamma) > 40 GeV)");
  leg_highPt->AddEntry(bdtout_cat0_tot,"both EB, both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat1_tot,"both EB, !both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat2_tot,"!both EB, both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat3_tot,"!both EB, !both R9>0.94");
  leg_highPt->Draw("hist");

  TLine* line_highPt[4];
  for (int i=0; i<4; i++) {
    line_highPt[i] = new TLine(boundaries[i],0.,boundaries[i],max);
    line_highPt[i]->SetLineColor(4);
    line_highPt[i]->SetLineWidth(2);
    line_highPt[i]->SetLineStyle(9);
    line_highPt[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout_highPt->SaveAs("categories_highPt_sig.png");



  TCanvas *c_bdtout_compareCiC_highPt = new TCanvas("c_bdtout_compareCiC_highPt","BDT output: di-photon pT>40 GeV, pass CiC supertight",1000,650);
  c_bdtout_compareCiC_highPt->Divide(2,2);

  bdtout_passCiC_highPt_cat0_tot->SetFillColor(2);
  bdtout_passCiC_highPt_cat1_tot->SetFillColor(3);
  bdtout_passCiC_highPt_cat2_tot->SetFillColor(4);
  bdtout_passCiC_highPt_cat3_tot->SetFillColor(6);

  bdtout_passCiC_highPt_cat0_tot->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat1_tot->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat2_tot->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat3_tot->SetFillStyle(3002);

  bdtout_highPt_cat0_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat1_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat2_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat3_tot->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat0_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat1_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat2_tot->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat3_tot->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC_highPt->cd(1);

  box_passCiC_highPt_cat0 = (TBox*)box->Clone();
  box_passCiC_highPt_cat0->SetY2(bdtout_highPt_cat0_tot->GetMaximum()*1.05);

  bdtout_highPt_cat0_tot->Draw("hist");
  box_passCiC_highPt_cat0->Draw("hist,same");
  bdtout_highPt_cat0_tot->Draw("hist,same");
  bdtout_passCiC_highPt_cat0_tot->Draw("hist,same");

  TLine* line_passCiC_highPt_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat0[i]->SetY2(bdtout_highPt_cat0_tot->GetMaximum()*1.05);
    line_passCiC_highPt_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(2);

  box_passCiC_highPt_cat1 = (TBox*)box->Clone();
  box_passCiC_highPt_cat1->SetY2(bdtout_highPt_cat1_tot->GetMaximum()*1.05);

  bdtout_highPt_cat1_tot->Draw("hist");
  box_passCiC_highPt_cat1->Draw("hist,same");
  bdtout_highPt_cat1_tot->Draw("hist,same");
  bdtout_passCiC_highPt_cat1_tot->Draw("hist,same");

  TLine* line_passCiC_highPt_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat1[i]->SetY2(bdtout_highPt_cat1_tot->GetMaximum()*1.05);
    line_passCiC_highPt_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(3);

  box_passCiC_highPt_cat2 = (TBox*)box->Clone();
  box_passCiC_highPt_cat2->SetY2(bdtout_highPt_cat2_tot->GetMaximum()*1.05);

  bdtout_highPt_cat2_tot->Draw("hist");
  box_passCiC_highPt_cat2->Draw("hist,same");
  bdtout_highPt_cat2_tot->Draw("hist,same");
  bdtout_passCiC_highPt_cat2_tot->Draw("hist,same");

  TLine* line_passCiC_highPt_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat2[i]->SetY2(bdtout_highPt_cat2_tot->GetMaximum()*1.05);
    line_passCiC_highPt_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(4);

  box_passCiC_highPt_cat3 = (TBox*)box->Clone();
  box_passCiC_highPt_cat3->SetY2(bdtout_highPt_cat3_tot->GetMaximum()*1.05);

  bdtout_highPt_cat3_tot->Draw("hist");
  box_passCiC_highPt_cat3->Draw("hist,same");
  bdtout_highPt_cat3_tot->Draw("hist,same");
  bdtout_passCiC_highPt_cat3_tot->Draw("hist,same");

  TLine* line_passCiC_highPt_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat3[i]->SetY2(bdtout_highPt_cat3_tot->GetMaximum()*1.05);
    line_passCiC_highPt_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->SaveAs("categories_compareCiC_highPt_sig.png");


  TCanvas *c_ptVsBdtout = new TCanvas("c_ptVsBdtout","pT(#gamma#gamma) vs MVA output",1000,650);
  c_ptVsBdtout->Divide(2,2);

  pt_vs_bdtout_cat0_all = (TH2*)pt_vs_bdtout_cat0_ggh_m125_8TeV->Clone();
  pt_vs_bdtout_cat1_all = (TH2*)pt_vs_bdtout_cat1_ggh_m125_8TeV->Clone();
  pt_vs_bdtout_cat2_all = (TH2*)pt_vs_bdtout_cat2_ggh_m125_8TeV->Clone();
  pt_vs_bdtout_cat3_all = (TH2*)pt_vs_bdtout_cat3_ggh_m125_8TeV->Clone();
  pt_vs_bdtout_cat4_all = (TH2*)pt_vs_bdtout_cat4_ggh_m125_8TeV->Clone();

  pt_vs_bdtout_cat0_all->Add(pt_vs_bdtout_cat0_vbf_m125_8TeV);
  pt_vs_bdtout_cat1_all->Add(pt_vs_bdtout_cat1_vbf_m125_8TeV);
  pt_vs_bdtout_cat2_all->Add(pt_vs_bdtout_cat2_vbf_m125_8TeV);
  pt_vs_bdtout_cat3_all->Add(pt_vs_bdtout_cat3_vbf_m125_8TeV);
  pt_vs_bdtout_cat4_all->Add(pt_vs_bdtout_cat4_vbf_m125_8TeV);

  pt_vs_bdtout_cat0_all->Add(pt_vs_bdtout_cat0_wzh_m125_8TeV);
  pt_vs_bdtout_cat1_all->Add(pt_vs_bdtout_cat1_wzh_m125_8TeV);
  pt_vs_bdtout_cat2_all->Add(pt_vs_bdtout_cat2_wzh_m125_8TeV);
  pt_vs_bdtout_cat3_all->Add(pt_vs_bdtout_cat3_wzh_m125_8TeV);
  pt_vs_bdtout_cat4_all->Add(pt_vs_bdtout_cat4_wzh_m125_8TeV);

  pt_vs_bdtout_cat0_all->Add(pt_vs_bdtout_cat0_tth_m125_8TeV);
  pt_vs_bdtout_cat1_all->Add(pt_vs_bdtout_cat1_tth_m125_8TeV);
  pt_vs_bdtout_cat2_all->Add(pt_vs_bdtout_cat2_tth_m125_8TeV);
  pt_vs_bdtout_cat3_all->Add(pt_vs_bdtout_cat3_tth_m125_8TeV);
  pt_vs_bdtout_cat4_all->Add(pt_vs_bdtout_cat4_tth_m125_8TeV);

  pt_vs_bdtout_cat0_all->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat1_all->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat2_all->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat3_all->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat0_all->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat1_all->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat2_all->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat3_all->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat0_all->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat1_all->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat2_all->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat3_all->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat0_all->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat1_all->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat2_all->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat3_all->GetYaxis()->SetTitleSize(0.05);

  c_ptVsBdtout->cd(1);
  pt_vs_bdtout_cat0_all->Draw("colz");
  box_ptVsBdtout = (TBox*)box->Clone();
  box_ptVsBdtout->SetY2(200.);
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat0_all->Draw("colz,same");
  TLine* line_ptVsBdtout[4];
  for (int i=0; i<4; i++) {
    line_ptVsBdtout[i] = (TLine*)line[i]->Clone();
    line_ptVsBdtout[i]->SetY2(200.);
    line_ptVsBdtout[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(2);
  pt_vs_bdtout_cat1_all->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat1_all->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(3);
  pt_vs_bdtout_cat2_all->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat2_all->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(4);
  pt_vs_bdtout_cat3_all->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat3_all->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_ptVsBdtout->SaveAs("ptVsBdtout_sig.png");


  TCanvas *c_2D = new TCanvas("c_2D","min(R9) vs max(eta), sublead eta vs lead eta",1250,500);
  c_2D->Divide(4,2);

  minR9_vs_maxEta_cat2_all = (TH2*)minR9_vs_maxEta_cat2_ggh_m125_8TeV->Clone();
  minR9_vs_maxEta_cat3_all = (TH2*)minR9_vs_maxEta_cat3_ggh_m125_8TeV->Clone();
  minR9_vs_maxEta_cat4_all = (TH2*)minR9_vs_maxEta_cat4_ggh_m125_8TeV->Clone();
  minR9_vs_maxEta_cat5_all = (TH2*)minR9_vs_maxEta_cat5_ggh_m125_8TeV->Clone();
  minR9_vs_maxEta_cat6_all = (TH2*)minR9_vs_maxEta_cat6_ggh_m125_8TeV->Clone();

  minR9_vs_maxEta_cat2_all->Add(minR9_vs_maxEta_cat2_vbf_m125_8TeV);
  minR9_vs_maxEta_cat3_all->Add(minR9_vs_maxEta_cat3_vbf_m125_8TeV);
  minR9_vs_maxEta_cat4_all->Add(minR9_vs_maxEta_cat4_vbf_m125_8TeV);
  minR9_vs_maxEta_cat5_all->Add(minR9_vs_maxEta_cat5_vbf_m125_8TeV);
  minR9_vs_maxEta_cat6_all->Add(minR9_vs_maxEta_cat6_vbf_m125_8TeV);

  minR9_vs_maxEta_cat2_all->Add(minR9_vs_maxEta_cat2_wzh_m125_8TeV);
  minR9_vs_maxEta_cat3_all->Add(minR9_vs_maxEta_cat3_wzh_m125_8TeV);
  minR9_vs_maxEta_cat4_all->Add(minR9_vs_maxEta_cat4_wzh_m125_8TeV);
  minR9_vs_maxEta_cat5_all->Add(minR9_vs_maxEta_cat5_wzh_m125_8TeV);
  minR9_vs_maxEta_cat6_all->Add(minR9_vs_maxEta_cat6_wzh_m125_8TeV);

  minR9_vs_maxEta_cat2_all->Add(minR9_vs_maxEta_cat2_tth_m125_8TeV);
  minR9_vs_maxEta_cat3_all->Add(minR9_vs_maxEta_cat3_tth_m125_8TeV);
  minR9_vs_maxEta_cat4_all->Add(minR9_vs_maxEta_cat4_tth_m125_8TeV);
  minR9_vs_maxEta_cat5_all->Add(minR9_vs_maxEta_cat5_tth_m125_8TeV);
  minR9_vs_maxEta_cat6_all->Add(minR9_vs_maxEta_cat6_tth_m125_8TeV);

  eta2_vs_eta1_cat2_all = (TH2*)eta2_vs_eta1_cat2_ggh_m125_8TeV->Clone();
  eta2_vs_eta1_cat3_all = (TH2*)eta2_vs_eta1_cat3_ggh_m125_8TeV->Clone();
  eta2_vs_eta1_cat4_all = (TH2*)eta2_vs_eta1_cat4_ggh_m125_8TeV->Clone();
  eta2_vs_eta1_cat5_all = (TH2*)eta2_vs_eta1_cat5_ggh_m125_8TeV->Clone();
  eta2_vs_eta1_cat6_all = (TH2*)eta2_vs_eta1_cat6_ggh_m125_8TeV->Clone();

  eta2_vs_eta1_cat2_all->Add(eta2_vs_eta1_cat2_vbf_m125_8TeV);
  eta2_vs_eta1_cat3_all->Add(eta2_vs_eta1_cat3_vbf_m125_8TeV);
  eta2_vs_eta1_cat4_all->Add(eta2_vs_eta1_cat4_vbf_m125_8TeV);
  eta2_vs_eta1_cat5_all->Add(eta2_vs_eta1_cat5_vbf_m125_8TeV);
  eta2_vs_eta1_cat6_all->Add(eta2_vs_eta1_cat6_vbf_m125_8TeV);

  eta2_vs_eta1_cat2_all->Add(eta2_vs_eta1_cat2_wzh_m125_8TeV);
  eta2_vs_eta1_cat3_all->Add(eta2_vs_eta1_cat3_wzh_m125_8TeV);
  eta2_vs_eta1_cat4_all->Add(eta2_vs_eta1_cat4_wzh_m125_8TeV);
  eta2_vs_eta1_cat5_all->Add(eta2_vs_eta1_cat5_wzh_m125_8TeV);
  eta2_vs_eta1_cat6_all->Add(eta2_vs_eta1_cat6_wzh_m125_8TeV);

  eta2_vs_eta1_cat2_all->Add(eta2_vs_eta1_cat2_tth_m125_8TeV);
  eta2_vs_eta1_cat3_all->Add(eta2_vs_eta1_cat3_tth_m125_8TeV);
  eta2_vs_eta1_cat4_all->Add(eta2_vs_eta1_cat4_tth_m125_8TeV);
  eta2_vs_eta1_cat5_all->Add(eta2_vs_eta1_cat5_tth_m125_8TeV);
  eta2_vs_eta1_cat6_all->Add(eta2_vs_eta1_cat6_tth_m125_8TeV);

  minR9_vs_maxEta_cat2_all->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat3_all->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat4_all->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat5_all->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat6_all->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat2_all->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat3_all->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat4_all->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat5_all->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat6_all->GetYaxis()->SetTitle("min(R9)");

  minR9_vs_maxEta_cat2_all->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_all->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_all->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_all->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_all->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat2_all->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_all->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_all->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_all->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_all->GetYaxis()->SetTitleSize(0.05);

  line_minr9 = new TLine(0.,0.94,2.5,0.94);
  line_maxeta = new TLine(1.479,0.5,1.479,1.);
  //line_minr9->SetLineColor(4);
  line_minr9->SetLineWidth(2);
  line_minr9->SetLineStyle(9);
  //line_maxeta->SetLineColor(4);
  line_maxeta->SetLineWidth(2);
  line_maxeta->SetLineStyle(9);

  c_2D->cd(1);
  minR9_vs_maxEta_cat2_all->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(2);
  minR9_vs_maxEta_cat3_all->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(3);
  minR9_vs_maxEta_cat4_all->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(4);
  minR9_vs_maxEta_cat5_all->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  //c_2D->cd(5);
  //minR9_vs_maxEta_cat6_all->Draw("colz");
  //line_minr9->Draw("hist");
  //line_maxeta->Draw("hist");

  eta2_vs_eta1_cat2_all->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat3_all->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat4_all->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat5_all->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat6_all->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat2_all->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat3_all->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat4_all->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat5_all->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat6_all->GetXaxis()->SetTitle("lead #eta");

  eta2_vs_eta1_cat2_all->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat3_all->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat4_all->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat5_all->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat6_all->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat2_all->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat3_all->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat4_all->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat5_all->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat6_all->GetYaxis()->SetTitleSize(0.05);

  c_2D->cd(5);
  eta2_vs_eta1_cat2_all->Draw("colz");
  c_2D->cd(6);
  eta2_vs_eta1_cat3_all->Draw("colz");
  c_2D->cd(7);
  eta2_vs_eta1_cat4_all->Draw("colz");
  c_2D->cd(8);
  eta2_vs_eta1_cat5_all->Draw("colz");
  //c_2D->cd(10);
  //eta2_vs_eta1_cat6_all->Draw("colz");

  c_2D->SaveAs("categories_2D_sig.png");

  /*
  TCanvas *c_2D_v2 = new TCanvas("c_2D_v2","min(R9) vs max(eta), |lead eta| vs |deltaEeta|",1500,500);
  c_2D_v2->Divide(5,2);

  minR9_vs_maxEta_cat2_tot->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat3_tot->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat4_tot->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat5_tot->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat6_tot->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat2_tot->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat3_tot->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat4_tot->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat5_tot->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat6_tot->GetYaxis()->SetTitle("min(R9)");

  minR9_vs_maxEta_cat2_tot->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_tot->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_tot->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_tot->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_tot->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat2_tot->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_tot->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_tot->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_tot->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_tot->GetYaxis()->SetTitleSize(0.05);

  c_2D_v2->cd(1);
  minR9_vs_maxEta_cat2_tot->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(2);
  minR9_vs_maxEta_cat3_tot->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(3);
  minR9_vs_maxEta_cat4_tot->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(4);
  minR9_vs_maxEta_cat5_tot->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(5);
  minR9_vs_maxEta_cat6_tot->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");

  eta1_vs_deltaEta_cat2_tot->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat3_tot->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat4_tot->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat5_tot->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat6_tot->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat2_tot->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat3_tot->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat4_tot->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat5_tot->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat6_tot->GetXaxis()->SetTitle("|#Delta#eta|");

  eta1_vs_deltaEta_cat2_tot->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat3_tot->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat4_tot->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat5_tot->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat6_tot->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat2_tot->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat3_tot->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat4_tot->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat5_tot->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat6_tot->GetYaxis()->SetTitleSize(0.05);

  c_2D_v2->cd(6);
  eta1_vs_deltaEta_cat2_tot->Draw("colz");
  c_2D_v2->cd(7);
  eta1_vs_deltaEta_cat3_tot->Draw("colz");
  c_2D_v2->cd(8);
  eta1_vs_deltaEta_cat4_tot->Draw("colz");
  c_2D_v2->cd(9);
  eta1_vs_deltaEta_cat5_tot->Draw("colz");
  c_2D_v2->cd(10);
  eta1_vs_deltaEta_cat6_tot->Draw("colz");

  c_2D_v2->SaveAs("categories_2D_v2_sig.png");
  */
}

float frac_err(float a, float b) {
  float a_err = sqrt(a);
  float b_err = sqrt(b);
  return ((1./(b*b)) * sqrt(b*b*a_err*a_err + a*a*b_err*b_err));
}
