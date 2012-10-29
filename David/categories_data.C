void categories_data() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TText *text = new TText();
  text->SetNDC();
  text->SetTextSize(0.05);

  TFile *file = TFile::Open("histograms_CMS-HGG_categories.root");
  file->cd();

  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output");
  c_bdtout->SetLogy();

  bdtout_all_Data = (TH1*)bdtout_cat0_Data->Clone();

  bdtout_all_Data->Add(bdtout_cat1_Data);
  bdtout_all_Data->Add(bdtout_cat2_Data);
  bdtout_all_Data->Add(bdtout_cat3_Data);

  cout << bdtout_all_Data->Integral() << endl;
  cout << bdtout_all_Data->Integral(49,100)+bdtout_all_Data->GetBinContent(48)/2. << endl;
  cout << bdtout_all_Data->GetBinLowEdge(49) << endl;
  cout << bdtout_all_Data->GetBinLowEdge(48) << endl;
  cout << (bdtout_all_Data->Integral(49,100)+bdtout_all_Data->GetBinContent(48)/2.)/bdtout_all_Data->Integral() << endl;

  bdtout_passCiC_all_Data = (TH1*)bdtout_passCiC_cat0_Data->Clone();
  bdtout_passCiC_all_Data->Add(bdtout_passCiC_cat1_Data);
  bdtout_passCiC_all_Data->Add(bdtout_passCiC_cat2_Data);
  bdtout_passCiC_all_Data->Add(bdtout_passCiC_cat3_Data);

  bdtout_failpresel_cat0_Data->Rebin(2);
  bdtout_failpresel_cat1_Data->Rebin(2);
  bdtout_failpresel_cat2_Data->Rebin(2);
  bdtout_failpresel_cat3_Data->Rebin(2);

  bdtout_failpresel_all_Data = (TH1*)bdtout_failpresel_cat0_Data->Clone();
  bdtout_failpresel_all_Data->Add(bdtout_failpresel_cat1_Data);
  bdtout_failpresel_all_Data->Add(bdtout_failpresel_cat2_Data);
  bdtout_failpresel_all_Data->Add(bdtout_failpresel_cat3_Data);

  float nFail_m100180_bdtout005 = bdtout_all_Data->Integral(1,52) + bdtout_all_Data->GetBinContent(53)/2.;
  float nFail_m100180 = bdtout_all_Data->Integral(1,100);
  float frac_bdtout005 = nFail_m100180_bdtout005/nFail_m100180;
  float frac_bdtout005_err = frac_err(nFail_m100180_bdtout005,nFail_m100180);
  cout << bdtout_all_Data->Integral(1,52) << " " <<  bdtout_all_Data->GetBinContent(53)/2. << " " << bdtout_all_Data->Integral(1,100) << " " << bdtout_all_Data->GetBinLowEdge(53) << " " << frac_bdtout005 << " ± " << frac_bdtout005_err << endl;

  bdtout_all_Data->SetLineColor(1);
  bdtout_cat0_Data->SetLineColor(2);
  bdtout_cat1_Data->SetLineColor(3);
  bdtout_cat2_Data->SetLineColor(4);
  bdtout_cat3_Data->SetLineColor(6);

  bdtout_all_Data->SetLineWidth(2);
  bdtout_cat0_Data->SetLineWidth(2);
  bdtout_cat1_Data->SetLineWidth(2);
  bdtout_cat2_Data->SetLineWidth(2);
  bdtout_cat3_Data->SetLineWidth(2);

  bdtout_cat0_Data->SetMarkerColor(2);
  bdtout_cat1_Data->SetMarkerColor(3);
  bdtout_cat2_Data->SetMarkerColor(4);
  bdtout_cat3_Data->SetMarkerColor(6);

  bdtout_all_Data->SetMarkerStyle(20);
  bdtout_cat0_Data->SetMarkerStyle(20);
  bdtout_cat1_Data->SetMarkerStyle(20);
  bdtout_cat2_Data->SetMarkerStyle(20);
  bdtout_cat3_Data->SetMarkerStyle(20);

  bdtout_all_Data->SetMarkerSize(0.7);
  bdtout_cat0_Data->SetMarkerSize(0.7);
  bdtout_cat1_Data->SetMarkerSize(0.7);
  bdtout_cat2_Data->SetMarkerSize(0.7);
  bdtout_cat3_Data->SetMarkerSize(0.7);

  bdtout_passCiC_all_Data->SetLineColor(1);
  bdtout_passCiC_cat0_Data->SetLineColor(2);
  bdtout_passCiC_cat1_Data->SetLineColor(3);
  bdtout_passCiC_cat2_Data->SetLineColor(4);
  bdtout_passCiC_cat3_Data->SetLineColor(6);

  bdtout_passCiC_lowPt_cat0_Data->SetLineColor(2);
  bdtout_passCiC_lowPt_cat1_Data->SetLineColor(3);
  bdtout_passCiC_lowPt_cat2_Data->SetLineColor(4);
  bdtout_passCiC_lowPt_cat3_Data->SetLineColor(6);

  bdtout_passCiC_highPt_cat0_Data->SetLineColor(2);
  bdtout_passCiC_highPt_cat1_Data->SetLineColor(3);
  bdtout_passCiC_highPt_cat2_Data->SetLineColor(4);
  bdtout_passCiC_highPt_cat3_Data->SetLineColor(6);

  bdtout_failpresel_cat0_Data->SetLineColor(2);
  bdtout_failpresel_cat1_Data->SetLineColor(3);
  bdtout_failpresel_cat2_Data->SetLineColor(4);
  bdtout_failpresel_cat3_Data->SetLineColor(6);


  bdtout_all_Data->GetXaxis()->SetTitle("di-photon MVA output");

  float boundaries[4] = {-0.05,0.49,0.79,.91};
  float max = 3.*bdtout_all_Data->GetMaximum();
  bdtout_all_Data->SetMaximum(max);

  TBox* box = new TBox(-1.,0.,boundaries[0],max);
  box->SetFillColor(38);
  box->SetFillStyle(3002);

  bdtout_all_Data->Draw("hist");
  box->Draw("hist,same");
  bdtout_all_Data->Draw("hist,same");
  bdtout_cat0_Data->Draw("hist,same");
  bdtout_cat1_Data->Draw("hist,same");
  bdtout_cat2_Data->Draw("hist,same");
  bdtout_cat3_Data->Draw("hist,same");

  TLegend *leg;
  leg = new TLegend(.12,.6,.44,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(bdtout_all_Data,"All");
  leg->AddEntry(bdtout_cat0_Data,"both EB, both R9>0.94");
  leg->AddEntry(bdtout_cat1_Data,"both EB, !both R9>0.94");
  leg->AddEntry(bdtout_cat2_Data,"!both EB, both R9>0.94");
  leg->AddEntry(bdtout_cat3_Data,"!both EB, !both R9>0.94");

  TLegend *leg2;
  leg2 = new TLegend(.37,.63,.61,.89);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(bdtout_all_Data,"All");
  leg2->AddEntry(bdtout_cat0_Data,"both EB, both R9>0.94");
  leg2->AddEntry(bdtout_cat1_Data,"both EB, !both R9>0.94");
  leg2->AddEntry(bdtout_cat2_Data,"!both EB, both R9>0.94");
  leg2->AddEntry(bdtout_cat3_Data,"!both EB, !both R9>0.94");
  leg2->Draw("hist");

  TLegend *leg3;
  leg3 = new TLegend(.37,.11,.61,.37);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(10);
  leg3->SetTextSize(.035);
  leg3->AddEntry(bdtout_all_Data,"All");
  leg3->AddEntry(bdtout_cat0_Data,"both EB, both R9>0.94");
  leg3->AddEntry(bdtout_cat1_Data,"both EB, !both R9>0.94");
  leg3->AddEntry(bdtout_cat2_Data,"!both EB, both R9>0.94");
  leg3->AddEntry(bdtout_cat3_Data,"!both EB, !both R9>0.94");

  TLine* line[4];
  for (int i=0; i<4; i++) {
    line[i] = new TLine(boundaries[i],0.,boundaries[i],max*1.05);
    line[i]->SetLineColor(4);
    line[i]->SetLineWidth(2);
    line[i]->SetLineStyle(9);
    line[i]->Draw("hist");
  }

  gPad->RedrawAxis();

  c_bdtout->SaveAs("categories.png");


  TCanvas *c_bdtout_passCiC = new TCanvas("c_bdtout_passCiC","BDT output, pass CiC");

  bdtout_passCiC_cat0_Data_clone = (TH1*)bdtout_passCiC_cat0_Data->Clone();
  bdtout_passCiC_cat1_Data_clone = (TH1*)bdtout_passCiC_cat1_Data->Clone();
  bdtout_passCiC_cat2_Data_clone = (TH1*)bdtout_passCiC_cat2_Data->Clone();
  bdtout_passCiC_cat3_Data_clone = (TH1*)bdtout_passCiC_cat3_Data->Clone();

  bdtout_passCiC_all_Data->SetLineWidth(2);
  bdtout_passCiC_cat0_Data_clone->SetLineWidth(2);
  bdtout_passCiC_cat1_Data_clone->SetLineWidth(2);
  bdtout_passCiC_cat2_Data_clone->SetLineWidth(2);
  bdtout_passCiC_cat3_Data_clone->SetLineWidth(2);

  bdtout_passCiC_all_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_passCiC_all_Data->GetYaxis()->SetTitle("Events/0.02");
  bdtout_passCiC_all_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_passCiC_all_Data->GetXaxis()->SetLabelSize(0.05);
  bdtout_passCiC_all_Data->GetYaxis()->SetTitleSize(0.05);
  bdtout_passCiC_all_Data->GetYaxis()->SetLabelSize(0.05);

  max = bdtout_passCiC_all_Data->GetMaximum()*1.15;
  TBox* box_passCiC = new TBox(-1.,0.,boundaries[0],max*1.05);
  box_passCiC->SetFillColor(38);
  box_passCiC->SetFillStyle(3002);

  bdtout_passCiC_all_Data->SetMaximum(max*1.05);
  bdtout_passCiC_all_Data->Draw("hist");
  box_passCiC->Draw("hist,same");
  bdtout_passCiC_all_Data->Draw("hist,same");
  bdtout_passCiC_cat0_Data_clone->Draw("hist,same");
  bdtout_passCiC_cat1_Data_clone->Draw("hist,same");
  bdtout_passCiC_cat2_Data_clone->Draw("hist,same");
  bdtout_passCiC_cat3_Data_clone->Draw("hist,same");

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

  c_bdtout_passCiC->SaveAs("categories_passCiC.png");

  /*
  TCanvas *c_bdtout_failpresel = new TCanvas("c_bdtout_failpresel","BDT output, pass CiC, fail presel");

  bdtout_failpresel_cat0_Data_clone = (TH1*)bdtout_failpresel_cat0_Data->Clone();
  bdtout_failpresel_cat1_Data_clone = (TH1*)bdtout_failpresel_cat1_Data->Clone();
  bdtout_failpresel_cat2_Data_clone = (TH1*)bdtout_failpresel_cat2_Data->Clone();
  bdtout_failpresel_cat3_Data_clone = (TH1*)bdtout_failpresel_cat3_Data->Clone();

  bdtout_failpresel_all_Data->SetLineWidth(2);
  bdtout_failpresel_cat0_Data_clone->SetLineWidth(2);
  bdtout_failpresel_cat1_Data_clone->SetLineWidth(2);
  bdtout_failpresel_cat2_Data_clone->SetLineWidth(2);
  bdtout_failpresel_cat3_Data_clone->SetLineWidth(2);

  bdtout_failpresel_all_Data->GetXaxis()->SetTitle("di-photon MVA output");

  max = bdtout_failpresel_all_Data->GetMaximum();
  TBox* box_failpresel = new TBox(-1.,0.,boundaries[0],max*1.05);
  box_failpresel->SetFillColor(38);
  box_failpresel->SetFillStyle(3002);

  bdtout_failpresel_all_Data->Draw("hist");
  box_failpresel->Draw("hist,same");
  bdtout_failpresel_all_Data->Draw("hist,same");
  bdtout_failpresel_cat0_Data_clone->Draw("hist,same");
  bdtout_failpresel_cat1_Data_clone->Draw("hist,same");
  bdtout_failpresel_cat2_Data_clone->Draw("hist,same");
  bdtout_failpresel_cat3_Data_clone->Draw("hist,same");

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

  c_bdtout_failpresel->SaveAs("categories_failpresel.png");
  */

  TCanvas *c_bdtout_compareCiC = new TCanvas("c_bdtout_compareCiC","BDT output: pass CiC supertight",1000,650);
  c_bdtout_compareCiC->Divide(2,2);

  bdtout_passCiC_cat0_Data->SetFillColor(2);
  bdtout_passCiC_cat1_Data->SetFillColor(3);
  bdtout_passCiC_cat2_Data->SetFillColor(4);
  bdtout_passCiC_cat3_Data->SetFillColor(6);

  bdtout_passCiC_cat0_Data->SetFillStyle(3002);
  bdtout_passCiC_cat1_Data->SetFillStyle(3002);
  bdtout_passCiC_cat2_Data->SetFillStyle(3002);
  bdtout_passCiC_cat3_Data->SetFillStyle(3002);

  bdtout_cat0_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat1_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat2_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat3_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat3_Data->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC->cd(1);

  box_passCiC_cat0 = (TBox*)box->Clone();
  box_passCiC_cat0->SetY2(bdtout_cat0_Data->GetMaximum()*1.05);

  bdtout_cat0_Data->Draw("hist");
  box_passCiC_cat0->Draw("hist,same");
  bdtout_cat0_Data->Draw("hist,same");
  bdtout_passCiC_cat0_Data->Draw("hist,same");

  TLine* line_passCiC_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat0[i]->SetY2(bdtout_cat0_Data->GetMaximum()*1.05);
    line_passCiC_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(2);

  box_passCiC_cat1 = (TBox*)box->Clone();
  box_passCiC_cat1->SetY2(bdtout_cat1_Data->GetMaximum()*1.05);

  bdtout_cat1_Data->Draw("hist");
  box_passCiC_cat1->Draw("hist,same");
  bdtout_cat1_Data->Draw("hist,same");
  bdtout_passCiC_cat1_Data->Draw("hist,same");

  TLine* line_passCiC_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat1[i]->SetY2(bdtout_cat1_Data->GetMaximum()*1.05);
    line_passCiC_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(3);

  box_passCiC_cat2 = (TBox*)box->Clone();
  box_passCiC_cat2->SetY2(bdtout_cat2_Data->GetMaximum()*1.05);

  bdtout_cat2_Data->Draw("hist");
  box_passCiC_cat2->Draw("hist,same");
  bdtout_cat2_Data->Draw("hist,same");
  bdtout_passCiC_cat2_Data->Draw("hist,same");

  TLine* line_passCiC_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat2[i]->SetY2(bdtout_cat2_Data->GetMaximum()*1.05);
    line_passCiC_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->cd(4);

  box_passCiC_cat3 = (TBox*)box->Clone();
  box_passCiC_cat3->SetY2(bdtout_cat3_Data->GetMaximum()*1.05);

  bdtout_cat3_Data->Draw("hist");
  box_passCiC_cat3->Draw("hist,same");
  bdtout_cat3_Data->Draw("hist,same");
  bdtout_passCiC_cat3_Data->Draw("hist,same");

  TLine* line_passCiC_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_cat3[i]->SetY2(bdtout_cat3_Data->GetMaximum()*1.05);
    line_passCiC_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC->SaveAs("categories_compareCiC.png");



  TCanvas *c_bdtout_lowPt = new TCanvas("c_bdtout_lowPt","BDT output: di-photon pT<40 GeV");
  c_bdtout_lowPt->SetLogy();

  bdtout_lowPt_all_Data = (TH1*)bdtout_lowPt_cat0_Data->Clone();
  bdtout_lowPt_all_Data->Add(bdtout_lowPt_cat1_Data);
  bdtout_lowPt_all_Data->Add(bdtout_lowPt_cat2_Data);
  bdtout_lowPt_all_Data->Add(bdtout_lowPt_cat3_Data);

  bdtout_lowPt_all_Data->SetLineColor(1);
  bdtout_lowPt_cat0_Data->SetLineColor(2);
  bdtout_lowPt_cat1_Data->SetLineColor(3);
  bdtout_lowPt_cat2_Data->SetLineColor(4);
  bdtout_lowPt_cat3_Data->SetLineColor(6);

  bdtout_lowPt_all_Data->SetLineWidth(2);
  bdtout_lowPt_cat0_Data->SetLineWidth(2);
  bdtout_lowPt_cat1_Data->SetLineWidth(2);
  bdtout_lowPt_cat2_Data->SetLineWidth(2);
  bdtout_lowPt_cat3_Data->SetLineWidth(2);

  bdtout_lowPt_all_Data->GetXaxis()->SetTitle("di-photon MVA output");

  max = 5.*bdtout_lowPt_all_Data->GetMaximum();
  bdtout_lowPt_all_Data->SetMaximum(max);
  TBox* box_lowPt = new TBox(-1.,0.,boundaries[0],max);
  box_lowPt->SetFillColor(38);
  box_lowPt->SetFillStyle(3002);

  bdtout_lowPt_all_Data->Draw("hist");
  box_lowPt->Draw("hist,same");
  bdtout_lowPt_all_Data->Draw("hist,same");
  bdtout_lowPt_cat0_Data->Draw("hist,same");
  bdtout_lowPt_cat1_Data->Draw("hist,same");
  bdtout_lowPt_cat2_Data->Draw("hist,same");
  bdtout_lowPt_cat3_Data->Draw("hist,same");

  TLegend *leg_lowPt = (TLegend*)leg2->Clone();
  leg_lowPt->Clear();
  leg_lowPt->AddEntry(bdtout_all_Data,"All (p_{T}(#gamma#gamma) < 40 GeV)");
  leg_lowPt->AddEntry(bdtout_cat0_Data,"both EB, both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat1_Data,"both EB, !both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat2_Data,"!both EB, both R9>0.94");
  leg_lowPt->AddEntry(bdtout_cat3_Data,"!both EB, !both R9>0.94");
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

  c_bdtout_lowPt->SaveAs("categories_lowPt.png");


  TCanvas *c_bdtout_compareCiC_lowPt = new TCanvas("c_bdtout_compareCiC_lowPt","BDT output: di-photon pT<40 GeV, pass CiC supertight",1000,650);
  c_bdtout_compareCiC_lowPt->Divide(2,2);

  bdtout_passCiC_lowPt_cat0_Data->SetFillColor(2);
  bdtout_passCiC_lowPt_cat1_Data->SetFillColor(3);
  bdtout_passCiC_lowPt_cat2_Data->SetFillColor(4);
  bdtout_passCiC_lowPt_cat3_Data->SetFillColor(6);

  bdtout_passCiC_lowPt_cat0_Data->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat1_Data->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat2_Data->SetFillStyle(3002);
  bdtout_passCiC_lowPt_cat3_Data->SetFillStyle(3002);

  bdtout_lowPt_cat0_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat1_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat2_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat3_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_lowPt_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_lowPt_cat3_Data->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC_lowPt->cd(1);

  box_passCiC_lowPt_cat0 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat0->SetY2(bdtout_lowPt_cat0_Data->GetMaximum()*1.05);

  bdtout_lowPt_cat0_Data->Draw("hist");
  box_passCiC_lowPt_cat0->Draw("hist,same");
  bdtout_lowPt_cat0_Data->Draw("hist,same");
  bdtout_passCiC_lowPt_cat0_Data->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat0[i]->SetY2(bdtout_lowPt_cat0_Data->GetMaximum()*1.05);
    line_passCiC_lowPt_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(2);

  box_passCiC_lowPt_cat1 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat1->SetY2(bdtout_lowPt_cat1_Data->GetMaximum()*1.05);

  bdtout_lowPt_cat1_Data->Draw("hist");
  box_passCiC_lowPt_cat1->Draw("hist,same");
  bdtout_lowPt_cat1_Data->Draw("hist,same");
  bdtout_passCiC_lowPt_cat1_Data->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat1[i]->SetY2(bdtout_lowPt_cat1_Data->GetMaximum()*1.05);
    line_passCiC_lowPt_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(3);

  box_passCiC_lowPt_cat2 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat2->SetY2(bdtout_lowPt_cat2_Data->GetMaximum()*1.05);

  bdtout_lowPt_cat2_Data->Draw("hist");
  box_passCiC_lowPt_cat2->Draw("hist,same");
  bdtout_lowPt_cat2_Data->Draw("hist,same");
  bdtout_passCiC_lowPt_cat2_Data->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat2[i]->SetY2(bdtout_lowPt_cat2_Data->GetMaximum()*1.05);
    line_passCiC_lowPt_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->cd(4);

  box_passCiC_lowPt_cat3 = (TBox*)box->Clone();
  box_passCiC_lowPt_cat3->SetY2(bdtout_lowPt_cat3_Data->GetMaximum()*1.05);

  bdtout_lowPt_cat3_Data->Draw("hist");
  box_passCiC_lowPt_cat3->Draw("hist,same");
  bdtout_lowPt_cat3_Data->Draw("hist,same");
  bdtout_passCiC_lowPt_cat3_Data->Draw("hist,same");

  TLine* line_passCiC_lowPt_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_lowPt_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_lowPt_cat3[i]->SetY2(bdtout_lowPt_cat3_Data->GetMaximum()*1.05);
    line_passCiC_lowPt_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_lowPt->SaveAs("categories_compareCiC_lowPt.png");


  TCanvas *c_bdtout_highPt = new TCanvas("c_bdtout_highPt","BDT output: di-photon pT>40 GeV");
  c_bdtout_highPt->SetLogy();

  bdtout_highPt_all_Data = (TH1*)bdtout_highPt_cat0_Data->Clone();
  bdtout_highPt_all_Data->Add(bdtout_highPt_cat1_Data);
  bdtout_highPt_all_Data->Add(bdtout_highPt_cat2_Data);
  bdtout_highPt_all_Data->Add(bdtout_highPt_cat3_Data);

  bdtout_highPt_all_Data->SetLineColor(1);
  bdtout_highPt_cat0_Data->SetLineColor(2);
  bdtout_highPt_cat1_Data->SetLineColor(3);
  bdtout_highPt_cat2_Data->SetLineColor(4);
  bdtout_highPt_cat3_Data->SetLineColor(6);

  bdtout_highPt_all_Data->SetLineWidth(2);
  bdtout_highPt_cat0_Data->SetLineWidth(2);
  bdtout_highPt_cat1_Data->SetLineWidth(2);
  bdtout_highPt_cat2_Data->SetLineWidth(2);
  bdtout_highPt_cat3_Data->SetLineWidth(2);

  bdtout_highPt_all_Data->GetXaxis()->SetTitle("di-photon MVA output");

  max = 3.*bdtout_highPt_all_Data->GetMaximum();
  bdtout_highPt_all_Data->SetMaximum(max);

  TBox* box_highPt = new TBox(-1.,0.,boundaries[0],max);
  box_highPt->SetFillColor(38);
  box_highPt->SetFillStyle(3002);

  bdtout_highPt_all_Data->Draw("hist");
  box_highPt->Draw("hist,same");
  bdtout_highPt_all_Data->Draw("hist,same");
  bdtout_highPt_cat0_Data->Draw("hist,same");
  bdtout_highPt_cat1_Data->Draw("hist,same");
  bdtout_highPt_cat2_Data->Draw("hist,same");
  bdtout_highPt_cat3_Data->Draw("hist,same");

  TLegend *leg_highPt = (TLegend*)leg2->Clone();
  leg_highPt->Clear();
  leg_highPt->AddEntry(bdtout_all_Data,"All (p_{T}(#gamma#gamma) > 40 GeV)");
  leg_highPt->AddEntry(bdtout_cat0_Data,"both EB, both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat1_Data,"both EB, !both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat2_Data,"!both EB, both R9>0.94");
  leg_highPt->AddEntry(bdtout_cat3_Data,"!both EB, !both R9>0.94");
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

  c_bdtout_highPt->SaveAs("categories_highPt.png");



  TCanvas *c_bdtout_compareCiC_highPt = new TCanvas("c_bdtout_compareCiC_highPt","BDT output: di-photon pT>40 GeV, pass CiC supertight",1000,650);
  c_bdtout_compareCiC_highPt->Divide(2,2);

  bdtout_passCiC_highPt_cat0_Data->SetFillColor(2);
  bdtout_passCiC_highPt_cat1_Data->SetFillColor(3);
  bdtout_passCiC_highPt_cat2_Data->SetFillColor(4);
  bdtout_passCiC_highPt_cat3_Data->SetFillColor(6);

  bdtout_passCiC_highPt_cat0_Data->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat1_Data->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat2_Data->SetFillStyle(3002);
  bdtout_passCiC_highPt_cat3_Data->SetFillStyle(3002);

  bdtout_highPt_cat0_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat1_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat2_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat3_Data->GetXaxis()->SetTitle("di-photon MVA output");
  bdtout_highPt_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_highPt_cat3_Data->GetXaxis()->SetTitleSize(0.05);

  c_bdtout_compareCiC_highPt->cd(1);

  box_passCiC_highPt_cat0 = (TBox*)box->Clone();
  box_passCiC_highPt_cat0->SetY2(bdtout_highPt_cat0_Data->GetMaximum()*1.05);

  bdtout_highPt_cat0_Data->Draw("hist");
  box_passCiC_highPt_cat0->Draw("hist,same");
  bdtout_highPt_cat0_Data->Draw("hist,same");
  bdtout_passCiC_highPt_cat0_Data->Draw("hist,same");

  TLine* line_passCiC_highPt_cat0[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat0[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat0[i]->SetY2(bdtout_highPt_cat0_Data->GetMaximum()*1.05);
    line_passCiC_highPt_cat0[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(2);

  box_passCiC_highPt_cat1 = (TBox*)box->Clone();
  box_passCiC_highPt_cat1->SetY2(bdtout_highPt_cat1_Data->GetMaximum()*1.05);

  bdtout_highPt_cat1_Data->Draw("hist");
  box_passCiC_highPt_cat1->Draw("hist,same");
  bdtout_highPt_cat1_Data->Draw("hist,same");
  bdtout_passCiC_highPt_cat1_Data->Draw("hist,same");

  TLine* line_passCiC_highPt_cat1[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat1[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat1[i]->SetY2(bdtout_highPt_cat1_Data->GetMaximum()*1.05);
    line_passCiC_highPt_cat1[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(3);

  box_passCiC_highPt_cat2 = (TBox*)box->Clone();
  box_passCiC_highPt_cat2->SetY2(bdtout_highPt_cat2_Data->GetMaximum()*1.05);

  bdtout_highPt_cat2_Data->Draw("hist");
  box_passCiC_highPt_cat2->Draw("hist,same");
  bdtout_highPt_cat2_Data->Draw("hist,same");
  bdtout_passCiC_highPt_cat2_Data->Draw("hist,same");

  TLine* line_passCiC_highPt_cat2[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat2[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat2[i]->SetY2(bdtout_highPt_cat2_Data->GetMaximum()*1.05);
    line_passCiC_highPt_cat2[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->cd(4);

  box_passCiC_highPt_cat3 = (TBox*)box->Clone();
  box_passCiC_highPt_cat3->SetY2(bdtout_highPt_cat3_Data->GetMaximum()*1.05);

  bdtout_highPt_cat3_Data->Draw("hist");
  box_passCiC_highPt_cat3->Draw("hist,same");
  bdtout_highPt_cat3_Data->Draw("hist,same");
  bdtout_passCiC_highPt_cat3_Data->Draw("hist,same");

  TLine* line_passCiC_highPt_cat3[4];
  for (int i=0; i<4; i++) {
    line_passCiC_highPt_cat3[i] = (TLine*)line[i]->Clone();
    line_passCiC_highPt_cat3[i]->SetY2(bdtout_highPt_cat3_Data->GetMaximum()*1.05);
    line_passCiC_highPt_cat3[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_bdtout_compareCiC_highPt->SaveAs("categories_compareCiC_highPt.png");


  TCanvas *c_ptVsBdtout = new TCanvas("c_ptVsBdtout","pT(#gamma#gamma) vs MVA output",1000,650);
  c_ptVsBdtout->Divide(2,2);

  pt_vs_bdtout_cat0_Data->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat1_Data->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat2_Data->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat3_Data->GetXaxis()->SetTitle("di-photon MVA output");
  pt_vs_bdtout_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat0_Data->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat1_Data->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat2_Data->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat3_Data->GetYaxis()->SetTitle("di-photon p_{T} (GeV)");
  pt_vs_bdtout_cat0_Data->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat1_Data->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat2_Data->GetYaxis()->SetTitleSize(0.05);
  pt_vs_bdtout_cat3_Data->GetYaxis()->SetTitleSize(0.05);

  c_ptVsBdtout->cd(1);
  pt_vs_bdtout_cat0_Data->Draw("colz");
  box_ptVsBdtout = (TBox*)box->Clone();
  box_ptVsBdtout->SetY2(200.);
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat0_Data->Draw("colz,same");
  TLine* line_ptVsBdtout[4];
  for (int i=0; i<4; i++) {
    line_ptVsBdtout[i] = (TLine*)line[i]->Clone();
    line_ptVsBdtout[i]->SetY2(200.);
    line_ptVsBdtout[i]->Draw("hist");
  }
  text->DrawText(0.15,0.75,"both EB, both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(2);
  pt_vs_bdtout_cat1_Data->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat1_Data->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"both EB, !both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(3);
  pt_vs_bdtout_cat2_Data->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat2_Data->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"!both EB, both R9>0.94");
  gPad->RedrawAxis();
  c_ptVsBdtout->cd(4);
  pt_vs_bdtout_cat3_Data->Draw("colz");
  box_ptVsBdtout->Draw("hist");
  pt_vs_bdtout_cat3_Data->Draw("colz,same");
  for (int i=0; i<4; i++) line_ptVsBdtout[i]->Draw("hist");
  text->DrawText(0.15,0.75,"!both EB, !both R9>0.94");
  gPad->RedrawAxis();

  c_ptVsBdtout->SaveAs("ptVsBdtout.png");


  TCanvas *c_2D = new TCanvas("c_2D","min(R9) vs max(eta), sublead eta vs lead eta",1200,500);
  c_2D->Divide(4,2);

  minR9_vs_maxEta_cat2_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat3_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat4_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat5_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat6_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat2_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat3_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat4_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat5_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat6_Data->GetYaxis()->SetTitle("min(R9)");

  minR9_vs_maxEta_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat2_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_Data->GetYaxis()->SetTitleSize(0.05);

  line_minr9 = new TLine(0.,0.94,2.5,0.94);
  line_maxeta = new TLine(1.479,0.5,1.479,1.);
  //line_minr9->SetLineColor(4);
  line_minr9->SetLineWidth(2);
  line_minr9->SetLineStyle(9);
  //line_maxeta->SetLineColor(4);
  line_maxeta->SetLineWidth(2);
  line_maxeta->SetLineStyle(9);

  c_2D->cd(1);
  minR9_vs_maxEta_cat2_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(2);
  minR9_vs_maxEta_cat3_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(3);
  minR9_vs_maxEta_cat4_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D->cd(4);
  minR9_vs_maxEta_cat5_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  //c_2D->cd(5);
  //minR9_vs_maxEta_cat6_Data->Draw("colz");
  //line_minr9->Draw("hist");
  //line_maxeta->Draw("hist");

  eta2_vs_eta1_cat2_Data->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat3_Data->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat4_Data->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat5_Data->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat6_Data->GetYaxis()->SetTitle("sublead #eta");
  eta2_vs_eta1_cat2_Data->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat3_Data->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat4_Data->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat5_Data->GetXaxis()->SetTitle("lead #eta");
  eta2_vs_eta1_cat6_Data->GetXaxis()->SetTitle("lead #eta");

  eta2_vs_eta1_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat5_Data->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat6_Data->GetXaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat2_Data->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat3_Data->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat4_Data->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat5_Data->GetYaxis()->SetTitleSize(0.05);
  eta2_vs_eta1_cat6_Data->GetYaxis()->SetTitleSize(0.05);

  c_2D->cd(5);
  eta2_vs_eta1_cat2_Data->Draw("colz");
  c_2D->cd(6);
  eta2_vs_eta1_cat3_Data->Draw("colz");
  c_2D->cd(7);
  eta2_vs_eta1_cat4_Data->Draw("colz");
  c_2D->cd(8);
  eta2_vs_eta1_cat5_Data->Draw("colz");
  //c_2D->cd(10);
  //eta2_vs_eta1_cat6_Data->Draw("colz");

  c_2D->SaveAs("categories_2D.png");

  /*
  TCanvas *c_2D_v2 = new TCanvas("c_2D_v2","min(R9) vs max(eta), |lead eta| vs |deltaEeta|",1500,500);
  c_2D_v2->Divide(5,2);

  minR9_vs_maxEta_cat2_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat3_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat4_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat5_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat6_Data->GetXaxis()->SetTitle("max(#eta)");
  minR9_vs_maxEta_cat2_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat3_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat4_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat5_Data->GetYaxis()->SetTitle("min(R9)");
  minR9_vs_maxEta_cat6_Data->GetYaxis()->SetTitle("min(R9)");

  minR9_vs_maxEta_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_Data->GetXaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat2_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat3_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat4_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat5_Data->GetYaxis()->SetTitleSize(0.05);
  minR9_vs_maxEta_cat6_Data->GetYaxis()->SetTitleSize(0.05);

  c_2D_v2->cd(1);
  minR9_vs_maxEta_cat2_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(2);
  minR9_vs_maxEta_cat3_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(3);
  minR9_vs_maxEta_cat4_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(4);
  minR9_vs_maxEta_cat5_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");
  c_2D_v2->cd(5);
  minR9_vs_maxEta_cat6_Data->Draw("colz");
  line_minr9->Draw("hist");
  line_maxeta->Draw("hist");

  eta1_vs_deltaEta_cat2_Data->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat3_Data->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat4_Data->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat5_Data->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat6_Data->GetYaxis()->SetTitle("|lead #eta|");
  eta1_vs_deltaEta_cat2_Data->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat3_Data->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat4_Data->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat5_Data->GetXaxis()->SetTitle("|#Delta#eta|");
  eta1_vs_deltaEta_cat6_Data->GetXaxis()->SetTitle("|#Delta#eta|");

  eta1_vs_deltaEta_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat5_Data->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat6_Data->GetXaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat2_Data->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat3_Data->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat4_Data->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat5_Data->GetYaxis()->SetTitleSize(0.05);
  eta1_vs_deltaEta_cat6_Data->GetYaxis()->SetTitleSize(0.05);

  c_2D_v2->cd(6);
  eta1_vs_deltaEta_cat2_Data->Draw("colz");
  c_2D_v2->cd(7);
  eta1_vs_deltaEta_cat3_Data->Draw("colz");
  c_2D_v2->cd(8);
  eta1_vs_deltaEta_cat4_Data->Draw("colz");
  c_2D_v2->cd(9);
  eta1_vs_deltaEta_cat5_Data->Draw("colz");
  c_2D_v2->cd(10);
  eta1_vs_deltaEta_cat6_Data->Draw("colz");

  c_2D_v2->SaveAs("categories_2D_v2.png");
  */
}

float frac_err(float a, float b) {
  float a_err = sqrt(a);
  float b_err = sqrt(b);
  return ((1./(b*b)) * sqrt(b*b*a_err*a_err + a*a*b_err*b_err));
}
