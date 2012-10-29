void BDT_Zee_bdtin() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  /*
  TFile *file_hgg = TFile::Open("histograms_CMS-HGG_nominal_hgg.root");
  file_hgg->cd();

  bdtout_Hgg90 = (TH1*)bdtout_cat0_tot->Clone();
  bdtout_Hgg90->Add(bdtout_cat0_tot);
  bdtoutEB_Hgg90 = (TH1*)bdtoutEB_cat0_tot->Clone();
  bdtoutEB_Hgg90->Add(bdtoutEB_cat0_tot);
  bdtoutEBEE_Hgg90 = (TH1*)bdtoutEBEE_cat0_tot->Clone();
  bdtoutEBEE_Hgg90->Add(bdtoutEBEE_cat0_tot);
  bdtoutEE_Hgg90 = (TH1*)bdtoutEE_cat0_tot->Clone();
  bdtoutEE_Hgg90->Add(bdtoutEE_cat0_tot);

  bdtout_Hgg90->SetLineColor(3);
  bdtoutEB_Hgg90->SetLineColor(3);
  bdtoutEBEE_Hgg90->SetLineColor(3);
  bdtoutEE_Hgg90->SetLineColor(3);
  bdtout_Hgg90->SetLineWidth(2);
  bdtoutEB_Hgg90->SetLineWidth(2);
  bdtoutEBEE_Hgg90->SetLineWidth(2);
  bdtoutEE_Hgg90->SetLineWidth(2);
  bdtout_Hgg90->Rebin(2);
  bdtoutEB_Hgg90->Rebin(2);
  bdtoutEBEE_Hgg90->Rebin(2);
  bdtoutEE_Hgg90->Rebin(4);
  */

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
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
  bdtoutEE_cat0_DYJetsToLL->Rebin(4);
  bdtout_cat0_Data->Rebin(2);
  bdtout_cat1_Data->Rebin(2);
  bdtout_cat2_Data->Rebin(2);
  bdtout_cat3_Data->Rebin(2);
  bdtout_cat4_Data->Rebin(2);
  bdtoutEB_cat0_Data->Rebin(2);
  bdtoutEBEE_cat0_Data->Rebin(2);
  bdtoutEE_cat0_Data->Rebin(4);

  bdtout_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output");
  bdtout_cat1_Data->GetXaxis()->SetTitle("di-photon BDT output, baseline cat0");
  bdtout_cat2_Data->GetXaxis()->SetTitle("di-photon BDT output, baseline cat1");
  bdtout_cat3_Data->GetXaxis()->SetTitle("di-photon BDT output, baseline cat2");
  bdtout_cat4_Data->GetXaxis()->SetTitle("di-photon BDT output, baseline cat3");
  bdtoutEB_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EB)");
  bdtoutEBEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (one EB, one EE)");
  bdtoutEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EE)");
  pho1_phoidMva_cat0_Data->GetXaxis()->SetTitle("lead photon ID MVA output");
  pho2_phoidMva_cat0_Data->GetXaxis()->SetTitle("sublead photon ID MVA output");
  sigmaMOverM_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex)");
  sigmaMOverM_wrongVtx_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex)");
  sigmaMOverM_EB_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex)");
  sigmaMOverM_wrongVtx_EB_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex)");
  sigmaMOverM_EE_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex)");
  sigmaMOverM_wrongVtx_EE_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex)");
  vtxProb_cat0_Data->GetXaxis()->SetTitle("Vertex probabilty");
  pho1_ptOverM_cat0_Data->GetXaxis()->SetTitle("lead p_{T}/M_{#gamma#gamma}");
  pho2_ptOverM_cat0_Data->GetXaxis()->SetTitle("sublead p_{T}/M_{#gamma#gamma}");
  pho1_eta_cat0_Data->GetXaxis()->SetTitle("lead #eta");
  pho2_eta_cat0_Data->GetXaxis()->SetTitle("sublead #eta");
  cosDeltaPhi_cat0_Data->GetXaxis()->SetTitle("cos(#Delta#phi)");

  bdtout_cat0_Data->GetYaxis()->SetTitle("");
  bdtout_cat1_Data->GetYaxis()->SetTitle("");
  bdtout_cat2_Data->GetYaxis()->SetTitle("");
  bdtout_cat3_Data->GetYaxis()->SetTitle("");
  bdtout_cat4_Data->GetYaxis()->SetTitle("");
  bdtoutEB_cat0_Data->GetYaxis()->SetTitle("");
  bdtoutEBEE_cat0_Data->GetYaxis()->SetTitle("");
  bdtoutEE_cat0_Data->GetYaxis()->SetTitle("");
  pho1_phoidMva_cat0_Data->GetYaxis()->SetTitle("");
  pho2_phoidMva_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_wrongVtx_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_EB_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_wrongVtx_EB_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_EE_cat0_Data->GetYaxis()->SetTitle("");
  sigmaMOverM_wrongVtx_EE_cat0_Data->GetYaxis()->SetTitle("");
  vtxProb_cat0_Data->GetYaxis()->SetTitle("");
  pho1_ptOverM_cat0_Data->GetYaxis()->SetTitle("");
  pho2_ptOverM_cat0_Data->GetYaxis()->SetTitle("");
  pho1_eta_cat0_Data->GetYaxis()->SetTitle("");
  pho2_eta_cat0_Data->GetYaxis()->SetTitle("");
  cosDeltaPhi_cat0_Data->GetYaxis()->SetTitle("");

  bdtout_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat1_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat2_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat3_Data->GetXaxis()->SetTitleSize(0.05);
  bdtout_cat4_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEB_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEBEE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  bdtoutEE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho1_phoidMva_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho2_phoidMva_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_wrongVtx_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_EB_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_wrongVtx_EB_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  sigmaMOverM_wrongVtx_EE_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  vtxProb_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho1_ptOverM_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho2_ptOverM_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho1_eta_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  pho2_eta_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  cosDeltaPhi_cat0_Data->GetXaxis()->SetTitleSize(0.05);

  bdtout_cat0_DYJetsToLL->SetFillColor(38);
  bdtout_cat1_DYJetsToLL->SetFillColor(38);
  bdtout_cat2_DYJetsToLL->SetFillColor(38);
  bdtout_cat3_DYJetsToLL->SetFillColor(38);
  bdtout_cat4_DYJetsToLL->SetFillColor(38);
  bdtoutEB_cat0_DYJetsToLL->SetFillColor(38);
  bdtoutEBEE_cat0_DYJetsToLL->SetFillColor(38);
  bdtoutEE_cat0_DYJetsToLL->SetFillColor(38);
  pho1_phoidMva_cat0_DYJetsToLL->SetFillColor(38);
  pho2_phoidMva_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_wrongVtx_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_EB_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_EE_cat0_DYJetsToLL->SetFillColor(38);
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->SetFillColor(38);
  vtxProb_cat0_DYJetsToLL->SetFillColor(38);
  pho1_ptOverM_cat0_DYJetsToLL->SetFillColor(38);
  pho2_ptOverM_cat0_DYJetsToLL->SetFillColor(38);
  pho1_eta_cat0_DYJetsToLL->SetFillColor(38);
  pho2_eta_cat0_DYJetsToLL->SetFillColor(38);
  cosDeltaPhi_cat0_DYJetsToLL->SetFillColor(38);

  bdtout_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEB_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEBEE_cat0_DYJetsToLL->SetLineColor(1);
  bdtoutEE_cat0_DYJetsToLL->SetLineColor(1);
  pho1_phoidMva_cat0_DYJetsToLL->SetLineColor(1);
  pho2_phoidMva_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_wrongVtx_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_EB_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_EE_cat0_DYJetsToLL->SetLineColor(1);
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->SetLineColor(1);
  vtxProb_cat0_DYJetsToLL->SetLineColor(1);
  pho1_ptOverM_cat0_DYJetsToLL->SetLineColor(1);
  pho2_ptOverM_cat0_DYJetsToLL->SetLineColor(1);
  pho1_eta_cat0_DYJetsToLL->SetLineColor(1);
  pho2_eta_cat0_DYJetsToLL->SetLineColor(1);
  cosDeltaPhi_cat0_DYJetsToLL->SetLineColor(1);

  bdtout_cat0_Data->SetMarkerStyle(20);
  bdtout_cat1_Data->SetMarkerStyle(20);
  bdtout_cat2_Data->SetMarkerStyle(20);
  bdtout_cat3_Data->SetMarkerStyle(20);
  bdtout_cat4_Data->SetMarkerStyle(20);
  bdtoutEB_cat0_Data->SetMarkerStyle(20);
  bdtoutEBEE_cat0_Data->SetMarkerStyle(20);
  bdtoutEE_cat0_Data->SetMarkerStyle(20);
  pho1_phoidMva_cat0_Data->SetMarkerStyle(20);
  pho2_phoidMva_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_wrongVtx_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_EB_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_wrongVtx_EB_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_EE_cat0_Data->SetMarkerStyle(20);
  sigmaMOverM_wrongVtx_EE_cat0_Data->SetMarkerStyle(20);
  vtxProb_cat0_Data->SetMarkerStyle(20);
  pho1_ptOverM_cat0_Data->SetMarkerStyle(20);
  pho2_ptOverM_cat0_Data->SetMarkerStyle(20);
  pho1_eta_cat0_Data->SetMarkerStyle(20);
  pho2_eta_cat0_Data->SetMarkerStyle(20);
  cosDeltaPhi_cat0_Data->SetMarkerStyle(20);

  bdtout_cat0_Data->SetMarkerSize(0.4);
  bdtout_cat1_Data->SetMarkerSize(0.4);
  bdtout_cat2_Data->SetMarkerSize(0.4);
  bdtout_cat3_Data->SetMarkerSize(0.4);
  bdtout_cat4_Data->SetMarkerSize(0.4);
  bdtoutEB_cat0_Data->SetMarkerSize(0.4);
  bdtoutEBEE_cat0_Data->SetMarkerSize(0.4);
  bdtoutEE_cat0_Data->SetMarkerSize(0.4);
  pho1_phoidMva_cat0_Data->SetMarkerSize(0.4);
  pho2_phoidMva_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_wrongVtx_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_EB_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_wrongVtx_EB_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_EE_cat0_Data->SetMarkerSize(0.4);
  sigmaMOverM_wrongVtx_EE_cat0_Data->SetMarkerSize(0.4);
  vtxProb_cat0_Data->SetMarkerSize(0.4);
  pho1_ptOverM_cat0_Data->SetMarkerSize(0.4);
  pho2_ptOverM_cat0_Data->SetMarkerSize(0.4);
  pho1_eta_cat0_Data->SetMarkerSize(0.4);
  pho2_eta_cat0_Data->SetMarkerSize(0.4);
  cosDeltaPhi_cat0_Data->SetMarkerSize(0.4);

  bdtout_cat0_Data->SetLineColor(1);
  bdtout_cat1_Data->SetLineColor(1);
  bdtout_cat2_Data->SetLineColor(1);
  bdtout_cat3_Data->SetLineColor(1);
  bdtout_cat4_Data->SetLineColor(1);
  bdtoutEB_cat0_Data->SetLineColor(1);
  bdtoutEBEE_cat0_Data->SetLineColor(1);
  bdtoutEE_cat0_Data->SetLineColor(1);
  pho1_phoidMva_cat0_Data->SetLineColor(1);
  pho2_phoidMva_cat0_Data->SetLineColor(1);
  sigmaMOverM_cat0_Data->SetLineColor(1);
  sigmaMOverM_wrongVtx_cat0_Data->SetLineColor(1);
  sigmaMOverM_EB_cat0_Data->SetLineColor(1);
  sigmaMOverM_wrongVtx_EB_cat0_Data->SetLineColor(1);
  sigmaMOverM_EE_cat0_Data->SetLineColor(1);
  sigmaMOverM_wrongVtx_EE_cat0_Data->SetLineColor(1);
  vtxProb_cat0_Data->SetLineColor(1);
  pho1_ptOverM_cat0_Data->SetLineColor(1);
  pho2_ptOverM_cat0_Data->SetLineColor(1);
  pho1_eta_cat0_Data->SetLineColor(1);
  pho2_eta_cat0_Data->SetLineColor(1);
  cosDeltaPhi_cat0_Data->SetLineColor(1);


  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(vtxProb_cat0_Data,"Data (12.2fb^{-1})");
  leg->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(vtxProb_cat0_Data,"Data (12.2fb^{-1})");
  leg2->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");
  //leg2->AddEntry(bdtout_Hgg90,"H#rightarrow#gamma#gamma, m_{H}=90 GeV");

  TLegend *leg3;
  leg3 = new TLegend(.45,.65,.87,.87);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(10);
  leg3->SetTextSize(.03);
  leg3->AddEntry(vtxProb_cat0_Data,"Data (12.2fb^{-1})");
  leg3->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg4;
  leg4 = new TLegend(.2,.65,.52,.87);
  leg4->SetBorderSize(0);
  leg4->SetFillColor(10);
  leg4->SetTextSize(.025);
  leg4->AddEntry(vtxProb_cat0_Data,"Data (12.2fb^{-1})");
  leg4->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TCanvas *c_bdtout = new TCanvas("c_bdtout","BDT output",2200,800);
  c_bdtout->Divide(4,2);

  c_bdtout->cd(1);
  float sf = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  bdtout_cat0_DYJetsToLL->Scale(sf);
  bdtout_cat0_Data->Draw("e");
  bdtout_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_cat0_Data->Draw("e,same");
  //bdtout_Hgg90->Scale(bdtout_cat0_Data->Integral()/bdtout_Hgg90->Integral());
  //bdtout_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtout->cd(2);
  sf = bdtoutEB_cat0_Data->Integral()/bdtoutEB_cat0_DYJetsToLL->Integral();
  bdtoutEB_cat0_DYJetsToLL->Scale(sf);
  bdtoutEB_cat0_Data->Draw("e");
  bdtoutEB_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEB_cat0_Data->Draw("e,same");
  //bdtoutEB_Hgg90->Scale(bdtoutEB_cat0_Data->Integral()/bdtoutEB_Hgg90->Integral());
  //bdtoutEB_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtout->cd(3);
  sf = bdtoutEBEE_cat0_Data->Integral()/bdtoutEBEE_cat0_DYJetsToLL->Integral();
  bdtoutEBEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEBEE_cat0_Data->Draw("e");
  bdtoutEBEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEBEE_cat0_Data->Draw("e,same");
  //bdtoutEBEE_Hgg90->Scale(bdtoutEBEE_cat0_Data->Integral()/bdtoutEBEE_Hgg90->Integral());
  //bdtoutEBEE_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtout->cd(4);
  sf = bdtoutEE_cat0_Data->Integral()/bdtoutEE_cat0_DYJetsToLL->Integral();
  bdtoutEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEE_cat0_Data->Draw("e");
  bdtoutEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEE_cat0_Data->Draw("e,same");
  //bdtoutEE_Hgg90->Scale(bdtoutEE_cat0_Data->Integral()/bdtoutEE_Hgg90->Integral());
  //bdtoutEE_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg->Draw();

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);

  c_bdtout->cd(5);
  gPad->SetGrid();
  ratio_bdtout = (TH1*)bdtout_cat0_Data->Clone();
  ratio_bdtout->Sumw2();
  ratio_bdtout->Divide(bdtout_cat0_DYJetsToLL);
  ratio_bdtout->SetMaximum(1.4);
  ratio_bdtout->SetMinimum(0.6);
  ratio_bdtout->Draw("e");
  line->Draw();

  c_bdtout->cd(6);
  gPad->SetGrid();
  ratio_bdtoutEB = (TH1*)bdtoutEB_cat0_Data->Clone();
  ratio_bdtoutEB->Sumw2();
  ratio_bdtoutEB->Divide(bdtoutEB_cat0_DYJetsToLL);
  ratio_bdtoutEB->SetMaximum(1.4);
  ratio_bdtoutEB->SetMinimum(0.6);
  ratio_bdtoutEB->Draw("e");
  line->Draw();

  c_bdtout->cd(7);
  gPad->SetGrid();
  ratio_bdtoutEBEE = (TH1*)bdtoutEBEE_cat0_Data->Clone();
  ratio_bdtoutEBEE->Sumw2();
  ratio_bdtoutEBEE->Divide(bdtoutEBEE_cat0_DYJetsToLL);
  ratio_bdtoutEBEE->SetMaximum(1.8);
  ratio_bdtoutEBEE->SetMinimum(0.2);
  ratio_bdtoutEBEE->Draw("e");
  line->Draw();

  c_bdtout->cd(8);
  gPad->SetGrid();
  ratio_bdtoutEE = (TH1*)bdtoutEE_cat0_Data->Clone();
  ratio_bdtoutEE->Sumw2();
  ratio_bdtoutEE->Divide(bdtoutEE_cat0_DYJetsToLL);
  ratio_bdtoutEE->SetMaximum(1.8);
  ratio_bdtoutEE->SetMinimum(0.2);
  ratio_bdtoutEE->Draw("e");
  line->Draw();

  c_bdtout->SaveAs("bdtout.png");


  TCanvas *c_bdtout_basecat = new TCanvas("c_bdtout_basecat","BDT output in baseline categories",2200,800);
  c_bdtout_basecat->Divide(4,2);

  c_bdtout_basecat->cd(1);
  float sf = bdtout_cat1_Data->Integral()/bdtout_cat1_DYJetsToLL->Integral();
  bdtout_cat1_DYJetsToLL->Scale(sf);
  bdtout_cat1_Data->Draw("e");
  bdtout_cat1_DYJetsToLL->Draw("hist,same");
  bdtout_cat1_Data->Draw("e,same");
  //bdtout_cat1_Hgg90->Scale(bdtout_cat1_Data->Integral()/bdtout_cat1_Hgg90->Integral());
  //bdtout_cat1_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtout_basecat->cd(2);
  sf = bdtout_cat2_Data->Integral()/bdtout_cat2_DYJetsToLL->Integral();
  bdtout_cat2_DYJetsToLL->Scale(sf);
  bdtout_cat2_Data->Draw("e");
  bdtout_cat2_DYJetsToLL->Draw("hist,same");
  bdtout_cat2_Data->Draw("e,same");
  //bdtout_cat2_Hgg90->Scale(bdtout_cat2_Data->Integral()/bdtout_cat2_Hgg90->Integral());
  //bdtout_cat2_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtout_basecat->cd(3);
  sf = bdtout_cat3_Data->Integral()/bdtout_cat3_DYJetsToLL->Integral();
  bdtout_cat3_DYJetsToLL->Scale(sf);
  bdtout_cat3_Data->Draw("e");
  bdtout_cat3_DYJetsToLL->Draw("hist,same");
  bdtout_cat3_Data->Draw("e,same");
  //bdtout_cat3_Hgg90->Scale(bdtout_cat3_Data->Integral()/bdtout_cat3_Hgg90->Integral());
  //bdtout_cat3_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtout_basecat->cd(4);
  sf = bdtout_cat4_Data->Integral()/bdtout_cat4_DYJetsToLL->Integral();
  bdtout_cat4_DYJetsToLL->Scale(sf);
  bdtout_cat4_Data->Draw("e");
  bdtout_cat4_DYJetsToLL->Draw("hist,same");
  bdtout_cat4_Data->Draw("e,same");
  //bdtout_cat4_Hgg90->Scale(bdtout_cat4_Data->Integral()/bdtout_cat4_Hgg90->Integral());
  //bdtout_cat4_Hgg90->Draw("hist,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtout_basecat->cd(5);
  gPad->SetGrid();
  ratio_bdtout = (TH1*)bdtout_cat1_Data->Clone();
  ratio_bdtout->Sumw2();
  ratio_bdtout->Divide(bdtout_cat1_DYJetsToLL);
  ratio_bdtout->SetMaximum(1.4);
  ratio_bdtout->SetMinimum(0.6);
  ratio_bdtout->Draw("e");
  line->Draw();

  c_bdtout_basecat->cd(6);
  gPad->SetGrid();
  ratio_bdtoutEB = (TH1*)bdtout_cat2_Data->Clone();
  ratio_bdtoutEB->Sumw2();
  ratio_bdtoutEB->Divide(bdtout_cat2_DYJetsToLL);
  ratio_bdtoutEB->SetMaximum(1.4);
  ratio_bdtoutEB->SetMinimum(0.6);
  ratio_bdtoutEB->Draw("e");
  line->Draw();

  c_bdtout_basecat->cd(7);
  gPad->SetGrid();
  ratio_bdtoutEBEE = (TH1*)bdtout_cat3_Data->Clone();
  ratio_bdtoutEBEE->Sumw2();
  ratio_bdtoutEBEE->Divide(bdtout_cat3_DYJetsToLL);
  ratio_bdtoutEBEE->SetMaximum(1.8);
  ratio_bdtoutEBEE->SetMinimum(0.2);
  ratio_bdtoutEBEE->Draw("e");
  line->Draw();

  c_bdtout_basecat->cd(8);
  gPad->SetGrid();
  ratio_bdtoutEE = (TH1*)bdtout_cat4_Data->Clone();
  ratio_bdtoutEE->Sumw2();
  ratio_bdtoutEE->Divide(bdtout_cat4_DYJetsToLL);
  ratio_bdtoutEE->SetMaximum(1.8);
  ratio_bdtoutEE->SetMinimum(0.2);
  ratio_bdtoutEE->Draw("e");
  line->Draw();

  c_bdtout_basecat->SaveAs("bdtout_basecat.png");



  TCanvas *c_bdtin_1 = new TCanvas("c_bdtin_1","BDT input variables",1600,800);
  c_bdtin_1->Divide(4,2);

  c_bdtin_1->cd(1);
  sf = pho1_phoidMva_cat0_Data->Integral()/pho1_phoidMva_cat0_DYJetsToLL->Integral();
  pho1_phoidMva_cat0_DYJetsToLL->Scale(sf);
  pho1_phoidMva_cat0_Data->GetXaxis()->SetRangeUser(-0.3,0.6);
  pho1_phoidMva_cat0_Data->Draw("e");
  pho1_phoidMva_cat0_DYJetsToLL->Draw("hist,same");
  pho1_phoidMva_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_1->cd(2);
  sf = pho2_phoidMva_cat0_Data->Integral()/pho2_phoidMva_cat0_DYJetsToLL->Integral();
  pho2_phoidMva_cat0_DYJetsToLL->Scale(sf);
  pho2_phoidMva_cat0_Data->GetXaxis()->SetRangeUser(-0.3,0.6);
  pho2_phoidMva_cat0_Data->Draw("e");
  pho2_phoidMva_cat0_DYJetsToLL->Draw("hist,same");
  pho2_phoidMva_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_1->cd(3);
  gPad->SetLogy();
  sf = vtxProb_cat0_Data->Integral()/vtxProb_cat0_DYJetsToLL->Integral();
  vtxProb_cat0_DYJetsToLL->Scale(sf);
  //vtxProb_cat0_Data->GetXaxis()->SetRangeUser(0.2,1.);
  vtxProb_cat0_Data->Draw("e");
  vtxProb_cat0_DYJetsToLL->Draw("hist,same");
  vtxProb_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtin_1->cd(4);
  gPad->SetLogy();
  sf = cosDeltaPhi_cat0_Data->Integral()/cosDeltaPhi_cat0_DYJetsToLL->Integral();
  cout << "sf_cosdphi " << sf << endl;
  cosDeltaPhi_cat0_DYJetsToLL->Scale(sf);
  cosDeltaPhi_cat0_Data->Draw("e");
  cosDeltaPhi_cat0_DYJetsToLL->Draw("hist,same");
  cosDeltaPhi_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_1->cd(5);
  gPad->SetGrid();
  ratio_pho1_phoidMva = (TH1*)pho1_phoidMva_cat0_Data->Clone();
  ratio_pho1_phoidMva->Sumw2();
  ratio_pho1_phoidMva->Divide(pho1_phoidMva_cat0_DYJetsToLL);
  ratio_pho1_phoidMva->SetMaximum(1.8);
  ratio_pho1_phoidMva->SetMinimum(0.2);
  ratio_pho1_phoidMva->Draw("e");
  line1 = (TLine*)line->Clone();
  line1->SetX1(-0.3);
  line1->SetX2(0.6);
  line1->Draw();

  c_bdtin_1->cd(6);
  gPad->SetGrid();
  ratio_pho2_phoidMva = (TH1*)pho2_phoidMva_cat0_Data->Clone();
  ratio_pho2_phoidMva->Sumw2();
  ratio_pho2_phoidMva->Divide(pho2_phoidMva_cat0_DYJetsToLL);
  ratio_pho2_phoidMva->SetMaximum(1.8);
  ratio_pho2_phoidMva->SetMinimum(0.2);
  ratio_pho2_phoidMva->Draw("e");
  line1->Draw();

  c_bdtin_1->cd(7);
  gPad->SetGrid();
  ratio_vtxProb = (TH1*)vtxProb_cat0_Data->Clone();
  ratio_vtxProb->Sumw2();
  ratio_vtxProb->Divide(vtxProb_cat0_DYJetsToLL);
  ratio_vtxProb->SetMaximum(1.8);
  ratio_vtxProb->SetMinimum(0.2);
  ratio_vtxProb->Draw("e");
  line2 = (TLine*)line->Clone();
  //line2->SetX1(0.2);
  line2->Draw();

  c_bdtin_1->cd(8);
  gPad->SetGrid();
  ratio_cosDeltaPhi = (TH1*)cosDeltaPhi_cat0_Data->Clone();
  ratio_cosDeltaPhi->Sumw2();
  ratio_cosDeltaPhi->Divide(cosDeltaPhi_cat0_DYJetsToLL);
  ratio_cosDeltaPhi->SetMaximum(1.8);
  ratio_cosDeltaPhi->SetMinimum(0.2);
  ratio_cosDeltaPhi->Draw("e");
  line->Draw();

  c_bdtin_1->SaveAs("bdtin_1.png");


  TCanvas *c_bdtin_2 = new TCanvas("c_bdtin_2","BDT input variables",1500,900);
  c_bdtin_2->Divide(4,3);
  c_bdtin_2->cd(1);

  sf = sigmaMOverM_cat0_Data->Integral()/sigmaMOverM_cat0_DYJetsToLL->Integral();
  sigmaMOverM_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_cat0_Data->Draw("e");
  sigmaMOverM_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(2);
  sf = sigmaMOverM_wrongVtx_cat0_Data->Integral()/sigmaMOverM_wrongVtx_cat0_DYJetsToLL->Integral();
  sigmaMOverM_wrongVtx_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_wrongVtx_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_bdtin_2->cd(3);
  sf = pho1_ptOverM_cat0_Data->Integral()/pho1_ptOverM_cat0_DYJetsToLL->Integral();
  pho1_ptOverM_cat0_DYJetsToLL->Scale(sf);
  pho1_ptOverM_cat0_Data->GetXaxis()->SetRangeUser(0.2,1.2);
  pho1_ptOverM_cat0_Data->Draw("e");
  pho1_ptOverM_cat0_DYJetsToLL->Draw("hist,same");
  pho1_ptOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(4);
  sf = pho2_ptOverM_cat0_Data->Integral()/pho2_ptOverM_cat0_DYJetsToLL->Integral();
  pho2_ptOverM_cat0_DYJetsToLL->Scale(sf);
  pho2_ptOverM_cat0_Data->GetXaxis()->SetRangeUser(0.2,.9);
  pho2_ptOverM_cat0_Data->Draw("e");
  pho2_ptOverM_cat0_DYJetsToLL->Draw("hist,same");
  pho2_ptOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(5);
  gPad->SetLogy();
  sigmaMOverM_cat0_Data->Draw("e");
  sigmaMOverM_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(6);
  gPad->SetLogy();
  sigmaMOverM_wrongVtx_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtin_2->cd(7);
  gPad->SetLogy();
  pho1_ptOverM_cat0_Data->Draw("e");
  pho1_ptOverM_cat0_DYJetsToLL->Draw("hist,same");
  pho1_ptOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(8);
  gPad->SetLogy();
  pho2_ptOverM_cat0_Data->Draw("e");
  pho2_ptOverM_cat0_DYJetsToLL->Draw("hist,same");
  pho2_ptOverM_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(9);
  gPad->SetGrid();
  ratio_sigmaMOverM = (TH1*)sigmaMOverM_cat0_Data->Clone();
  ratio_sigmaMOverM->Sumw2();
  ratio_sigmaMOverM->Divide(sigmaMOverM_cat0_DYJetsToLL);
  ratio_sigmaMOverM->SetMaximum(1.8);
  ratio_sigmaMOverM->SetMinimum(0.2);
  ratio_sigmaMOverM->Draw("e");
  line3 = (TLine*)line->Clone();
  line3->SetX1(0.);
  line3->SetX2(0.06);
  line3->Draw();

  c_bdtin_2->cd(10);
  gPad->SetGrid();
  ratio_sigmaMOverM_wrongVtx = (TH1*)sigmaMOverM_wrongVtx_cat0_Data->Clone();
  ratio_sigmaMOverM_wrongVtx->Sumw2();
  ratio_sigmaMOverM_wrongVtx->Divide(sigmaMOverM_wrongVtx_cat0_DYJetsToLL);
  ratio_sigmaMOverM_wrongVtx->SetMaximum(1.8);
  ratio_sigmaMOverM_wrongVtx->SetMinimum(0.2);
  ratio_sigmaMOverM_wrongVtx->Draw("e");
  line3->Draw();

  c_bdtin_2->cd(11);
  gPad->SetGrid();
  ratio_pho1_ptOverM = (TH1*)pho1_ptOverM_cat0_Data->Clone();
  ratio_pho1_ptOverM->Sumw2();
  ratio_pho1_ptOverM->Divide(pho1_ptOverM_cat0_DYJetsToLL);
  ratio_pho1_ptOverM->SetMaximum(1.8);
  ratio_pho1_ptOverM->SetMinimum(0.2);
  ratio_pho1_ptOverM->Draw("e");
  line4 = (TLine*)line->Clone();
  line4->SetX1(0.2);
  line4->SetX2(1.2);
  line4->Draw();

  c_bdtin_2->cd(12);
  gPad->SetGrid();
  ratio_pho2_ptOverM = (TH1*)pho2_ptOverM_cat0_Data->Clone();
  ratio_pho2_ptOverM->Sumw2();
  ratio_pho2_ptOverM->Divide(pho2_ptOverM_cat0_DYJetsToLL);
  ratio_pho2_ptOverM->SetMaximum(1.8);
  ratio_pho2_ptOverM->SetMinimum(0.2);
  ratio_pho2_ptOverM->Draw("e");
  line5 = (TLine*)line->Clone();
  line5->SetX1(0.2);
  line5->SetX2(0.9);
  line5->Draw();

  c_bdtin_2->SaveAs("bdtin_2.png");


  TCanvas *c_bdtin_3 = new TCanvas("c_bdtin_3","BDT input variables",1200,800);
  c_bdtin_3->Divide(2,2);
  c_bdtin_3->cd(1);

  c_bdtin_3->cd(1);
  sf = pho1_eta_cat0_Data->Integral()/pho1_eta_cat0_DYJetsToLL->Integral();
  pho1_eta_cat0_DYJetsToLL->Scale(sf);
  pho1_eta_cat0_Data->Draw("e");
  pho1_eta_cat0_DYJetsToLL->Draw("hist,same");
  pho1_eta_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtin_3->cd(2);
  sf = pho2_eta_cat0_Data->Integral()/pho2_eta_cat0_DYJetsToLL->Integral();
  pho2_eta_cat0_DYJetsToLL->Scale(sf);
  pho2_eta_cat0_Data->Draw("e");
  pho2_eta_cat0_DYJetsToLL->Draw("hist,same");
  pho2_eta_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtin_3->cd(3);
  gPad->SetGrid();
  ratio_pho1_eta = (TH1*)pho1_eta_cat0_Data->Clone();
  ratio_pho1_eta->Sumw2();
  ratio_pho1_eta->Divide(pho1_eta_cat0_DYJetsToLL);
  ratio_pho1_eta->SetMaximum(1.8);
  ratio_pho1_eta->SetMinimum(0.2);
  ratio_pho1_eta->Draw("e");
  line6 = (TLine*)line->Clone();
  line6->SetX1(-2.5);
  line6->SetX2(2.5);
  line6->Draw();

  c_bdtin_3->cd(4);
  gPad->SetGrid();
  ratio_pho2_eta = (TH1*)pho2_eta_cat0_Data->Clone();
  ratio_pho2_eta->Sumw2();
  ratio_pho2_eta->Divide(pho2_eta_cat0_DYJetsToLL);
  ratio_pho2_eta->SetMaximum(1.8);
  ratio_pho2_eta->SetMinimum(0.2);
  ratio_pho2_eta->Draw("e");
  line6->Draw();

  c_bdtin_3->SaveAs("bdtin_3.png");


  TCanvas *c_sigmam_ebee = new TCanvas("c_sigmam_ebee","BDT input variables",1500,900);
  c_sigmam_ebee->Divide(4,3);
  c_sigmam_ebee->cd(1);

  sf = sigmaMOverM_EB_cat0_Data->Integral()/sigmaMOverM_EB_cat0_DYJetsToLL->Integral();
  sigmaMOverM_EB_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_EB_cat0_Data->Draw("e");
  sigmaMOverM_EB_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_sigmam_ebee->cd(2);
  sf = sigmaMOverM_wrongVtx_EB_cat0_Data->Integral()/sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->Integral();
  sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_wrongVtx_EB_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_sigmam_ebee->cd(3);
  sf = sigmaMOverM_EE_cat0_Data->Integral()/sigmaMOverM_EE_cat0_DYJetsToLL->Integral();
  sigmaMOverM_EE_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_EE_cat0_Data->Draw("e");
  sigmaMOverM_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_sigmam_ebee->cd(4);
  sf = sigmaMOverM_wrongVtx_EE_cat0_Data->Integral()/sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Integral();
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg2->Draw();

  c_sigmam_ebee->cd(5);
  gPad->SetLogy();
  sigmaMOverM_EB_cat0_Data->Draw("e");
  sigmaMOverM_EB_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_sigmam_ebee->cd(6);
  gPad->SetLogy();
  sigmaMOverM_wrongVtx_EB_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_sigmam_ebee->cd(7);
  gPad->SetLogy();
  sigmaMOverM_EE_cat0_Data->Draw("e");
  sigmaMOverM_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_sigmam_ebee->cd(8);
  gPad->SetLogy();
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_sigmam_ebee->cd(9);
  gPad->SetGrid();
  ratio_sigmaMOverM_EB_cat0 = (TH1*)sigmaMOverM_EB_cat0_Data->Clone();
  ratio_sigmaMOverM_EB_cat0->Sumw2();
  ratio_sigmaMOverM_EB_cat0->Divide(sigmaMOverM_EB_cat0_DYJetsToLL);
  ratio_sigmaMOverM_EB_cat0->SetMaximum(1.8);
  ratio_sigmaMOverM_EB_cat0->SetMinimum(0.2);
  ratio_sigmaMOverM_EB_cat0->Draw("e");
  line3->Draw();

  c_sigmam_ebee->cd(10);
  gPad->SetGrid();
  ratio_sigmaMOverM_wrongVtx_EB_cat0 = (TH1*)sigmaMOverM_wrongVtx_EB_cat0_Data->Clone();
  ratio_sigmaMOverM_wrongVtx_EB_cat0->Sumw2();
  ratio_sigmaMOverM_wrongVtx_EB_cat0->Divide(sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL);
  ratio_sigmaMOverM_wrongVtx_EB_cat0->SetMaximum(1.8);
  ratio_sigmaMOverM_wrongVtx_EB_cat0->SetMinimum(0.2);
  ratio_sigmaMOverM_wrongVtx_EB_cat0->Draw("e");
  line3->Draw();

  c_sigmam_ebee->cd(11);
  gPad->SetGrid();
  ratio_sigmaMOverM_EE_cat0 = (TH1*)sigmaMOverM_EE_cat0_Data->Clone();
  ratio_sigmaMOverM_EE_cat0->Sumw2();
  ratio_sigmaMOverM_EE_cat0->Divide(sigmaMOverM_EE_cat0_DYJetsToLL);
  ratio_sigmaMOverM_EE_cat0->SetMaximum(1.8);
  ratio_sigmaMOverM_EE_cat0->SetMinimum(0.2);
  ratio_sigmaMOverM_EE_cat0->Draw("e");
  line3->Draw();

  c_sigmam_ebee->cd(12);
  gPad->SetGrid();
  ratio_sigmaMOverM_wrongVtx_EE_cat0 = (TH1*)sigmaMOverM_wrongVtx_EE_cat0_Data->Clone();
  ratio_sigmaMOverM_wrongVtx_EE_cat0->Sumw2();
  ratio_sigmaMOverM_wrongVtx_EE_cat0->Divide(sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL);
  ratio_sigmaMOverM_wrongVtx_EE_cat0->SetMaximum(1.8);
  ratio_sigmaMOverM_wrongVtx_EE_cat0->SetMinimum(0.2);
  ratio_sigmaMOverM_wrongVtx_EE_cat0->Draw("e");
  line3->Draw();

  c_sigmam_ebee->SaveAs("sigmam_ebee.png");


}
