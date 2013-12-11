void diphomva_inputs(bool passMVAcut=false, bool equalArea=true) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString preselNorm_str = "";
  if (!equalArea) preselNorm_str = "_preselNorm";

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
  bdtoutEE_cat0_DYJetsToLL->Rebin(2);
  bdtout_cat0_Data->Rebin(2);
  bdtout_cat1_Data->Rebin(2);
  bdtout_cat2_Data->Rebin(2);
  bdtout_cat3_Data->Rebin(2);
  bdtout_cat4_Data->Rebin(2);
  bdtoutEB_cat0_Data->Rebin(2);
  bdtoutEBEE_cat0_Data->Rebin(2);
  bdtoutEE_cat0_Data->Rebin(2);

  bdtout_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output");
  bdtout_cat1_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat0");
  bdtout_cat2_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat1");
  bdtout_cat3_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat2");
  bdtout_cat4_Data->GetXaxis()->SetTitle("di-photon BDT output, CiC cat3");
  bdtoutEB_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EB)");
  bdtoutEBEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (one EB, one EE)");
  bdtoutEE_cat0_Data->GetXaxis()->SetTitle("di-photon BDT output (both EE)");
  pho1_phoidMva_cat0_Data->GetXaxis()->SetTitle("lead photon ID MVA output");
  pho2_phoidMva_cat0_Data->GetXaxis()->SetTitle("sublead photon ID MVA output");
  sigmaMOverM_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex)");
  sigmaMOverM_wrongVtx_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex)");
  sigmaMOverM_EB_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex, both EB)");
  sigmaMOverM_wrongVtx_EB_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex, both EB)");
  sigmaMOverM_EE_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (right vertex, both EE)");
  sigmaMOverM_wrongVtx_EE_cat0_Data->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} (wrong vertex, both EE)");
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

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(vtxProb_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(vtxProb_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(vtxProb_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  cout << "sf_presel " << sf_presel << endl;

  //------------------------------------------------------------------------------

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
  bdtout_cat0_Data->Draw("e");
  leg2->Draw();
  bdtout_cat0_DYJetsToLL->Draw("hist,same");
  bdtout_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout->cd(2);
  if (passMVAcut) {
    sf = bdtoutEB_cat0_Data->Integral(48,100)/bdtoutEB_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEB_cat0_Data->Integral()/bdtoutEB_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEB_cat0_DYJetsToLL->Scale(sf);
  bdtoutEB_cat0_Data->Draw("e");
  leg2->Draw();
  bdtoutEB_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout->cd(3);
  if (passMVAcut) {
    sf = bdtoutEBEE_cat0_Data->Integral(48,100)/bdtoutEBEE_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEBEE_cat0_Data->Integral()/bdtoutEBEE_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEBEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEBEE_cat0_Data->Draw("e");
  leg->Draw();
  bdtoutEBEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEBEE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout->cd(4);
  if (passMVAcut) {
    sf = bdtoutEE_cat0_Data->Integral(48,100)/bdtoutEE_cat0_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtoutEE_cat0_Data->Integral()/bdtoutEE_cat0_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtoutEE_cat0_DYJetsToLL->Scale(sf);
  bdtoutEE_cat0_Data->Draw("e");
  leg->Draw();
  bdtoutEE_cat0_DYJetsToLL->Draw("hist,same");
  bdtoutEE_cat0_Data->Draw("e,same");
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
  plotRatio(bdtout_cat0_Data,bdtout_cat0_DYJetsToLL);
  line->Draw();

  c_bdtout->cd(6);
  plotRatio(bdtoutEB_cat0_Data,bdtoutEB_cat0_DYJetsToLL);
  line->Draw();

  c_bdtout->cd(7);
  plotRatio(bdtoutEBEE_cat0_Data,bdtoutEBEE_cat0_DYJetsToLL);
  line->Draw();

  c_bdtout->cd(8);
  plotRatio(bdtoutEE_cat0_Data,bdtoutEE_cat0_DYJetsToLL);
  line->Draw();

  c_bdtout->SaveAs("bdtout"+preselNorm_str+".png");

  //------------------------------------------------------------------------------

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
  bdtout_cat1_Data->Draw("e");
  leg2->Draw();
  bdtout_cat1_DYJetsToLL->Draw("hist,same");
  bdtout_cat1_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(2);
  if (passMVAcut) {
    sf = bdtout_cat2_Data->Integral(48,100)/bdtout_cat2_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat2_Data->Integral()/bdtout_cat2_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat2_DYJetsToLL->Scale(sf);
  bdtout_cat2_Data->Draw("e");
  leg2->Draw();
  bdtout_cat2_DYJetsToLL->Draw("hist,same");
  bdtout_cat2_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(3);
  if (passMVAcut) {
    sf = bdtout_cat3_Data->Integral(48,100)/bdtout_cat3_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat3_Data->Integral()/bdtout_cat3_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat3_DYJetsToLL->Scale(sf);
  bdtout_cat3_Data->Draw("e");
  leg->Draw();
  bdtout_cat3_DYJetsToLL->Draw("hist,same");
  bdtout_cat3_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(4);
  if (passMVAcut) {
    sf = bdtout_cat4_Data->Integral(48,100)/bdtout_cat4_DYJetsToLL->Integral(48,100);
  } else {
    sf = bdtout_cat4_Data->Integral()/bdtout_cat4_DYJetsToLL->Integral();
    if (!equalArea) sf = sf_presel;
  }
  bdtout_cat4_DYJetsToLL->Scale(sf);
  bdtout_cat4_Data->Draw("e");
  leg->Draw();
  bdtout_cat4_DYJetsToLL->Draw("hist,same");
  bdtout_cat4_Data->Draw("e,same");
  gPad->RedrawAxis();

  c_bdtout_basecat->cd(5);
  plotRatio(bdtout_cat1_Data,bdtout_cat1_DYJetsToLL);
  line->Draw();

  c_bdtout_basecat->cd(6);
  plotRatio(bdtout_cat2_Data,bdtout_cat2_DYJetsToLL);
  line->Draw();

  c_bdtout_basecat->cd(7);
  plotRatio(bdtout_cat3_Data,bdtout_cat3_DYJetsToLL);
  line->Draw();

  c_bdtout_basecat->cd(8);
  plotRatio(bdtout_cat4_Data,bdtout_cat4_DYJetsToLL);
  line->Draw();

  c_bdtout_basecat->SaveAs("bdtout_basecat"+preselNorm_str+".png");

  //------------------------------------------------------------------------------

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
  cosDeltaPhi_cat0_DYJetsToLL->Scale(sf);
  cosDeltaPhi_cat0_Data->Draw("e");
  cosDeltaPhi_cat0_DYJetsToLL->Draw("hist,same");
  cosDeltaPhi_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_1->cd(5);
  plotRatio(pho1_phoidMva_cat0_Data,pho1_phoidMva_cat0_DYJetsToLL);
  line1 = (TLine*)line->Clone();
  line1->SetX1(-0.3);
  line1->SetX2(0.6);
  line1->Draw();

  c_bdtin_1->cd(6);
  plotRatio(pho2_phoidMva_cat0_Data,pho2_phoidMva_cat0_DYJetsToLL);
  line1->Draw();

  c_bdtin_1->cd(7);
  plotRatio(vtxProb_cat0_Data,vtxProb_cat0_DYJetsToLL);
  line->Draw();

  c_bdtin_1->cd(8);
  plotRatio(cosDeltaPhi_cat0_Data,cosDeltaPhi_cat0_DYJetsToLL);
  line->Draw();

  c_bdtin_1->SaveAs("bdtin_1.png");

  //------------------------------------------------------------------------------

  TCanvas *c_bdtin_2 = new TCanvas("c_bdtin_2","BDT input variables",1500,900);
  c_bdtin_2->Divide(4,3);
  c_bdtin_2->cd(1);

  sf = sigmaMOverM_EE_cat0_Data->Integral()/sigmaMOverM_EE_cat0_DYJetsToLL->Integral();
  sigmaMOverM_EE_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_EE_cat0_Data->Draw("e");
  sigmaMOverM_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(2);
  sf = sigmaMOverM_wrongVtx_EE_cat0_Data->Integral()/sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Integral();
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Scale(sf);
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e,same");
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
  sigmaMOverM_EE_cat0_Data->Draw("e");
  sigmaMOverM_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_EE_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_bdtin_2->cd(6);
  gPad->SetLogy();
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e");
  sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL->Draw("hist,same");
  sigmaMOverM_wrongVtx_EE_cat0_Data->Draw("e,same");
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
  plotRatio(sigmaMOverM_EE_cat0_Data,sigmaMOverM_EE_cat0_DYJetsToLL);
  line3 = (TLine*)line->Clone();
  line3->SetX1(0.);
  line3->SetX2(0.06);
  line3->Draw();

  c_bdtin_2->cd(10);
  plotRatio(sigmaMOverM_wrongVtx_EE_cat0_Data,sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL);
  line3->Draw();

  c_bdtin_2->cd(11);
  plotRatio(pho1_ptOverM_cat0_Data,pho1_ptOverM_cat0_DYJetsToLL);
  line4 = (TLine*)line->Clone();
  line4->SetX1(0.2);
  line4->SetX2(1.2);
  line4->Draw();

  c_bdtin_2->cd(12);
  plotRatio(pho2_ptOverM_cat0_Data,pho2_ptOverM_cat0_DYJetsToLL);
  line5 = (TLine*)line->Clone();
  line5->SetX1(0.2);
  line5->SetX2(0.9);
  line5->Draw();

  c_bdtin_2->SaveAs("bdtin_2.png");

  //------------------------------------------------------------------------------

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
  plotRatio(pho1_eta_cat0_Data,pho1_eta_cat0_DYJetsToLL);
  line6 = (TLine*)line->Clone();
  line6->SetX1(-2.5);
  line6->SetX2(2.5);
  line6->Draw();

  c_bdtin_3->cd(4);
  plotRatio(pho2_eta_cat0_Data,pho2_eta_cat0_DYJetsToLL);
  line6->Draw();

  c_bdtin_3->SaveAs("bdtin_3.png");

  //------------------------------------------------------------------------------

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
  leg->Draw();

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
  plotRatio(sigmaMOverM_EB_cat0_Data,
	    sigmaMOverM_EB_cat0_DYJetsToLL);
  line3->Draw();

  c_sigmam_ebee->cd(10);
  plotRatio(sigmaMOverM_wrongVtx_EB_cat0_Data,
	    sigmaMOverM_wrongVtx_EB_cat0_DYJetsToLL);
  line3->Draw();

  c_sigmam_ebee->cd(11);
  plotRatio(sigmaMOverM_EE_cat0_Data,
	    sigmaMOverM_EE_cat0_DYJetsToLL);
  line3->Draw();

  c_sigmam_ebee->cd(12);
  plotRatio(sigmaMOverM_wrongVtx_EE_cat0_Data,
	    sigmaMOverM_wrongVtx_EE_cat0_DYJetsToLL);
  line3->Draw();

  c_sigmam_ebee->SaveAs("sigmam_ebee.png");


}

void plotRatio(TH1* hist_Data, TH1* hist_MC)
{
  gPad->SetGrid();
  ratio = (TH1*)hist_Data->Clone();
  ratio->Divide(hist_MC);
  ratio->SetMaximum(1.8);
  ratio->SetMinimum(0.2);
  ratio->Draw("e");
}
