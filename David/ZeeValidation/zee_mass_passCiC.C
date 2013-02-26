void zee_mass_passCiC(bool equalArea=true) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString preselNorm_str = "";
  if (!equalArea) preselNorm_str = "_preselNorm";

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  //TFile *file = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/zee_trees/tree_zee_moriond_preapproval.root");
  file->cd();

  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  cout << "sf_presel " << sf_presel << endl;

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(mass_passCiC_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(mass_passCiC_cat0_DYJetsToLL,"DYJetsToLL MC","F");


  mass_passCiC_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_basecat_passCiC_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_basecat_passCiC_cat1_Data->GetXaxis()->SetTitleSize(0.06);
  mass_basecat_passCiC_cat2_Data->GetXaxis()->SetTitleSize(0.06);
  mass_basecat_passCiC_cat3_Data->GetXaxis()->SetTitleSize(0.06);
  mass_passCiC_pt0to20_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_passCiC_pt20to40_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_passCiC_pt40to60_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_passCiC_pt60to100_cat0_Data->GetXaxis()->SetTitleSize(0.06);
  mass_passCiC_pt100up_cat0_Data->GetXaxis()->SetTitleSize(0.06);

  TCanvas *c_mass = new TCanvas("c_mass","di-photon mass",1200,800);
  c_mass->Divide(2,2);

  mass_passCiC_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}");
  mass_passCiC_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_cat0_Data->SetLineColor(1);
  mass_passCiC_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_cat0_DYJetsToLL->SetLineColor(1);
  double sf_all = mass_passCiC_cat0_Data->Integral()/mass_passCiC_cat0_DYJetsToLL->Integral();
  mass_passCiC_cat0_DYJetsToLL->Scale(sf_all);
  cout << sf_all << endl;

  c_mass->cd(1);
  mass_Data = (TH1*)mass_passCiC_cat0_Data->Clone();
  mass_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_Data->Draw("e");
  mass_passCiC_cat0_DYJetsToLL->Draw("hist,same");
  mass_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass->cd(2);
  c_mass_2->SetLogy();
  mass_passCiC_cat0_Data->GetXaxis()->SetRangeUser(70.,180.);
  mass_passCiC_cat0_Data->Draw("e");
  mass_passCiC_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass->cd(3);
  c_mass_3->SetGrid();
  mass_ratio = (TH1F*)mass_passCiC_cat0_Data->Clone();
  mass_ratio->Divide(mass_passCiC_cat0_DYJetsToLL);
  mass_ratio->SetMaximum(1.8);
  mass_ratio->SetMinimum(0.2);
  mass_ratio->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio->GetXaxis()->SetTitle("m_{e^{+}e^{-}}");
  mass_ratio->Draw("e");

  txt3 = new TLatex();
  txt3->SetNDC();
  txt3->SetTextSize(0.04);
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");

  TLine *line_mass = new TLine(60.,1.,120.,1.);
  line_mass->SetLineColor(4);
  line_mass->SetLineWidth(2);
  line_mass->Draw();

  c_mass->SaveAs("mass_passCiC.png");



  TCanvas *c_mass_basecat_passCiC = new TCanvas("c_mass_basecat_passCiC","di-photon mass",1600,600);
  c_mass_basecat_passCiC->Divide(4,2);

  /*
  mass_basecat_passCiC_cat0_Data->Rebin(2);
  mass_basecat_passCiC_cat1_Data->Rebin(2);
  mass_basecat_passCiC_cat2_Data->Rebin(2);
  mass_basecat_passCiC_cat3_Data->Rebin(2);
  mass_basecat_passCiC_cat0_DYJetsToLL->Rebin(2);
  mass_basecat_passCiC_cat1_DYJetsToLL->Rebin(2);
  mass_basecat_passCiC_cat2_DYJetsToLL->Rebin(2);
  mass_basecat_passCiC_cat3_DYJetsToLL->Rebin(2);
  */

  mass_basecat_passCiC_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat0");
  mass_basecat_passCiC_cat0_Data->GetYaxis()->SetTitle("");
  mass_basecat_passCiC_cat0_Data->SetMarkerStyle(20);
  mass_basecat_passCiC_cat0_Data->SetMarkerSize(.4);
  mass_basecat_passCiC_cat0_Data->SetLineColor(1);
  mass_basecat_passCiC_cat0_DYJetsToLL->SetFillColor(38);
  mass_basecat_passCiC_cat0_DYJetsToLL->SetLineColor(1);
  float sf = mass_basecat_passCiC_cat0_Data->Integral()/mass_basecat_passCiC_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_basecat_passCiC_cat0_DYJetsToLL->Scale(sf);

  c_mass_basecat_passCiC->cd(1);
  mass_basecat_passCiC_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_basecat_passCiC_cat0_Data->Draw("e");
  mass_basecat_passCiC_cat0_DYJetsToLL->Draw("hist,same");
  mass_basecat_passCiC_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_basecat_passCiC->cd(5);
  c_mass_basecat_passCiC_5->SetGrid();
  mass_ratio_basecat_cat0 = (TH1F*)mass_basecat_passCiC_cat0_Data->Clone();
  mass_ratio_basecat_cat0->Divide(mass_basecat_passCiC_cat0_DYJetsToLL);
  mass_ratio_basecat_cat0->SetMaximum(1.8);
  mass_ratio_basecat_cat0->SetMinimum(0.2);
  mass_ratio_basecat_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_basecat_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat0");
  mass_ratio_basecat_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_basecat_passCiC_cat1_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat1");
  mass_basecat_passCiC_cat1_Data->GetYaxis()->SetTitle("");
  mass_basecat_passCiC_cat1_Data->SetMarkerStyle(20);
  mass_basecat_passCiC_cat1_Data->SetMarkerSize(.4);
  mass_basecat_passCiC_cat1_Data->SetLineColor(1);
  mass_basecat_passCiC_cat1_DYJetsToLL->SetFillColor(38);
  mass_basecat_passCiC_cat1_DYJetsToLL->SetLineColor(1);
  sf = mass_basecat_passCiC_cat1_Data->Integral()/mass_basecat_passCiC_cat1_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_basecat_passCiC_cat1_DYJetsToLL->Scale(sf);

  c_mass_basecat_passCiC->cd(2);
  mass_basecat_passCiC_cat1_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_basecat_passCiC_cat1_Data->Draw("e");
  mass_basecat_passCiC_cat1_DYJetsToLL->Draw("hist,same");
  mass_basecat_passCiC_cat1_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_basecat_passCiC->cd(6);
  c_mass_basecat_passCiC_6->SetGrid();
  mass_ratio_basecat_cat1 = (TH1F*)mass_basecat_passCiC_cat1_Data->Clone();
  mass_ratio_basecat_cat1->Divide(mass_basecat_passCiC_cat1_DYJetsToLL);
  mass_ratio_basecat_cat1->SetMaximum(1.8);
  mass_ratio_basecat_cat1->SetMinimum(0.2);
  mass_ratio_basecat_cat1->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_basecat_cat1->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat1");
  mass_ratio_basecat_cat1->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_basecat_passCiC_cat2_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat2");
  mass_basecat_passCiC_cat2_Data->GetYaxis()->SetTitle("");
  mass_basecat_passCiC_cat2_Data->SetMarkerStyle(20);
  mass_basecat_passCiC_cat2_Data->SetMarkerSize(.4);
  mass_basecat_passCiC_cat2_Data->SetLineColor(1);
  mass_basecat_passCiC_cat2_DYJetsToLL->SetFillColor(38);
  mass_basecat_passCiC_cat2_DYJetsToLL->SetLineColor(1);
  sf = mass_basecat_passCiC_cat2_Data->Integral()/mass_basecat_passCiC_cat2_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_basecat_passCiC_cat2_DYJetsToLL->Scale(sf);

  c_mass_basecat_passCiC->cd(3);
  mass_basecat_passCiC_cat2_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_basecat_passCiC_cat2_Data->Draw("e");
  mass_basecat_passCiC_cat2_DYJetsToLL->Draw("hist,same");
  mass_basecat_passCiC_cat2_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_basecat_passCiC->cd(7);
  c_mass_basecat_passCiC_7->SetGrid();
  mass_ratio_basecat_cat2 = (TH1F*)mass_basecat_passCiC_cat2_Data->Clone();
  mass_ratio_basecat_cat2->Divide(mass_basecat_passCiC_cat2_DYJetsToLL);
  mass_ratio_basecat_cat2->SetMaximum(1.8);
  mass_ratio_basecat_cat2->SetMinimum(0.2);
  mass_ratio_basecat_cat2->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_basecat_cat2->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat2");
  mass_ratio_basecat_cat2->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_basecat_passCiC_cat3_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat3");
  mass_basecat_passCiC_cat3_Data->GetYaxis()->SetTitle("");
  mass_basecat_passCiC_cat3_Data->SetMarkerStyle(20);
  mass_basecat_passCiC_cat3_Data->SetMarkerSize(.4);
  mass_basecat_passCiC_cat3_Data->SetLineColor(1);
  mass_basecat_passCiC_cat3_DYJetsToLL->SetFillColor(38);
  mass_basecat_passCiC_cat3_DYJetsToLL->SetLineColor(1);
  sf = mass_basecat_passCiC_cat3_Data->Integral()/mass_basecat_passCiC_cat3_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_basecat_passCiC_cat3_DYJetsToLL->Scale(sf);

  c_mass_basecat_passCiC->cd(4);
  mass_basecat_passCiC_cat3_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_basecat_passCiC_cat3_Data->Draw("e");
  mass_basecat_passCiC_cat3_DYJetsToLL->Draw("hist,same");
  mass_basecat_passCiC_cat3_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_basecat_passCiC->cd(8);
  c_mass_basecat_passCiC_8->SetGrid();
  mass_ratio_basecat_cat3 = (TH1F*)mass_basecat_passCiC_cat3_Data->Clone();
  mass_ratio_basecat_cat3->Divide(mass_basecat_passCiC_cat3_DYJetsToLL);
  mass_ratio_basecat_cat3->SetMaximum(1.8);
  mass_ratio_basecat_cat3->SetMinimum(0.2);
  mass_ratio_basecat_cat3->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_basecat_cat3->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, CiC cat3");
  mass_ratio_basecat_cat3->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  c_mass_basecat_passCiC->SaveAs("mass_basecat_passCiC"+preselNorm_str+".png");


  TCanvas *c_mass_passCiC_pt = new TCanvas("c_mass_passCiC_pt","di-photon mass",1600,600);
  c_mass_passCiC_pt->Divide(5,2);

  /*
  mass_passCiC_pt0to20_cat0_Data->Rebin(2);
  mass_passCiC_pt20to40_cat0_Data->Rebin(2);
  mass_passCiC_pt40to60_cat0_Data->Rebin(2);
  mass_passCiC_pt60to100_cat0_Data->Rebin(2);
  mass_passCiC_pt100up_cat0_Data->Rebin(2);
  mass_passCiC_pt0to20_cat0_DYJetsToLL->Rebin(2);
  mass_passCiC_pt20to40_cat0_DYJetsToLL->Rebin(2);
  mass_passCiC_pt40to60_cat0_DYJetsToLL->Rebin(2);
  mass_passCiC_pt60to100_cat0_DYJetsToLL->Rebin(2);
  mass_passCiC_pt100up_cat0_DYJetsToLL->Rebin(2);
  */

  mass_passCiC_pt0to20_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 0<p_{T}(#gamma#gamma)<20 GeV");
  mass_passCiC_pt0to20_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_pt0to20_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_pt0to20_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_pt0to20_cat0_Data->SetLineColor(1);
  mass_passCiC_pt0to20_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_pt0to20_cat0_DYJetsToLL->SetLineColor(1);
  sf = mass_passCiC_pt0to20_cat0_Data->Integral()/mass_passCiC_pt0to20_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_passCiC_pt0to20_cat0_DYJetsToLL->Scale(sf);

  c_mass_passCiC_pt->cd(1);
  mass_passCiC_pt0to20_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_passCiC_pt0to20_cat0_Data->Draw("e");
  mass_passCiC_pt0to20_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_pt0to20_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_passCiC_pt->cd(6);
  c_mass_passCiC_pt_6->SetGrid();
  mass_ratio_pt0to20_cat0 = (TH1F*)mass_passCiC_pt0to20_cat0_Data->Clone();
  mass_ratio_pt0to20_cat0->Divide(mass_passCiC_pt0to20_cat0_DYJetsToLL);
  mass_ratio_pt0to20_cat0->SetMaximum(1.8);
  mass_ratio_pt0to20_cat0->SetMinimum(0.2);
  mass_ratio_pt0to20_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_pt0to20_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 0<p_{T}(#gamma#gamma)<20 GeV");
  mass_ratio_pt0to20_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_passCiC_pt20to40_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 20<p_{T}(#gamma#gamma)<40 GeV");
  mass_passCiC_pt20to40_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_pt20to40_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_pt20to40_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_pt20to40_cat0_Data->SetLineColor(1);
  mass_passCiC_pt20to40_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_pt20to40_cat0_DYJetsToLL->SetLineColor(1);
  sf = mass_passCiC_pt20to40_cat0_Data->Integral()/mass_passCiC_pt20to40_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_passCiC_pt20to40_cat0_DYJetsToLL->Scale(sf);

  c_mass_passCiC_pt->cd(2);
  mass_passCiC_pt20to40_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_passCiC_pt20to40_cat0_Data->Draw("e");
  mass_passCiC_pt20to40_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_pt20to40_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_passCiC_pt->cd(7);
  c_mass_passCiC_pt_7->SetGrid();
  mass_ratio_pt20to40_cat0 = (TH1F*)mass_passCiC_pt20to40_cat0_Data->Clone();
  mass_ratio_pt20to40_cat0->Divide(mass_passCiC_pt20to40_cat0_DYJetsToLL);
  mass_ratio_pt20to40_cat0->SetMaximum(1.8);
  mass_ratio_pt20to40_cat0->SetMinimum(0.2);
  mass_ratio_pt20to40_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_pt20to40_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 20<p_{T}(#gamma#gamma)<40 GeV");
  mass_ratio_pt20to40_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_passCiC_pt40to60_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 40<p_{T}(#gamma#gamma)<60 GeV");
  mass_passCiC_pt40to60_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_pt40to60_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_pt40to60_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_pt40to60_cat0_Data->SetLineColor(1);
  mass_passCiC_pt40to60_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_pt40to60_cat0_DYJetsToLL->SetLineColor(1);
  sf = mass_passCiC_pt40to60_cat0_Data->Integral()/mass_passCiC_pt40to60_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_passCiC_pt40to60_cat0_DYJetsToLL->Scale(sf);

  c_mass_passCiC_pt->cd(3);
  mass_passCiC_pt40to60_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_passCiC_pt40to60_cat0_Data->Draw("e");
  mass_passCiC_pt40to60_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_pt40to60_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_passCiC_pt->cd(8);
  c_mass_passCiC_pt_8->SetGrid();
  mass_ratio_pt40to60_cat0 = (TH1F*)mass_passCiC_pt40to60_cat0_Data->Clone();
  mass_ratio_pt40to60_cat0->Divide(mass_passCiC_pt40to60_cat0_DYJetsToLL);
  mass_ratio_pt40to60_cat0->SetMaximum(1.8);
  mass_ratio_pt40to60_cat0->SetMinimum(0.2);
  mass_ratio_pt40to60_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_pt40to60_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 40<p_{T}(#gamma#gamma)<60 GeV");
  mass_ratio_pt40to60_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_passCiC_pt60to100_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 60<p_{T}(#gamma#gamma)<100 GeV");
  mass_passCiC_pt60to100_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_pt60to100_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_pt60to100_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_pt60to100_cat0_Data->SetLineColor(1);
  mass_passCiC_pt60to100_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_pt60to100_cat0_DYJetsToLL->SetLineColor(1);
  sf = mass_passCiC_pt60to100_cat0_Data->Integral()/mass_passCiC_pt60to100_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_passCiC_pt60to100_cat0_DYJetsToLL->Scale(sf);

  c_mass_passCiC_pt->cd(4);
  mass_passCiC_pt60to100_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_passCiC_pt60to100_cat0_Data->Draw("e");
  mass_passCiC_pt60to100_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_pt60to100_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_passCiC_pt->cd(9);
  c_mass_passCiC_pt_9->SetGrid();
  mass_ratio_pt60to100_cat0 = (TH1F*)mass_passCiC_pt60to100_cat0_Data->Clone();
  mass_ratio_pt60to100_cat0->Divide(mass_passCiC_pt60to100_cat0_DYJetsToLL);
  mass_ratio_pt60to100_cat0->SetMaximum(1.8);
  mass_ratio_pt60to100_cat0->SetMinimum(0.2);
  mass_ratio_pt60to100_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_pt60to100_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, 60<p_{T}(#gamma#gamma)<100 GeV");
  mass_ratio_pt60to100_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  mass_passCiC_pt100up_cat0_Data->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, p_{T}(#gamma#gamma)>100 GeV");
  mass_passCiC_pt100up_cat0_Data->GetYaxis()->SetTitle("");
  mass_passCiC_pt100up_cat0_Data->SetMarkerStyle(20);
  mass_passCiC_pt100up_cat0_Data->SetMarkerSize(.4);
  mass_passCiC_pt100up_cat0_Data->SetLineColor(1);
  mass_passCiC_pt100up_cat0_DYJetsToLL->SetFillColor(38);
  mass_passCiC_pt100up_cat0_DYJetsToLL->SetLineColor(1);
  sf = mass_passCiC_pt100up_cat0_Data->Integral()/mass_passCiC_pt100up_cat0_DYJetsToLL->Integral();
  if (!equalArea) sf = sf_presel;
  mass_passCiC_pt100up_cat0_DYJetsToLL->Scale(sf);

  c_mass_passCiC_pt->cd(5);
  mass_passCiC_pt100up_cat0_Data->GetXaxis()->SetRangeUser(70.,119.5);
  mass_passCiC_pt100up_cat0_Data->Draw("e");
  mass_passCiC_pt100up_cat0_DYJetsToLL->Draw("hist,same");
  mass_passCiC_pt100up_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_mass_passCiC_pt->cd(10);
  c_mass_passCiC_pt_10->SetGrid();
  mass_ratio_pt100up_cat0 = (TH1F*)mass_passCiC_pt100up_cat0_Data->Clone();
  mass_ratio_pt100up_cat0->Divide(mass_passCiC_pt100up_cat0_DYJetsToLL);
  mass_ratio_pt100up_cat0->SetMaximum(1.8);
  mass_ratio_pt100up_cat0->SetMinimum(0.2);
  mass_ratio_pt100up_cat0->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio_pt100up_cat0->GetXaxis()->SetTitle("m_{e^{+}e^{-}}, p_{T}(#gamma#gamma)>100 GeV");
  mass_ratio_pt100up_cat0->Draw("e");
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
  line_mass->Draw();

  c_mass_passCiC_pt->SaveAs("mass_passCiC_pt"+preselNorm_str+".png");

}
