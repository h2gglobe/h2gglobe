void zee_mass(bool equalArea=true, bool includeDijetCats=false) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString preselNorm_str = "";
  if (!equalArea) preselNorm_str = "_preselNorm";

  int nMvaCats=4;
  if (includeDijetCats) nMvaCats=6;

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();
  cout << "sf_presel " << sf_presel << endl;

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  txt1 = new TLatex();
  txt1->SetNDC();
  txt1->SetTextSize(0.06);
  txt1->SetTextAlign(12);

  TH1* mass_Data[7];
  TH1* mass_MC[7];
  TH1* mass_ratio[7];
  TH1* mass_basecat_Data[4];
  TH1* mass_basecat_MC[4];
  TH1* mass_pt_Data[5];
  TH1* mass_pt_MC[5];
  TH1* mass_nvtx_Data[3];
  TH1* mass_nvtx_MC[3];
  TH1* mass_passMVA_nvtx_Data[3];
  TH1* mass_passMVA_nvtx_MC[3];

  for (int i=0; i<7; i++) {
    TString iStr;
    iStr+=i;
    mass_Data[i] = (TH1*)(file->Get("mass_cat"+iStr+"_Data"))->Clone();
    mass_MC[i] = (TH1*)(file->Get("mass_cat"+iStr+"_DYJetsToLL"))->Clone();
    mass_Data[i]->GetXaxis()->SetTitleSize(0.06);
  }

  for (int i=0; i<4; i++) {
    TString iStr;
    iStr+=i;
    mass_basecat_Data[i] = (TH1*)(file->Get("mass_basecat_cat"+iStr+"_Data"))->Clone();
    mass_basecat_MC[i] = (TH1*)(file->Get("mass_basecat_cat"+iStr+"_DYJetsToLL"))->Clone();
    mass_basecat_Data[i]->GetXaxis()->SetTitleSize(0.06);
  }

  mass_pt_Data[0] = (TH1*)(file->Get("mass_pt0to20_cat0_Data"))->Clone();
  mass_pt_Data[1] = (TH1*)(file->Get("mass_pt20to40_cat0_Data"))->Clone();
  mass_pt_Data[2] = (TH1*)(file->Get("mass_pt40to60_cat0_Data"))->Clone();
  mass_pt_Data[3] = (TH1*)(file->Get("mass_pt60to100_cat0_Data"))->Clone();
  mass_pt_Data[4] = (TH1*)(file->Get("mass_pt100up_cat0_Data"))->Clone();
  mass_pt_MC[0] = (TH1*)(file->Get("mass_pt0to20_cat0_DYJetsToLL"))->Clone();
  mass_pt_MC[1] = (TH1*)(file->Get("mass_pt20to40_cat0_DYJetsToLL"))->Clone();
  mass_pt_MC[2] = (TH1*)(file->Get("mass_pt40to60_cat0_DYJetsToLL"))->Clone();
  mass_pt_MC[3] = (TH1*)(file->Get("mass_pt60to100_cat0_DYJetsToLL"))->Clone();
  mass_pt_MC[4] = (TH1*)(file->Get("mass_pt100up_cat0_DYJetsToLL"))->Clone();
  for (int i=0; i<5; i++) mass_pt_Data[i]->GetXaxis()->SetTitleSize(0.06);

  TString iStr="0";
  mass_nvtx_Data[0] = (TH1*)(file->Get("mass_nvtx0to13_cat"+iStr+"_Data"))->Clone();
  mass_nvtx_MC[0] = (TH1*)(file->Get("mass_nvtx0to13_cat"+iStr+"_DYJetsToLL"))->Clone();
  mass_nvtx_Data[1] = (TH1*)(file->Get("mass_nvtx14to18_cat"+iStr+"_Data"))->Clone();
  mass_nvtx_MC[1] = (TH1*)(file->Get("mass_nvtx14to18_cat"+iStr+"_DYJetsToLL"))->Clone();
  mass_nvtx_Data[2] = (TH1*)(file->Get("mass_nvtx19up_cat"+iStr+"_Data"))->Clone();
  mass_nvtx_MC[2] = (TH1*)(file->Get("mass_nvtx19up_cat"+iStr+"_DYJetsToLL"))->Clone();
  mass_passMVA_nvtx_Data[0] = (TH1*)(file->Get("mass_passMVA_nvtx0to13_cat"+iStr+"_Data"))->Clone();
  mass_passMVA_nvtx_MC[0] = (TH1*)(file->Get("mass_passMVA_nvtx0to13_cat"+iStr+"_DYJetsToLL"))->Clone();
  mass_passMVA_nvtx_Data[1] = (TH1*)(file->Get("mass_passMVA_nvtx14to18_cat"+iStr+"_Data"))->Clone();
  mass_passMVA_nvtx_MC[1] = (TH1*)(file->Get("mass_passMVA_nvtx14to18_cat"+iStr+"_DYJetsToLL"))->Clone();
  mass_passMVA_nvtx_Data[2] = (TH1*)(file->Get("mass_passMVA_nvtx19up_cat"+iStr+"_Data"))->Clone();
  mass_passMVA_nvtx_MC[2] = (TH1*)(file->Get("mass_passMVA_nvtx19up_cat"+iStr+"_DYJetsToLL"))->Clone();

  for (int i=0; i<3; i++) {
    mass_nvtx_Data[i]->GetXaxis()->SetTitleSize(0.06);
    mass_passMVA_nvtx_Data[i]->GetXaxis()->SetTitleSize(0.06);
  }

  TLegend *leg;
  leg = new TLegend(.55,.45,.82,.67);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.045);
  //leg->AddEntry(mass_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(mass_Data[0],"Data");
  leg->AddEntry(mass_MC[0],"Z#rightarrow e^{+}e^{-} MC","F");

  //------------------------------------------------------------------------------

  TCanvas *c_mass = new TCanvas("c_mass","di-photon mass",1200,800);
  c_mass->Divide(2,2);

  mass_Data[0]->GetXaxis()->SetTitle("m_{e^{+}e^{-}}");
  mass_Data[0]->GetYaxis()->SetTitle("");
  mass_Data[0]->SetMarkerStyle(20);
  mass_Data[0]->SetMarkerSize(.4);
  mass_Data[0]->SetLineColor(1);
  mass_MC[0]->SetFillColor(38);
  mass_MC[0]->SetLineColor(1);
  double sf_all = mass_Data[0]->Integral()/mass_MC[0]->Integral();
  mass_MC[0]->Scale(sf_all);
  cout << sf_all << endl;

  c_mass->cd(1);
  mass_Data_clone = (TH1*)mass_Data[0]->Clone();
  mass_Data_clone->GetXaxis()->SetRangeUser(70.,119.5);
  mass_Data_clone->Draw("e");
  mass_MC[0]->Draw("hist,same");
  mass_Data_clone->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();
  txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

  c_mass->cd(2);
  c_mass_2->SetLogy();
  mass_Data[0]->GetXaxis()->SetRangeUser(70.,180.);
  mass_Data[0]->Draw("e");
  mass_MC[0]->Draw("hist,same");
  mass_Data[0]->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();
  txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

  c_mass->cd(3);
  c_mass_3->SetGrid();
  mass_ratio[0] = (TH1F*)mass_Data[0]->Clone();
  mass_ratio[0]->Divide(mass_MC[0]);
  mass_ratio[0]->SetMaximum(1.8);
  mass_ratio[0]->SetMinimum(0.2);
  mass_ratio[0]->GetXaxis()->SetRangeUser(70.,119.5);
  mass_ratio[0]->GetXaxis()->SetTitle("m_{e^{+}e^{-}}");
  mass_ratio[0]->Draw("e");

  txt3 = new TLatex();
  txt3->SetNDC();
  txt3->SetTextSize(0.06);
  txt3->DrawLatex(0.35,0.25, "Data/MC ratio");

  TLine *line_mass = new TLine(70.,1.,120.,1.);
  line_mass->SetLineColor(4);
  line_mass->SetLineWidth(2);
  line_mass->Draw();

  c_mass->SaveAs("mass.png");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_cat = new TCanvas("c_mass_cat","di-photon mass",1600,600);
  c_mass_cat->Divide(nMvaCats,2);

  if (includeDijetCats) {
    c_mass_cat_7->SetGrid();
    c_mass_cat_8->SetGrid();
    c_mass_cat_9->SetGrid();
    c_mass_cat_10->SetGrid();
    c_mass_cat_11->SetGrid();
    c_mass_cat_12->SetGrid();
  } else {
    c_mass_cat_5->SetGrid();
    c_mass_cat_6->SetGrid();
    c_mass_cat_7->SetGrid();
    c_mass_cat_8->SetGrid();
  }

  TString title[7];
  title[1] = "m_{e^{+}e^{-}}, cat0";
  title[2] = "m_{e^{+}e^{-}}, cat1";
  title[3] = "m_{e^{+}e^{-}}, cat2";
  title[4] = "m_{e^{+}e^{-}}, cat3";
  title[5] = "m_{e^{+}e^{-}}, dijet cat1";
  title[6] = "m_{e^{+}e^{-}}, dijet cat2";

  for (int i=1; i<nMvaCats+1; i++) {

    mass_Data[i]->GetXaxis()->SetTitle(title[i]);
    mass_Data[i]->GetYaxis()->SetTitle("");
    mass_Data[i]->SetMarkerStyle(20);
    mass_Data[i]->SetMarkerSize(.4);
    mass_Data[i]->SetLineColor(1);
    mass_MC[i]->SetFillColor(38);
    mass_MC[i]->SetLineColor(1);
    float sf = mass_Data[i]->Integral()/mass_MC[i]->Integral();
    if (!equalArea) sf = sf_presel;
    mass_MC[i]->Scale(sf);

    c_mass_cat->cd(i);
    mass_Data[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_Data[i]->Draw("e");
    mass_MC[i]->Draw("hist,same");
    mass_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_cat->cd(nMvaCats+i);
    mass_ratio[i] = (TH1F*)mass_Data[i]->Clone();
    mass_ratio[i]->Divide(mass_MC[i]);
    mass_ratio[i]->SetMaximum(1.8);
    mass_ratio[i]->SetMinimum(0.2);
    mass_ratio[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_ratio[i]->GetXaxis()->SetTitle(title[i]);
    mass_ratio[i]->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_cat->SaveAs("mass_cat"+preselNorm_str+".png");
  c_mass_cat->SaveAs("mass_cat"+preselNorm_str+".pdf");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_cat_2x2 = new TCanvas("c_mass_cat_2x2","di-photon mass",1000,1400);
  c_mass_cat_2x2->Divide(2,4);

  c_mass_cat_2x2_3->SetGrid();
  c_mass_cat_2x2_4->SetGrid();
  c_mass_cat_2x2_7->SetGrid();
  c_mass_cat_2x2_8->SetGrid();

  for (int i=1; i<5; i++) {

    c_mass_cat_2x2->cd(i<3 ? i : i+2);
    mass_Data[i]->Draw("e");
    mass_MC[i]->Draw("hist,same");
    mass_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_cat_2x2->cd(i<3 ? i+2 : i+4);
    mass_ratio[i]->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_cat_2x2->SaveAs("mass_cat"+preselNorm_str+"_2x2.png");
  c_mass_cat_2x2->SaveAs("mass_cat"+preselNorm_str+"_2x2.pdf");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_basecat = new TCanvas("c_mass_basecat","di-photon mass",1600,600);
  c_mass_basecat->Divide(4,2);

  c_mass_basecat_5->SetGrid();
  c_mass_basecat_6->SetGrid();
  c_mass_basecat_7->SetGrid();
  c_mass_basecat_8->SetGrid();

  title[0] = "m_{e^{+}e^{-}}, CiC cat0";
  title[1] = "m_{e^{+}e^{-}}, CiC cat1";
  title[2] = "m_{e^{+}e^{-}}, CiC cat2";
  title[3] = "m_{e^{+}e^{-}}, CiC cat3";

  for (int i=0; i<4; i++) {

    mass_basecat_Data[i]->GetXaxis()->SetTitle(title[i]);
    mass_basecat_Data[i]->GetYaxis()->SetTitle("");
    mass_basecat_Data[i]->SetMarkerStyle(20);
    mass_basecat_Data[i]->SetMarkerSize(.4);
    mass_basecat_Data[i]->SetLineColor(1);
    mass_basecat_MC[i]->SetFillColor(38);
    mass_basecat_MC[i]->SetLineColor(1);
    sf = mass_basecat_Data[i]->Integral()/mass_basecat_MC[i]->Integral();
    if (!equalArea) sf = sf_presel;
    mass_basecat_MC[i]->Scale(sf);

    c_mass_basecat->cd(i+1);
    mass_basecat_Data[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_basecat_Data[i]->Draw("e");
    mass_basecat_MC[i]->Draw("hist,same");
    mass_basecat_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_basecat->cd(i+5);
    mass_ratio_basecat = (TH1F*)mass_basecat_Data[i]->Clone();
    mass_ratio_basecat->Divide(mass_basecat_MC[i]);
    mass_ratio_basecat->SetMaximum(1.8);
    mass_ratio_basecat->SetMinimum(0.2);
    mass_ratio_basecat->GetXaxis()->SetRangeUser(70.,119.5);
    mass_ratio_basecat->GetXaxis()->SetTitle(title[i]);
    mass_ratio_basecat->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_basecat->SaveAs("mass_basecat"+preselNorm_str+".png");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_pt = new TCanvas("c_mass_pt","di-photon mass",1600,600);
  c_mass_pt->Divide(5,2);

  c_mass_pt_6->SetGrid();
  c_mass_pt_7->SetGrid();
  c_mass_pt_8->SetGrid();
  c_mass_pt_9->SetGrid();
  c_mass_pt_10->SetGrid();

  title[0] = "m_{e^{+}e^{-}}, 0<p_{T}(#gamma#gamma)<20 GeV";
  title[1] = "m_{e^{+}e^{-}}, 20<p_{T}(#gamma#gamma)<40 GeV";
  title[2] = "m_{e^{+}e^{-}}, 40<p_{T}(#gamma#gamma)<60 GeV";
  title[3] = "m_{e^{+}e^{-}}, 60<p_{T}(#gamma#gamma)<100 GeV";
  title[4] = "m_{e^{+}e^{-}}, p_{T}(#gamma#gamma)>100 GeV";

  for (int i=0; i<5; i++) {

    mass_pt_Data[i]->GetYaxis()->SetTitle(title[i]);
    mass_pt_Data[i]->SetMarkerStyle(20);
    mass_pt_Data[i]->SetMarkerSize(.4);
    mass_pt_Data[i]->SetLineColor(1);
    mass_pt_MC[i]->SetFillColor(38);
    mass_pt_MC[i]->SetLineColor(1);
    sf = mass_pt_Data[i]->Integral()/mass_pt_MC[i]->Integral();
    if (!equalArea) sf = sf_presel;
    mass_pt_MC[i]->Scale(sf);

    c_mass_pt->cd(i+1);
    mass_pt_Data[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_pt_Data[i]->Draw("e");
    mass_pt_MC[i]->Draw("hist,same");
    mass_pt_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_pt->cd(i+6);
    mass_ratio_pt = (TH1F*)mass_pt_Data[i]->Clone();
    mass_ratio_pt->Divide(mass_pt_MC[i]);
    mass_ratio_pt->SetMaximum(1.8);
    mass_ratio_pt->SetMinimum(0.2);
    mass_ratio_pt->GetXaxis()->SetRangeUser(70.,119.5);
    mass_ratio_pt->GetXaxis()->SetTitle(title[i]);
    mass_ratio_pt->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_pt->SaveAs("mass_pt"+preselNorm_str+".png");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_nvtx = new TCanvas("c_mass_nvtx","di-photon mass",1600,600);
  c_mass_nvtx->Divide(3,2);

  c_mass_nvtx_4->SetGrid();
  c_mass_nvtx_5->SetGrid();
  c_mass_nvtx_6->SetGrid();

  title[0] = "m_{e^{+}e^{-}}, nvtx<=13";
  title[1] = "m_{e^{+}e^{-}}, 14<=nvtx<=18";
  title[2] = "m_{e^{+}e^{-}}, nvtx>=19";

  for (int i=0; i<3; i++) {

    mass_nvtx_Data[i]->GetXaxis()->SetTitle(title[i]);
    mass_nvtx_Data[i]->GetYaxis()->SetTitle("");
    mass_nvtx_Data[i]->SetMarkerStyle(20);
    mass_nvtx_Data[i]->SetMarkerSize(.4);
    mass_nvtx_Data[i]->SetLineColor(1);
    mass_nvtx_MC[i]->SetFillColor(38);
    mass_nvtx_MC[i]->SetLineColor(1);
    sf = mass_nvtx_Data[i]->Integral()/mass_nvtx_MC[i]->Integral();
    if (!equalArea) sf = sf_presel;
    mass_nvtx_MC[i]->Scale(sf);

    c_mass_nvtx->cd(i+1);
    mass_nvtx_Data[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_nvtx_Data[i]->Draw("e");
    mass_nvtx_MC[i]->Draw("hist,same");
    mass_nvtx_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_nvtx->cd(i+4);
    mass_ratio_nvtx = (TH1F*)mass_nvtx_Data[i]->Clone();
    mass_ratio_nvtx->Divide(mass_nvtx_MC[i]);
    mass_ratio_nvtx->SetMaximum(1.8);
    mass_ratio_nvtx->SetMinimum(0.2);
    mass_ratio_nvtx->GetXaxis()->SetRangeUser(70.,119.5);
    mass_ratio_nvtx->GetXaxis()->SetTitle(title[i]);
    mass_ratio_nvtx->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_nvtx->SaveAs("mass_nvtx"+preselNorm_str+".png");

  //------------------------------------------------------------------------------

  TCanvas *c_mass_passMVA_nvtx = new TCanvas("c_mass_passMVA_nvtx","di-photon mass",1600,600);
  c_mass_passMVA_nvtx->Divide(3,2);

  c_mass_passMVA_nvtx_4->SetGrid();
  c_mass_passMVA_nvtx_5->SetGrid();
  c_mass_passMVA_nvtx_6->SetGrid();

  for (int i=0; i<3; i++) {

    mass_passMVA_nvtx_Data[i]->GetXaxis()->SetTitle(title[i]);
    mass_passMVA_nvtx_Data[i]->GetYaxis()->SetTitle("");
    mass_passMVA_nvtx_Data[i]->SetMarkerStyle(20);
    mass_passMVA_nvtx_Data[i]->SetMarkerSize(.4);
    mass_passMVA_nvtx_Data[i]->SetLineColor(1);
    mass_passMVA_nvtx_MC[i]->SetFillColor(38);
    mass_passMVA_nvtx_MC[i]->SetLineColor(1);
    sf = mass_passMVA_nvtx_Data[i]->Integral()/mass_passMVA_nvtx_MC[i]->Integral();
    if (!equalArea) sf = sf_presel;
    mass_passMVA_nvtx_MC[i]->Scale(sf);

    c_mass_passMVA_nvtx->cd(i+1);
    mass_passMVA_nvtx_Data[i]->GetXaxis()->SetRangeUser(70.,119.5);
    mass_passMVA_nvtx_Data[i]->Draw("e");
    mass_passMVA_nvtx_MC[i]->Draw("hist,same");
    mass_passMVA_nvtx_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    leg->Draw();
    txt1->DrawLatex(0.55,0.8,"#scale[0.8]{#splitline{CMS Preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");

    c_mass_passMVA_nvtx->cd(i+4);
    mass_ratio_nvtx = (TH1F*)mass_passMVA_nvtx_Data[i]->Clone();
    mass_ratio_nvtx->Divide(mass_passMVA_nvtx_MC[i]);
    mass_ratio_nvtx->SetMaximum(1.8);
    mass_ratio_nvtx->SetMinimum(0.2);
    mass_ratio_nvtx->GetXaxis()->SetRangeUser(70.,119.5);
    mass_ratio_nvtx->GetXaxis()->SetTitle(title[i]);
    mass_ratio_nvtx->Draw("e");
    txt3->DrawLatex(0.35,0.25, "Data/MC ratio");
    line_mass->Draw();

  }

  c_mass_passMVA_nvtx->SaveAs("mass_passMVA_nvtx"+preselNorm_str+".png");
  c_mass_passMVA_nvtx->SaveAs("mass_passMVA_nvtx"+preselNorm_str+".pdf");

}
