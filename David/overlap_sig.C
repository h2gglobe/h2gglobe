void overlap_sig() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetPalette(1);
  
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetPaintTextFormat("5.2f");

  TFile *file = TFile::Open("histograms_CMS-HGG_categories_sig.root");
  file->cd();

  TCanvas *c_overlap = new TCanvas("c_overlap","passCiC vs passMVA",750,500);

  passCiC_vs_passMVA = (TH1*)passCiC_vs_passMVA_cat0_ggh_m125_8TeV->Clone();
  passCiC_vs_passMVA->Add(passCiC_vs_passMVA_cat0_vbf_m125_8TeV);
  passCiC_vs_passMVA->Add(passCiC_vs_passMVA_cat0_wzh_m125_8TeV);
  passCiC_vs_passMVA->Add(passCiC_vs_passMVA_cat0_tth_m125_8TeV);
  passCiC_vs_passMVA->SetBinContent(1,1,0.);
  passCiC_vs_passMVA->SetMarkerSize(4.);
  passCiC_vs_passMVA->GetXaxis()->SetNdivisions(2);
  passCiC_vs_passMVA->GetYaxis()->SetNdivisions(2);
  passCiC_vs_passMVA->GetXaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA->GetYaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA->GetXaxis()->SetTitle("Selected in MVA analysis");
  passCiC_vs_passMVA->GetYaxis()->SetTitle("Selected in cut based analysis");
  passCiC_vs_passMVA->Draw("colz,text");

  c_overlap->SaveAs("passCiC_vs_passMVA.png");

  TCanvas *c_overlap_vbf = new TCanvas("c_overlap_vbf","passCiC vs passMVA: Dijet tagged",80,500,750,500);

  passVBFCiC_vs_passVBFMVA = (TH1*)passVBFCiC_vs_passVBFMVA_cat0_ggh_m125_8TeV->Clone();
  passVBFCiC_vs_passVBFMVA->Add(passVBFCiC_vs_passVBFMVA_cat0_vbf_m125_8TeV);
  passVBFCiC_vs_passVBFMVA->Add(passVBFCiC_vs_passVBFMVA_cat0_wzh_m125_8TeV);
  passVBFCiC_vs_passVBFMVA->Add(passVBFCiC_vs_passVBFMVA_cat0_tth_m125_8TeV);
  passVBFCiC_vs_passVBFMVA->SetBinContent(1,1,0.);
  passVBFCiC_vs_passVBFMVA->SetMarkerSize(4.);
  passVBFCiC_vs_passVBFMVA->GetXaxis()->SetNdivisions(2);
  passVBFCiC_vs_passVBFMVA->GetYaxis()->SetNdivisions(2);
  passVBFCiC_vs_passVBFMVA->GetXaxis()->SetTitleSize(0.05);
  passVBFCiC_vs_passVBFMVA->GetYaxis()->SetTitleSize(0.05);
  passVBFCiC_vs_passVBFMVA->GetXaxis()->SetTitle("Dijet tagged in MVA analysis");
  passVBFCiC_vs_passVBFMVA->GetYaxis()->SetTitle("Dijet tagged in cut based analysis");
  passVBFCiC_vs_passVBFMVA->Draw("colz,text");

  c_overlap_vbf->SaveAs("passVBFCiC_vs_passVBFMVA.png");

  TCanvas *c_overlap_cat = new TCanvas("c_overlap_cat","passCiC vs passMVA in categories",900,20,750,500);

  passCiC_vs_passMVA_cat = (TH1*)passCiC_vs_passMVA_cat_cat0_ggh_m125_8TeV->Clone();
  passCiC_vs_passMVA_cat->Add(passCiC_vs_passMVA_cat_cat0_vbf_m125_8TeV);
  passCiC_vs_passMVA_cat->Add(passCiC_vs_passMVA_cat_cat0_wzh_m125_8TeV);
  passCiC_vs_passMVA_cat->Add(passCiC_vs_passMVA_cat_cat0_tth_m125_8TeV);
  passCiC_vs_passMVA_cat->SetBinContent(1,1,0.);
  passCiC_vs_passMVA_cat->SetBinContent(2,1,0.);
  passCiC_vs_passMVA_cat->SetMarkerSize(2.);
  passCiC_vs_passMVA_cat->GetXaxis()->SetNdivisions(6);
  passCiC_vs_passMVA_cat->GetYaxis()->SetNdivisions(5);
  passCiC_vs_passMVA_cat->GetXaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA_cat->GetYaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA_cat->GetXaxis()->SetRangeUser(-1.5,3.5);
  passCiC_vs_passMVA_cat->GetXaxis()->SetTitle("Selected in MVA analysis");
  passCiC_vs_passMVA_cat->GetYaxis()->SetTitle("Selected in cut based analysis");
  passCiC_vs_passMVA_cat->Draw("colz,text");

  c_overlap_cat->SaveAs("passCiC_vs_passMVA_cat.png");

  TCanvas *c_passlevel = new TCanvas("c_passlevel","pass acceptance, presel and dipho mva, in diphoton categories",900,500,750,500);

  diphocat_vs_passMVA = (TH1*)diphocat_vs_passMVA_cat0_ggh_m125_8TeV->Clone();
  diphocat_vs_passMVA->Add(diphocat_vs_passMVA_cat0_vbf_m125_8TeV);
  diphocat_vs_passMVA->Add(diphocat_vs_passMVA_cat0_wzh_m125_8TeV);
  diphocat_vs_passMVA->Add(diphocat_vs_passMVA_cat0_tth_m125_8TeV);
  diphocat_vs_passMVA->SetMarkerSize(2.);
  diphocat_vs_passMVA->GetXaxis()->SetRangeUser(0.5,3.5);
  diphocat_vs_passMVA->GetXaxis()->SetNdivisions(3);
  diphocat_vs_passMVA->GetYaxis()->SetNdivisions(4);
  diphocat_vs_passMVA->GetXaxis()->SetTitleSize(0.05);
  diphocat_vs_passMVA->GetYaxis()->SetTitleSize(0.05);
  diphocat_vs_passMVA->GetXaxis()->SetTitle("MVA selection level");
  diphocat_vs_passMVA->GetYaxis()->SetTitle("Di-photon category");
  diphocat_vs_passMVA->Draw("colz,text");

  c_passlevel->SaveAs("diphocat_vs_passMVA.png");

}
