void overlap() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetPalette(1);
  
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);

  TFile *file = TFile::Open("histograms_CMS-HGG_categories.root");
  file->cd();

  TCanvas *c_overlap = new TCanvas("c_overlap","passCiC vs passMVA",750,500);

  passCiC_vs_passMVA_cat0_Data->SetBinContent(1,1,0.);
  passCiC_vs_passMVA_cat0_Data->SetMarkerSize(4.);
  passCiC_vs_passMVA_cat0_Data->Draw("colz,text");
  passCiC_vs_passMVA_cat0_Data->GetXaxis()->SetNdivisions(2);
  passCiC_vs_passMVA_cat0_Data->GetYaxis()->SetNdivisions(2);
  passCiC_vs_passMVA_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA_cat0_Data->GetYaxis()->SetTitleSize(0.05);
  passCiC_vs_passMVA_cat0_Data->GetXaxis()->SetTitle("Selected in MVA analysis");
  passCiC_vs_passMVA_cat0_Data->GetYaxis()->SetTitle("Selected in cut based analysis");

  c_overlap->SaveAs("passCiC_vs_passMVA.png");

  TCanvas *c_overlap_vbf = new TCanvas("c_overlap_vbf","passCiC vs passMVA: Dijet tagged",80,500,750,500);

  passVBFCiC_vs_passVBFMVA_cat0_Data->SetBinContent(1,1,0.);
  passVBFCiC_vs_passVBFMVA_cat0_Data->SetMarkerSize(4.);
  passVBFCiC_vs_passVBFMVA_cat0_Data->Draw("colz,text");
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetXaxis()->SetNdivisions(2);
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetYaxis()->SetNdivisions(2);
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetYaxis()->SetTitleSize(0.05);
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetXaxis()->SetTitle("Dijet tagged in MVA analysis");
  passVBFCiC_vs_passVBFMVA_cat0_Data->GetYaxis()->SetTitle("Dijet tagged in cut based analysis");

  c_overlap_vbf->SaveAs("passVBFCiC_vs_passVBFMVA.png");

  TCanvas *c_overlap_cat = new TCanvas("c_overlap_cat","passCiC vs passMVA in categories",900,20,750,500);

  passCiC_vs_passMVA_cat_cat0_Data->SetBinContent(1,1,0.);
  passCiC_vs_passMVA_cat_cat0_Data->SetBinContent(2,1,0.);
  passCiC_vs_passMVA_cat_cat0_Data->SetMarkerSize(2.);
  passCiC_vs_passMVA_cat_cat0_Data->Draw("colz,text");
  passCiC_vs_passMVA_cat_cat0_Data->GetXaxis()->SetNdivisions(6);
  passCiC_vs_passMVA_cat_cat0_Data->GetYaxis()->SetNdivisions(6);
  passCiC_vs_passMVA_cat_cat0_Data->GetXaxis()->SetRangeUser(-1.5,3.5);
  passCiC_vs_passMVA_cat_cat0_Data->GetXaxis()->SetTitle("Selected in MVA analysis");
  passCiC_vs_passMVA_cat_cat0_Data->GetYaxis()->SetTitle("Selected in cut based analysis");

  c_overlap_cat->SaveAs("passCiC_vs_passMVA_cat.png");


  TCanvas *c_passlevel = new TCanvas("c_passlevel","pass acceptance, presel and dipho mva, in diphoton categories",900,500,750,500);

  diphocat_vs_passMVA_cat0_Data->SetMarkerSize(2.);
  diphocat_vs_passMVA_cat0_Data->GetXaxis()->SetRangeUser(0.5,3.5);
  diphocat_vs_passMVA_cat0_Data->GetXaxis()->SetNdivisions(3);
  diphocat_vs_passMVA_cat0_Data->GetYaxis()->SetNdivisions(4);
  diphocat_vs_passMVA_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  diphocat_vs_passMVA_cat0_Data->GetYaxis()->SetTitleSize(0.05);
  diphocat_vs_passMVA_cat0_Data->GetXaxis()->SetTitle("MVA selection level");
  diphocat_vs_passMVA_cat0_Data->GetYaxis()->SetTitle("Di-photon category");
  diphocat_vs_passMVA_cat0_Data->Draw("colz,text");

  c_passlevel->SaveAs("diphocat_vs_passMVA.png");
}
