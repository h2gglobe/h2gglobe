void soverb_ratio() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetPalette(1);
  
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetPaintTextFormat("6.3f");

  TFile *file_sig = TFile::Open("histograms_CMS-HGG_categories_sig.root");
  file_sig->cd();
  hist_sig = (TH1*)passCiC_vs_passMVA_cat_cat0_tot->Clone();

  TFile *file = TFile::Open("histograms_CMS-HGG_categories.root");
  file->cd();
  hist_data = (TH1*)passCiC_vs_passMVA_cat_cat0_Data->Clone();

  TCanvas *c_overlap_cat_ratio = new TCanvas("c_overlap_cat_ratio","passCiC vs passMVA in categories: S/B",900,20,750,500);

  hist_sig->SetBinContent(1,1,0.);
  hist_sig->SetBinContent(2,1,0.);
  hist_sig->SetMarkerSize(2.);
  hist_sig->GetXaxis()->SetNdivisions(6);
  hist_sig->GetYaxis()->SetNdivisions(6);
  hist_sig->GetXaxis()->SetTitle("Selected in HIG-12-001");
  hist_sig->GetYaxis()->SetTitle("Selected in HIG-11-033");
  hist_sig->Divide(hist_data);
  hist_sig->Draw("colz,text");

  c_overlap_cat_ratio->SaveAs("passCiC_vs_passMVA_cat_ratio.png");
}
