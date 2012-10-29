void rhoVsNvtx() {

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

  rhoVsNvtx_cat0_Data->GetXaxis()->SetTitle("Number of vertices");
  rhoVsNvtx_cat0_Data->GetYaxis()->SetTitle("rho");
  rhoVsNvtx_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  rhoVsNvtx_cat0_Data->GetYaxis()->SetTitleSize(0.05);
  rhoVsNvtx_cat0_Data->SetMarkerStyle(20);
  rhoVsNvtx_cat0_Data->SetLineColor(1);

  rhoVsNvtx_cat0_DYJetsToLL->GetXaxis()->SetTitle("Number of vertices");
  rhoVsNvtx_cat0_DYJetsToLL->GetYaxis()->SetTitle("rho");
  rhoVsNvtx_cat0_DYJetsToLL->GetXaxis()->SetTitleSize(0.05);
  rhoVsNvtx_cat0_DYJetsToLL->GetYaxis()->SetTitleSize(0.05);
  rhoVsNvtx_cat0_DYJetsToLL->SetMarkerStyle(20);
  rhoVsNvtx_cat0_DYJetsToLL->SetLineColor(38);
  rhoVsNvtx_cat0_DYJetsToLL->SetMarkerColor(38);

  nvtx_cat0_Data->GetXaxis()->SetTitle("Number of vertices");
  nvtx_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  nvtx_cat0_DYJetsToLL->SetFillColor(38);
  nvtx_cat0_DYJetsToLL->SetLineColor(1);
  nvtx_cat0_Data->SetMarkerStyle(20);
  nvtx_cat0_Data->SetMarkerSize(0.4);
  nvtx_cat0_Data->SetLineColor(1);

  rho_EB_cat0_Data->GetXaxis()->SetTitle("rho");
  rho_EB_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  rho_EB_cat0_DYJetsToLL->SetFillColor(38);
  rho_EB_cat0_DYJetsToLL->SetLineColor(1);
  rho_EB_cat0_Data->SetMarkerStyle(20);
  rho_EB_cat0_Data->SetMarkerSize(0.4);
  rho_EB_cat0_Data->SetLineColor(1);

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(nvtx_cat0_Data,"Data (12.2fb^{-1})");
  leg->AddEntry(nvtx_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(rhoVsNvtx_cat0_Data,"Data (12.2fb^{-1})");
  leg2->AddEntry(rhoVsNvtx_cat0_DYJetsToLL,"DYJetsToLL MC");

  TLine *line = new TLine(-1.,1.,1.,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);

  TCanvas *c_rhoVsNvtx = new TCanvas("c_rhoVsNvtx","rho vs nVtx",1500,900);
  c_rhoVsNvtx->Divide(3,2);

  c_rhoVsNvtx->cd(1);
  sf = nvtx_cat0_Data->Integral()/nvtx_cat0_DYJetsToLL->Integral();
  nvtx_cat0_DYJetsToLL->Scale(sf);
  nvtx_cat0_Data->Draw("e");
  nvtx_cat0_DYJetsToLL->Draw("hist,same");
  nvtx_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_rhoVsNvtx->cd(2);
  sf = rho_EB_cat0_Data->Integral()/rho_EB_cat0_DYJetsToLL->Integral();
  rho_EB_cat0_DYJetsToLL->Scale(sf);
  rho_EB_cat0_Data->Draw("e");
  rho_EB_cat0_DYJetsToLL->Draw("hist,same");
  rho_EB_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw("colz");

  c_rhoVsNvtx->cd(4);
  rhoVsNvtx_cat0_Data->Draw("colz");
  TText *t = new TText(0.2, 0.8, "Data");
  t->SetNDC();
  t->Draw();

  c_rhoVsNvtx->cd(5);
  rhoVsNvtx_cat0_DYJetsToLL->Draw("colz");
  TText *t = new TText(0.2, 0.8, "DYJetsToLL MC");
  t->SetNDC();
  t->Draw();

  c_rhoVsNvtx->cd(6);
  rhoVsNvtx_cat0_DYJetsToLL->ProfileX()->Draw();
  rhoVsNvtx_cat0_Data->ProfileX()->Draw("same");
  leg2->Draw();

  c_rhoVsNvtx->SaveAs("rhoVsNvtx.png");

}
