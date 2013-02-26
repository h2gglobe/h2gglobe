void ptRescale(bool setWeights=false) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");

  pt_cat0_Data_rebin = (TH1*)pt_cat0_Data->Clone();
  pt_cat0_DYJetsToLL_rebin = (TH1*)pt_cat0_DYJetsToLL->Clone();

  if (setWeights) {
    Double_t xbins[21] = {0.,2.,4.,8.,12.,16.,20.,24.,28.,32.,36.,40.,50.,60.,70.,80.,95.,110.,130.,160.,200.};
    pt_cat0_Data_rebin->Rebin(20,"pt_cat0_Data_rebin",xbins);
    pt_cat0_DYJetsToLL_rebin->Rebin(20,"pt_cat0_DYJetsToLL_rebin",xbins);
  }

  //pt_cat0_Data_rebin->Rebin(4);
  //pt_cat0_DYJetsToLL_rebin->Rebin(4);

  pt_cat0_Data_rebin->GetXaxis()->SetTitleSize(0.05);
  pt_cat0_Data_rebin->GetYaxis()->SetTitle("");
  pt_cat0_Data_rebin->GetXaxis()->SetTitle("p_{T}(e^{+}e^{-}) (GeV)");
  pt_cat0_DYJetsToLL_rebin->SetFillColor(38);
  pt_cat0_DYJetsToLL_rebin->SetLineColor(1);
  pt_cat0_Data_rebin->SetMarkerStyle(20);
  pt_cat0_Data_rebin->SetMarkerSize(0.4);
  pt_cat0_Data_rebin->SetLineColor(1);

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(pt_cat0_Data_rebin,"Data (19.6fb^{-1})");
  leg->AddEntry(pt_cat0_DYJetsToLL_rebin,"DYJetsToLL MC","F");

  TLine *line = new TLine(0.,1.,200.,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);

  TCanvas *c_pt = new TCanvas("c_pt","di-photon p_{T}",600,800);
  c_pt->Divide(1,2);

  c_pt->cd(1);
  gPad->SetLogy();
  float sf = pt_cat0_Data_rebin->Integral()/pt_cat0_DYJetsToLL_rebin->Integral();
  pt_cat0_DYJetsToLL_rebin->Scale(sf);
  pt_cat0_Data_rebin->Draw("e");
  pt_cat0_DYJetsToLL_rebin->Draw("hist,same");
  pt_cat0_Data_rebin->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_pt->cd(2);
  gPad->SetGrid();
  ratio_pt = (TH1*)pt_cat0_Data_rebin->Clone();
  ratio_pt->Divide(pt_cat0_DYJetsToLL_rebin);
  //ratio_pt->SetMaximum(1.8);
  //ratio_pt->SetMinimum(0.2);
  ratio_pt->Draw("e");
  line->Draw();

  if (setWeights) {
    for (int i=1;i<21;i++) cout << "pt_reweight[" << i-1 << "] = " << ratio_pt->GetBinContent(i) << ";" << endl;
    cout << endl;
    for (int i=1;i<21;i++) cout << "pt_reweight_bin[" << i-1 << "] = " << ratio_pt->GetBinLowEdge(i) << ".;" << endl;
    cout << endl;
    cout << "zeePtWeight=";
    for (int i=1;i<21;i++) cout << ratio_pt->GetBinContent(i) << ",";
    cout << endl;
  }

  c_pt->SaveAs("pt.png");


  nvtx_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  nvtx_cat0_Data->GetYaxis()->SetTitle("");
  nvtx_cat0_Data->GetXaxis()->SetTitle("number of primary vertices");
  nvtx_cat0_Data->GetXaxis()->SetRangeUser(0.,40.);
  nvtx_cat0_DYJetsToLL->SetFillColor(38);
  nvtx_cat0_DYJetsToLL->SetLineColor(1);
  nvtx_cat0_Data->SetMarkerStyle(20);
  nvtx_cat0_Data->SetMarkerSize(0.4);
  nvtx_cat0_Data->SetLineColor(1);

  TLine *line1 = new TLine(0.,1.,40.,1.);
  line1->SetLineColor(4);
  line1->SetLineWidth(2);

  TCanvas *c_nvtx = new TCanvas("c_nvtx","no. of primmary vertices",800,900);
  c_nvtx->Divide(1,2);

  c_nvtx->cd(1);
  sf = nvtx_cat0_Data->Integral()/nvtx_cat0_DYJetsToLL->Integral();
  nvtx_cat0_DYJetsToLL->Scale(sf);
  nvtx_cat0_Data->GetXaxis()->SetRangeUser(0.,39.5);
  nvtx_cat0_Data->Draw("e");
  nvtx_cat0_DYJetsToLL->Draw("hist,same");
  nvtx_cat0_Data->Draw("e,same");
  gPad->RedrawAxis();
  leg->Draw();

  c_nvtx->cd(2);
  gPad->SetGrid();
  ratio_nvtx = (TH1*)nvtx_cat0_Data->Clone();
  ratio_nvtx->Divide(nvtx_cat0_DYJetsToLL);
  ratio_nvtx->SetMaximum(1.8);
  ratio_nvtx->SetMinimum(0.2);
  ratio_nvtx->Draw("e");
  line1->Draw();

  //for (int i=1;i<26;i++) cout << "nvtx_reweight[" << i-1 << "] = " << ratio_nvtx->GetBinContent(i) << ";" << endl;
  //cout << endl;

  c_nvtx->SaveAs("nvtx.png");

}
