#include <algorithm>

//argument var can be "bdtout", "idmva_EB" or "idmva_EE"

void diphomva_nvtx(TString var="bdtout") {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);
  gStyle->SetLineColor(1);

  TString nvtxbin_str[3] = {"0to13","14to18","19up"};
  TString nvtxbin_label[3] = {"0 <= nvtx <= 13","14 <= nvtx <= 18","nvtx >= 19"};

  TH1* hist_Data[3];
  TH1* hist_MC[3];
  TH1* hist_MC_up[3];
  TH1* hist_MC_down[3];
  TH1* hist_syst[3];

  TFile *file_data = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file_data->cd();

  for (int nvtxbin=0; nvtxbin<3; nvtxbin++) {
    hist_Data[nvtxbin] = (TH1*)(file_data->Get(var+"_nvtx"+nvtxbin_str[nvtxbin]+"_cat0_Data"))->Clone();
  }

  TFile *file_MC = TFile::Open("histograms_CMS-HGG_zeevalidation_envelopes.root");
  file_MC->cd();

  for (int nvtxbin=0; nvtxbin<3; nvtxbin++) {

    hist_MC[nvtxbin] = (TH1*)(file_MC->Get(var+"_nvtx"+nvtxbin_str[nvtxbin]+"_cat0_DYJetsToLL"))->Clone();
    hist_MC_up[nvtxbin] = (TH1*)(file_MC->Get(var+"_nvtx"+nvtxbin_str[nvtxbin]+"_cat0_DYJetsToLL_top"))->Clone();
    hist_MC_down[nvtxbin] = (TH1*)(file_MC->Get(var+"_nvtx"+nvtxbin_str[nvtxbin]+"_cat0_DYJetsToLL_bottom"))->Clone();

    if (var=="bdtout") {
      hist_MC[nvtxbin]->Rebin(2);
      hist_Data[nvtxbin]->Rebin(2);
      hist_MC_up[nvtxbin]->Rebin(2);
      hist_MC_down[nvtxbin]->Rebin(2);
      hist_Data[nvtxbin]->GetXaxis()->SetTitle("diphoton BDT output ("+nvtxbin_label[nvtxbin]+")");
    } else if (var=="idmva_EB"){
      hist_Data[nvtxbin]->GetXaxis()->SetTitle("photon ID MVA output (barrel, "+nvtxbin_label[nvtxbin]+")");
    } else {
      hist_Data[nvtxbin]->GetXaxis()->SetTitle("photon ID MVA output (endcaps, "+nvtxbin_label[nvtxbin]+")");
    }
    hist_Data[nvtxbin]->GetYaxis()->SetTitle("");
    hist_Data[nvtxbin]->GetXaxis()->SetTitleSize(0.045);
    hist_MC[nvtxbin]->SetFillColor(38);
    hist_MC[nvtxbin]->SetLineColor(1);
    hist_MC_up[nvtxbin]->SetLineColor(2);
    hist_MC_up[nvtxbin]->SetMarkerStyle(0);
    hist_MC_down[nvtxbin]->SetLineColor(2);
    hist_Data[nvtxbin]->SetMarkerStyle(20);
    hist_Data[nvtxbin]->SetMarkerSize(.7);
    hist_Data[nvtxbin]->SetLineColor(1);

    if (var=="idmva_EB" || var=="idmva_EE") hist_Data[nvtxbin]->GetXaxis()->SetRangeUser(-0.2,0.45);

    hist_syst[nvtxbin] = (TH1*)hist_MC_up[nvtxbin]->Clone();
    hist_syst[nvtxbin]->Add(hist_MC_down[nvtxbin]);
    hist_syst[nvtxbin]->Scale(0.5);
    for (int i=1; i<hist_syst[nvtxbin]->GetNbinsX()+1; i++) {
      float up = hist_MC_up[nvtxbin]->GetBinContent(i);
      float down = hist_MC_down[nvtxbin]->GetBinContent(i);
      hist_syst[nvtxbin]->SetBinError(i,fabs(up-down)/2.);
    }

  }

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.045);

  TLegend *leg;
  if (var=="bdtout") {
    leg = new TLegend(.15,.65,.4,.87);
  } else {
    leg = new TLegend(.12,.65,.4,.87);
  }
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_Data[0],"Data");
  leg->AddEntry(hist_MC[0],"Z#rightarrow e^{+}e^{-} MC","F");
  leg->AddEntry(hist_syst[0],"MC systematic uncertainty","F");

  TCanvas *c = new TCanvas("c","BDT output",2200,500);
  c->Divide(3,1);

  for (int nvtxbin=0; nvtxbin<3; nvtxbin++) {

    c->cd(nvtxbin+1);

    float sf = hist_Data[nvtxbin]->Integral()/hist_MC[nvtxbin]->Integral();
    hist_MC[nvtxbin]->Scale(sf);
    hist_MC_up[nvtxbin]->Scale(sf);
    hist_MC_down[nvtxbin]->Scale(sf);
    hist_syst[nvtxbin]->Scale(sf);
    hist_Data[nvtxbin]->Draw("e");
    leg->Draw();
    if (var=="bdtout") {
      txt->DrawLatex(0.45,0.8,"#scale[0.8]{#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");
    } else {
      txt->DrawLatex(0.63,0.8,"#scale[0.8]{#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 19.6 fb^{-1}}}");
    }

    hist_MC[nvtxbin]->Draw("hist,same");
    hist_syst[nvtxbin]->SetFillStyle(3013);
    hist_syst[nvtxbin]->SetFillColor(2);
    hist_syst[nvtxbin]->Draw("e2,same");
    hist_MC_up[nvtxbin]->Draw("hist,same");
    hist_MC_down[nvtxbin]->Draw("hist,same");
    hist_Data[nvtxbin]->Draw("e,same");
    gPad->RedrawAxis();

    if (var=="bdtout") {
      float eff_data = hist_Data[nvtxbin]->Integral(96,200)/hist_Data[nvtxbin]->Integral();
      float eff_mc = hist_MC[nvtxbin]->Integral(96,200)/hist_MC[nvtxbin]->Integral();
      cout << nvtxbin_str[nvtxbin] << ": " << hist_Data[nvtxbin]->GetBinLowEdge(96) << " " << eff_data << " " << eff_mc << " " << eff_data/eff_mc << endl;
    }

  }

  if (var=="bdtout") {
    c->SaveAs("diphomva_nvtx.png");
    c->SaveAs("diphomva_nvtx.pdf");
  } else {
    c->SaveAs(var+"_nvtx.png");
    c->SaveAs(var+"_nvtx.pdf");
  }

}
