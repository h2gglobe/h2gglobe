#include <algorithm>

void bdtout_etaR9bins(TString r9bin, bool preselNorm=false, bool passcut=false) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TString r9_str;
  if (r9bin=="lowR9") {r9_str="both R9<0.9";}
  else if (r9bin=="mixedR9") {r9_str="one R9<0.9, one R9>0.9";}
  else if (r9bin=="highR9") {r9_str="both R9>0.9";}
  else {
    cout << "1st argument must be lowR9, midR9 or highR9" << endl;
    return;
  }

  TString passcut_str="";
  if (passcut) passcut_str="_passcut";

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  TH1* hist_Data[4];
  TH1* hist_MC[4];
  TH1* hist_MC_up[4];
  TH1* hist_MC_down[4];
  TH1* hist_syst[4];

  TString eta_str[4] = {"Inclusive","EBEB","EBEE","EEEE"};
  
  float sf_presel = bdtout_cat0_Data->Integral()/bdtout_cat0_DYJetsToLL->Integral();

  for (int i=0; i<4; i++) {
    TString iStr;
    iStr+=i;
    hist_Data[i] = (TH1*)(file->Get("bdtout_"+r9bin+"_cat"+iStr+"_Data"))->Clone();
    hist_MC[i] = (TH1*)(file->Get("bdtout_"+r9bin+"_cat"+iStr+"_DYJetsToLL"))->Clone();
    hist_MC_up[i] = (TH1*)(file->Get("bdtout_"+r9bin+"_up_cat"+iStr+"_DYJetsToLL"))->Clone();
    hist_MC_down[i] = (TH1*)(file->Get("bdtout_"+r9bin+"_down_cat"+iStr+"_DYJetsToLL"))->Clone();

    if (passcut) {
      hist_Data[i]->Rebin(2);
      hist_MC[i]->Rebin(2);
      hist_MC_up[i]->Rebin(2);
      hist_MC_down[i]->Rebin(2);
    } else {
      hist_Data[i]->Rebin(4);
      hist_MC[i]->Rebin(4);
      hist_MC_up[i]->Rebin(4);
      hist_MC_down[i]->Rebin(4);
    }

    hist_MC_up[i]->SetLineColor(2);
    hist_MC_down[i]->SetLineColor(2);

    hist_Data[i]->GetYaxis()->SetTitle("");
    hist_Data[i]->GetXaxis()->SetTitle("diphoMVA output ("+r9_str+", "+eta_str[i]+")");
    hist_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_MC[i]->SetFillColor(38);
    hist_Data[i]->SetMarkerStyle(20);
    hist_Data[i]->SetMarkerSize(.5);

    hist_syst[i] = (TH1*)hist_MC_up[i]->Clone();
    hist_syst[i]->Add(hist_MC_down[i]);
    hist_syst[i]->Scale(0.5);
    for (int j=1; j<hist_syst[i]->GetNbinsX()+1; j++) {
      float up = hist_MC_up[i]->GetBinContent(j);
      float down = hist_MC_down[i]->GetBinContent(j);
      hist_syst[i]->SetBinError(j,fabs(up-down)/2.);
    }
  }

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_Data[0],"Data (19.6fb^{-1})");
  leg->AddEntry(hist_MC[0],"DYJetsToLL MC","F");
  leg->AddEntry(hist_syst[0],"MC with idmva±0.01","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(hist_Data[0],"Data (19.6fb^{-1})");
  leg2->AddEntry(hist_MC[0],"DYJetsToLL MC","F");
  leg2->AddEntry(hist_syst[0],"MC with idmva±0.01","F");

  TLine *line;
  if (passcut) {
    line = new TLine(-0.05,1.,1.,1.);
  } else {
    line = new TLine(-1.,1.,1.,1.);
  }
  line->SetLineColor(4);
  line->SetLineWidth(2);

  TCanvas *c = new TCanvas("c","",2200,800);
  c->Divide(4,2);

  c->cd(1);
  plotDataMC(hist_Data[0],hist_MC[0],hist_MC_up[0],hist_MC_down[0], hist_syst[0],passcut,preselNorm,sf_presel);
  gPad->RedrawAxis();
  if (r9bin=="highR9" || !passcut) leg2->Draw();

  c->cd(2);
  plotDataMC(hist_Data[1],hist_MC[1],hist_MC_up[1],hist_MC_down[1], hist_syst[1],passcut,preselNorm,sf_presel);
  gPad->RedrawAxis();
  if (r9bin=="highR9" || !passcut) leg2->Draw();

  c->cd(3);
  plotDataMC(hist_Data[2],hist_MC[2],hist_MC_up[2],hist_MC_down[2], hist_syst[2],passcut,preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(4);
  plotDataMC(hist_Data[3],hist_MC[3],hist_MC_up[3],hist_MC_down[3], hist_syst[3],passcut,preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(5);
  gPad->SetGrid();
  plotRatio(hist_Data[0],hist_MC[0],hist_MC_up[0],hist_MC_down[0],passcut);
  line->Draw();

  c->cd(6);
  gPad->SetGrid();
  plotRatio(hist_Data[1],hist_MC[1],hist_MC_up[1],hist_MC_down[1],passcut);
  line->Draw();

  c->cd(7);
  gPad->SetGrid();
  plotRatio(hist_Data[2],hist_MC[2],hist_MC_up[2],hist_MC_down[2],passcut);
  line->Draw();

  c->cd(8);
  gPad->SetGrid();
  plotRatio(hist_Data[3],hist_MC[3],hist_MC_up[3],hist_MC_down[3],passcut);
  line->Draw();

  if (preselNorm) {
    c->SaveAs("bdtout_"+r9bin+passcut_str+"_preselNorm.png");
  } else {
    c->SaveAs("bdtout_"+r9bin+passcut_str+"_equalArea.png");
  }

}

void plotDataMC(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, TH1* hist_syst, bool passcut, bool preselNorm, float sf_presel)
{

  float sf;
  if (passcut) {
    sf = hist_Data->Integral(48,100)/hist_MC->Integral(48,100);
  } else {
    sf = hist_Data->Integral()/hist_MC->Integral();
  }
  if (preselNorm) sf=sf_presel;
  hist_MC->Scale(sf);
  hist_MC_down->Scale(sf);
  hist_MC_up->Scale(sf);
  hist_syst->Scale(sf);
  hist_MC_line = (TH1*)hist_MC->Clone();
  hist_MC_line->SetFillColor(0);
  if (passcut) {
    hist_Data->GetXaxis()->SetRangeUser(-0.05,1.);
    hist_MC_up->GetXaxis()->SetRangeUser(-0.05,1.);
  }
  float max = hist_Data->GetMaximum();
  hist_Data->SetMaximum(max*1.1);
  hist_Data->Draw("e");
  hist_MC->Draw("hist,same");
  hist_syst->SetFillStyle(3013);
  hist_syst->SetFillColor(2);
  hist_syst->Draw("e2,same");
  hist_MC_up->Draw("hist,same");
  hist_MC_down->Draw("hist,same");
  hist_MC_line->Draw("hist,same");
  hist_Data->Draw("e,same");

}

void plotRatio(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, bool passcut)
{

  ratio = (TH1*)hist_Data->Clone();
  ratio_up = (TH1*)hist_MC_up->Clone();
  ratio_down = (TH1*)hist_MC_down->Clone();
  ratio->Divide(hist_MC);
  ratio_up->Divide(hist_MC);
  ratio_down->Divide(hist_MC);
  ratio_up->SetLineColor(2);
  ratio_down->SetLineColor(2);

  ratio_syst = (TH1*)ratio_up->Clone();
  ratio_syst->Add(ratio_down);
  ratio_syst->Scale(0.5);
  for (int i=1; i<ratio_syst->GetNbinsX()+1; i++) {
    float up = max(ratio_up->GetBinContent(i),ratio_down->GetBinContent(i));
    float down = min(ratio_up->GetBinContent(i),ratio_down->GetBinContent(i));
    ratio_syst->SetBinError(i,fabs(up-down)/2.);
  }

  ratio_syst->SetFillStyle(3013);
  ratio_syst->SetFillColor(2);
  ratio_syst->SetMarkerStyle(0);
  ratio_syst->SetMaximum(1.8);
  ratio_syst->SetMinimum(0.2);
  ratio_syst->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio_syst->GetYaxis()->SetTitleSize(0.05);
  if (passcut) ratio_syst->GetXaxis()->SetRangeUser(-0.05,1.);
  ratio_syst->Draw("e2");
  ratio_up->Draw("hist,same");
  ratio_down->Draw("hist,same");
  ratio->Draw("e,same");

}
