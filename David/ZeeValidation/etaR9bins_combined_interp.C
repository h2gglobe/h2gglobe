#include <algorithm>

void etaR9bins_combined_interp(TString r9bin, bool lead, bool idmva, bool combined) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TString combinedStr = "";
  if (!combined) combinedStr = "_posNegEta";

  TString labelStr = "#sigmaE/E";
  if (idmva) labelStr = "ID MVA output";

  TString lead_str = "1";
  if (!lead) lead_str="2";
  TString lead_str2="_lead";
  if (!lead) lead_str2="_sublead";

  TString r9_str;
  if (r9bin=="lowR9") {r9_str="R9<0.88";}
  else if (r9bin=="midR9") {r9_str="0.88<R9<0.94";}
  else if (r9bin=="highR9") {r9_str="R9>0.94";}
  else {
    cout << "1st argument must be lowR9, midR9 or highR9" << endl;
    return;
  }

  TFile *file_data = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  TFile *file_MC = TFile::Open("histograms_CMS-HGG_zeevalidation_envelopes.root");

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  TH1* hist_Data[8];
  TH1* hist_MC[8];
  TH1* hist_MC_up[8];
  TH1* hist_MC_down[8];
  TH1* hist_syst[8];

  TString eta_str[8] = {"-2.5<#eta<-2.0","-2.0<#eta<-1.5","-1.5<#eta<-0.9","-0.9<#eta<0","0<#eta<0.9","0.9<#eta<1.5","1.5<#eta<2.0","2.0<#eta<2.5"};

  for (int i=0; i<8; i++) {
    TString iStr;
    iStr+=i;
    if (!idmva) {
      hist_Data[i] = (TH1*)(file_data->Get("pho"+lead_str+"_sigmaEOverE_"+r9bin+"_cat"+iStr+"_Data"))->Clone();
      hist_MC[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_sigmaEOverE_"+r9bin+"_cat"+iStr+"_DYJetsToLL"))->Clone();
      hist_MC_up[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_sigmaEOverE_"+r9bin+"_cat"+iStr+"_DYJetsToLL_top"))->Clone();
      hist_MC_down[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_sigmaEOverE_"+r9bin+"_cat"+iStr+"_DYJetsToLL_bottom"))->Clone();
    } else {
      hist_Data[i] = (TH1*)(file_data->Get("pho"+lead_str+"_phoidMva_"+r9bin+"_cat"+iStr+"_Data"))->Clone();
      hist_MC[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_phoidMva_"+r9bin+"_cat"+iStr+"_DYJetsToLL"))->Clone();
      hist_MC_up[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_phoidMva_"+r9bin+"_cat"+iStr+"_DYJetsToLL_top"))->Clone();
      hist_MC_down[i] = (TH1*)(file_MC->Get("pho"+lead_str+"_phoidMva_"+r9bin+"_cat"+iStr+"_DYJetsToLL_bottom"))->Clone();
    }

//     if (r9bin!="highR9") {
//       hist_Data[i]->Rebin(2);
//       hist_MC[i]->Rebin(2);
//       hist_MC_up[i]->Rebin(2);
//       hist_MC_down[i]->Rebin(2);
//     }

    hist_MC_up[i]->SetLineColor(2);
    hist_MC_down[i]->SetLineColor(2);

    hist_Data[i]->GetYaxis()->SetTitle("");
    hist_Data[i]->GetXaxis()->SetTitle(labelStr+" ("+r9_str+", "+eta_str[i]+")");
    hist_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_MC[i]->SetFillColor(38);
    hist_MC_up[i]->SetMarkerStyle(0);
    hist_Data[i]->SetMarkerStyle(20);
    hist_Data[i]->SetMarkerSize(.5);
  }

  if (combined) {
    for (int i=0; i<4; i++) {
      hist_Data[i+4]->Add(hist_Data[3-i]);
      hist_MC[i+4]->Add(hist_MC[3-i]);
      hist_MC_up[i+4]->Add(hist_MC_up[3-i]);
      hist_MC_down[i+4]->Add(hist_MC_down[3-i]);
    }
  }

  for (int i=0; i<8; i++) {
    hist_syst[i] = (TH1*)hist_MC_up[i]->Clone();
    hist_syst[i]->Add(hist_MC_down[i]);
    hist_syst[i]->Scale(0.5);
    for (int j=1; j<hist_syst[i]->GetNbinsX()+1; j++) {
      float up = hist_MC_up[i]->GetBinContent(j);
      float down = hist_MC_down[i]->GetBinContent(j);
      hist_syst[i]->SetBinError(j,fabs(up-down)/2.);
    }
    hist_syst[i]->SetFillStyle(3013);
    hist_syst[i]->SetFillColor(2);
  }

  TLegend *leg2;
  leg2 = new TLegend(.5,.55,.87,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(hist_Data[0],"Data (19.6fb^{-1})");
  leg2->AddEntry(hist_MC[0],"DYJetsToLL MC","F");
  if (idmva) {
    leg2->AddEntry(hist_syst[0],"MC, "+labelStr+" scale factor ±0.1","F");
  } else {
    leg2->AddEntry(hist_syst[0],"MC, "+labelStr+" scale factor ±10%","F");
  }

  TLine *line;
  if (!idmva) {
    line = new TLine(0.,1.,0.06,1.);
  } else {
    line = new TLine(-0.2,1.,0.5,1.);
  }
  line->SetLineColor(4);
  line->SetLineWidth(2);

  int canvasHeight = combined ? 800:1600;
  int shift = combined ? 4:0;
  TCanvas *c = new TCanvas("c","",2200,canvasHeight);
  if (combined) {
    c->Divide(4,2);
  } else {
    c->Divide(4,4);
  }

  if (!combined) {
    c->cd(4);
    plotDataMC(hist_Data[0],hist_MC[0],hist_MC_up[0],hist_MC_down[0], hist_syst[0],idmva, r9bin);
    gPad->RedrawAxis();
    if (r9bin=="midR9" && !idmva) leg2->Draw();

    c->cd(3);
    plotDataMC(hist_Data[1],hist_MC[1],hist_MC_up[1],hist_MC_down[1], hist_syst[1],idmva, r9bin);
    gPad->RedrawAxis();
    if (r9bin=="midR9" && lead && !idmva) leg2->Draw();

    c->cd(2);
    plotDataMC(hist_Data[2],hist_MC[2],hist_MC_up[2],hist_MC_down[2], hist_syst[2],idmva, r9bin);
    gPad->RedrawAxis();
    if (r9bin!="lowR9" && !idmva) leg2->Draw();

    c->cd(1);
    plotDataMC(hist_Data[3],hist_MC[3],hist_MC_up[3],hist_MC_down[3], hist_syst[3],idmva, r9bin);
    gPad->RedrawAxis();
    if (!idmva) leg2->Draw();

    c->cd(8);
    gPad->SetGrid();
    plotRatio(hist_Data[0],hist_MC[0],hist_MC_up[0],hist_MC_down[0],idmva, r9bin);
    line->Draw();

    c->cd(7);
    gPad->SetGrid();
    plotRatio(hist_Data[1],hist_MC[1],hist_MC_up[1],hist_MC_down[1],idmva, r9bin);
    line->Draw();

    c->cd(6);
    gPad->SetGrid();
    plotRatio(hist_Data[2],hist_MC[2],hist_MC_up[2],hist_MC_down[2],idmva, r9bin);
    line->Draw();

    c->cd(5);
    gPad->SetGrid();
    plotRatio(hist_Data[3],hist_MC[3],hist_MC_up[3],hist_MC_down[3],idmva, r9bin);
    line->Draw();
  }

  c->cd(combined ? 1:9);
  plotDataMC(hist_Data[4],hist_MC[4],hist_MC_up[4],hist_MC_down[4], hist_syst[4],idmva, r9bin);
  gPad->RedrawAxis();
  if (!idmva) leg2->Draw();

  c->cd(combined ? 2:10);
  plotDataMC(hist_Data[5],hist_MC[5],hist_MC_up[5],hist_MC_down[5], hist_syst[5],idmva, r9bin);
  gPad->RedrawAxis();
  if (r9bin!="lowR9" && !idmva) leg2->Draw();

  c->cd(combined ? 3:11);
  plotDataMC(hist_Data[6],hist_MC[6],hist_MC_up[6],hist_MC_down[6], hist_syst[6],idmva, r9bin);
  gPad->RedrawAxis();
  if (r9bin=="midR9" && lead && !idmva) leg2->Draw();


//   TFile *fout = new TFile("sigmaEoverE_interpolation.root","RECREATE");
//   fout->cd();
//   hist_MC[6]->Write("sigmaEoverE");
//   hist_MC_up[6]->Write("sigmaEoverE_up");
//   hist_MC_down[6]->Write("sigmaEoverE_down");
//   fout->Close();


  c->cd(combined ? 4:12);
  plotDataMC(hist_Data[7],hist_MC[7],hist_MC_up[7],hist_MC_down[7], hist_syst[7],idmva, r9bin);
  gPad->RedrawAxis();
  if (r9bin=="midR9" && !idmva) leg2->Draw();

  c->cd(combined ? 5:13);
  gPad->SetGrid();
  plotRatio(hist_Data[4],hist_MC[4],hist_MC_up[4],hist_MC_down[4],idmva, r9bin);
  line->Draw();

  c->cd(combined ? 6:14);
  gPad->SetGrid();
  plotRatio(hist_Data[5],hist_MC[5],hist_MC_up[5],hist_MC_down[5],idmva, r9bin);
  line->Draw();

  c->cd(combined ? 7:15);
  gPad->SetGrid();
  plotRatio(hist_Data[6],hist_MC[6],hist_MC_up[6],hist_MC_down[6],idmva, r9bin);
  line->Draw();

  c->cd(combined ? 8:16);
  gPad->SetGrid();
  plotRatio(hist_Data[7],hist_MC[7],hist_MC_up[7],hist_MC_down[7],idmva, r9bin);
  line->Draw();

  if (!idmva) {
    c->SaveAs("sigmaE_"+r9bin+lead_str2+combinedStr+".png");
  } else {
    c->SaveAs("idmva_"+r9bin+lead_str2+combinedStr+".png");
  }

}

void plotDataMC(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, TH1* hist_syst, bool idmva, TString r9bin)
{

  float sf = hist_Data->Integral()/hist_MC->Integral();
  hist_MC->Scale(sf);
  hist_MC_down->Scale(sf);
  hist_MC_up->Scale(sf);
  hist_syst->Scale(sf);
  hist_MC_line = (TH1*)hist_MC->Clone();
  hist_MC_line->SetFillColor(0);
  if (!idmva) {
    if (r9bin!="highR9") {
      hist_Data->GetXaxis()->SetRangeUser(0.,0.06);
      hist_MC_up->GetXaxis()->SetRangeUser(0.,0.06);
    } else {
      hist_Data->GetXaxis()->SetRangeUser(0.,0.03);
      hist_MC_up->GetXaxis()->SetRangeUser(0.,0.03);
    }
  } else {
    hist_Data->GetXaxis()->SetRangeUser(-0.2,0.5);
    hist_MC_up->GetXaxis()->SetRangeUser(-0.2,0.5);
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

void plotRatio(TH1* hist_Data, TH1* hist_MC, TH1* hist_MC_up, TH1* hist_MC_down, bool idmva, TString r9bin)
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
  if (!idmva) {
    if (r9bin!="highR9") {
      ratio_syst->GetXaxis()->SetRangeUser(0.,0.06);
    } else {
      ratio_syst->GetXaxis()->SetRangeUser(0.,0.03);
    }
  } else {
    ratio_syst->GetXaxis()->SetRangeUser(-0.2,0.5);
  }
  ratio_syst->Draw("e2");
  ratio_up->Draw("hist,same");
  ratio_down->Draw("hist,same");
  ratio->SetLineColor(1);
  ratio->Draw("e,same");

}
