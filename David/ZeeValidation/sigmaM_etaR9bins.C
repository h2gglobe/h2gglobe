#include <algorithm>

void sigmaM_etaR9bins(bool wrongVtx=false, TString r9bin, bool preselNorm=false) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TString wrongVtx_str="";
  TString rightWrong="right";
  if (wrongVtx) {
    wrongVtx_str="wrongVtx_";
    rightWrong="wrong";
  }

  TString r9_str;
  if (r9bin=="lowR9") {r9_str="both R9<0.9";}
  else if (r9bin=="mixedR9") {r9_str="one R9<0.9, one R9>0.9";}
  else if (r9bin=="highR9") {r9_str="both R9>0.9";}
  else {
    cout << "1st argument must be lowR9, midR9 or highR9" << endl;
    return;
  }

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  TH1* hist_Data[4];
  TH1* hist_MC[4];

  TString eta_str[4] = {"Inclusive","EBEB","EBEE","EEEE"};
  
  float sf_presel = sigmaMOverM_cat0_Data->Integral()/sigmaMOverM_cat0_DYJetsToLL->Integral();

  for (int i=0; i<4; i++) {
    TString iStr;
    iStr+=i;
    hist_Data[i] = (TH1*)(file->Get("sigmaMOverM_"+wrongVtx_str+r9bin+"_cat"+iStr+"_Data"))->Clone();
    hist_MC[i] = (TH1*)(file->Get("sigmaMOverM_"+wrongVtx_str+r9bin+"_cat"+iStr+"_DYJetsToLL"))->Clone();

    //hist_Data[i]->Rebin(2);
    //hist_MC[i]->Rebin(2);

    hist_Data[i]->GetYaxis()->SetTitle("");
    hist_Data[i]->GetXaxis()->SetTitle("#sigma_{M}/M_{#gamma#gamma} ("+rightWrong+" vertex) ("+r9_str+", "+eta_str[i]+")");
    hist_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_MC[i]->SetFillColor(38);
    hist_Data[i]->SetMarkerStyle(20);
    hist_Data[i]->SetMarkerSize(.5);
  }

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_Data[0],"Data (19.6fb^{-1})");
  leg->AddEntry(hist_MC[0],"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(hist_Data[0],"Data (19.6fb^{-1})");
  leg2->AddEntry(hist_MC[0],"DYJetsToLL MC","F");

  TLine *line= new TLine(0.,1.,0.06,1.);
  line->SetLineColor(4);
  line->SetLineWidth(2);

  TCanvas *c = new TCanvas("c","",2200,800);
  c->Divide(4,2);

  c->cd(1);
  plotDataMC(hist_Data[0],hist_MC[0],preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(2);
  plotDataMC(hist_Data[1],hist_MC[1],preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(3);
  plotDataMC(hist_Data[2],hist_MC[2],preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(4);
  plotDataMC(hist_Data[3],hist_MC[3],preselNorm,sf_presel);
  gPad->RedrawAxis();
  leg->Draw();

  c->cd(5);
  gPad->SetGrid();
  plotRatio(hist_Data[0],hist_MC[0]);
  line->Draw();

  c->cd(6);
  gPad->SetGrid();
  plotRatio(hist_Data[1],hist_MC[1]);
  line->Draw();

  c->cd(7);
  gPad->SetGrid();
  plotRatio(hist_Data[2],hist_MC[2]);
  line->Draw();

  c->cd(8);
  gPad->SetGrid();
  plotRatio(hist_Data[3],hist_MC[3]);
  line->Draw();

  if (preselNorm) {
    c->SaveAs("sigmaMOverM_"+wrongVtx_str+r9bin+"_preselNorm.png");
  } else {
    c->SaveAs("sigmaMOverM_"+wrongVtx_str+r9bin+"_equalArea.png");
  }

}

void plotDataMC(TH1* hist_Data, TH1* hist_MC, bool preselNorm, float sf_presel)
{

  float sf = hist_Data->Integral()/hist_MC->Integral();
  if (preselNorm) sf=sf_presel;
  hist_MC->Scale(sf);
  float max = hist_Data->GetMaximum();
  hist_Data->SetMaximum(max*1.1);
  hist_Data->Draw("e");
  hist_MC->Draw("hist,same");
  hist_Data->Draw("e,same");

}

void plotRatio(TH1* hist_Data, TH1* hist_MC)
{

  ratio = (TH1*)hist_Data->Clone();
  ratio->Divide(hist_MC);
  ratio->SetLineColor(1);
  ratio->SetMaximum(1.8);
  ratio->SetMinimum(0.2);
  ratio->GetYaxis()->SetTitle("Data/MC Ratio");
  ratio->GetYaxis()->SetTitleSize(0.05);
  ratio->Draw("e");

}
