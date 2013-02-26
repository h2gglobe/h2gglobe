void correlations(TString var) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TString xvar,yvar;
  if (var=="pho1_idmvaVsEta") {
    xvar = "#eta(lead)";
    yvar = "ID MVA output (lead)";
  } else if (var=="pho2_idmvaVsEta") {
    xvar = "#eta (sublead)";
    yvar = "ID MVA output (sublead)";
  } else if (var=="pho1_idmvaVsPtOverM_EB" || var=="pho1_idmvaVsPtOverM_EE") {
    xvar = "lead p_{T}/m_{#gamma#gamma}";
    yvar = "ID MVA output (lead)";
  } else if (var=="pho2_idmvaVsPtOverM_EB" || var=="pho2_idmvaVsPtOverM_EE") {
    xvar = "sublead p_{T}/m_{#gamma#gamma}";
    yvar = "ID MVA output (sublead)";
  } else if (var=="pho1_sigmaEOverEVsEta") {
    xvar = "#eta (lead)";
    yvar = "#sigmaE/E (lead)";
  } else if (var=="pho2_sigmaEOverEVsEta") {
    xvar = "#eta (sublead)";
    yvar = "#sigmaE/E (sublead)";
  } else if (var=="pho1_sigmaEOverEVsIdmva_EB" || var=="pho1_sigmaEOverEVsIdmva_EE") {
    xvar = "ID MVA output (lead)";
    yvar = "#sigmaE/E (lead)";
  } else if (var=="pho2_sigmaEOverEVsIdmva_EB" || var=="pho2_sigmaEOverEVsIdmva_EE") {
    xvar = "ID MVA output (sublead)";
    yvar = "#sigmaE/E (sublead)";
  } else if (var=="pho1_sigmaEOverEVsPtOverM_EB" || var=="pho1_sigmaEOverEVsPtOverM_EE") {
    xvar = "lead p_{T}/m_{#gamma#gamma}";
    yvar = "#sigmaE/E (lead)";
  } else if (var=="pho2_sigmaEOverEVsPtOverM_EB" || var=="pho2_sigmaEOverEVsPtOverM_EE") {
    xvar = "sublead p_{T}/m_{#gamma#gamma}";
    yvar = "#sigmaE/E (sublead)";
  } else if (var=="idmvaVsRho_EB" || var=="idmvaVsRho_EE") {
    xvar = "Rho";
    yvar = "ID MVA output";
  } else if (var=="pfphotoniso03VsRho_EB" || var=="pfphotoniso03VsRho_EE") {
    xvar = "Rho";
    yvar = "PF photon isolation";
  } else {
    cout << "unknown histogram " << var << endl;
  }

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  hist_Data = (TH2*)(file->Get(var+"_cat0_Data"))->Clone();
  hist_MC = (TH2*)(file->Get(var+"_cat0_DYJetsToLL"))->Clone();

  hist_Data->GetXaxis()->SetTitle(xvar);
  hist_Data->GetYaxis()->SetTitle(yvar);
  hist_Data->GetXaxis()->SetTitleSize(0.05);
  hist_Data->GetYaxis()->SetTitleSize(0.05);
  hist_Data->SetMarkerStyle(20);
  hist_Data->SetLineColor(1);

  hist_MC->GetXaxis()->SetTitle(xvar);
  hist_MC->GetYaxis()->SetTitle(yvar);
  hist_MC->GetXaxis()->SetTitleSize(0.05);
  hist_MC->GetYaxis()->SetTitleSize(0.05);
  hist_MC->SetMarkerStyle(20);
  hist_MC->SetLineColor(38);
  hist_MC->SetMarkerColor(38);

  TLegend *leg;
  leg = new TLegend(.6,.9,.87,1.);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(hist_MC,"DYJetsToLL MC");

  TCanvas *c = new TCanvas("c","",1200,900);
  c->Divide(2,2);

  c->cd(1);
  if (var=="pfphotoniso03VsRho_EB" || var=="pfphotoniso03VsRho_EE") gPad->SetLogz();
  hist_Data->Draw("colz");
  TText *t = new TText(0.2, .93, "Data");
  t->SetNDC();
  t->Draw();

  c->cd(2);
  if (var=="pfphotoniso03VsRho_EB" || var=="pfphotoniso03VsRho_EE") gPad->SetLogz();
  hist_MC->Draw("colz");
  TText *t = new TText(0.2, .93, "DY MC");
  t->SetNDC();
  t->Draw();

  c->cd(3);
  hist_MC->ProfileX()->Draw();
  hist_Data->ProfileX()->Draw("same");
  TText *t = new TText(0.2, .93, "ProfileX");
  t->SetNDC();
  t->Draw();
  leg->Draw();

  c->cd(4);
  hist_MC->ProfileY()->Draw();
  hist_Data->ProfileY()->Draw("same");
  TText *t = new TText(0.2, .93, "ProfileY");
  t->SetNDC();
  t->Draw();
  leg->Draw();

  c->SaveAs(var+".png");

}
