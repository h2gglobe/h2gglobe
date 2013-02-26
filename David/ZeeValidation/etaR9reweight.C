void etaR9reweight() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  gStyle->SetPalette(1);

  TFile *file_zee = TFile::Open("../../AnalysisScripts/zee_8Feb/zee_phoPD_runABCD/histograms_CMS-HGG.root");

  zee_lead = (TH2*)pho1_r9VsEta_cat0_DYJetsToLL->Clone();
  zee_sublead = (TH2*)pho2_r9VsEta_cat0_DYJetsToLL->Clone();

  TFile *file_hgg = TFile::Open("../../AnalysisScripts/zee_sig/histograms_CMS-HGG_zeevalidation.root");
  
  hgg_lead = (TH2*)pho1_r9VsEta_cat0_tot->Clone();
  hgg_sublead = (TH2*)pho2_r9VsEta_cat0_tot->Clone();

  zee_lead->RebinX(2);
  zee_sublead->RebinX(2);
  hgg_lead->RebinX(2);
  hgg_sublead->RebinX(2);

  Double_t r9bins[14] = {0.,0.5,0.7,0.76,0.82,0.86,0.9,0.92,0.94,0.95,0.96,0.97,0.98,1.};
  TH2F* zee_lead_varbins = new TH2F("zee_lead_varbins", "zee_lead_varbins", 25, -2.5, 2.5, 13, r9bins);
  TH2F* zee_sublead_varbins = new TH2F("zee_sublead_varbins", "zee_sublead_varbins", 25, -2.5, 2.5, 13, r9bins);
  TH2F* hgg_lead_varbins = new TH2F("hgg_lead_varbins", "hgg_lead_varbins", 25, -2.5, 2.5, 13, r9bins);
  TH2F* hgg_sublead_varbins = new TH2F("hgg_sublead_varbins", "hgg_sublead_varbins", 25, -2.5, 2.5, 13, r9bins);

  for (int i=1; i<26; i++) {
    for (int j=1; j<101; j++) {
      zee_lead_varbins->Fill(zee_lead->GetXaxis()->GetBinCenter(i),zee_lead->GetYaxis()->GetBinCenter(j),zee_lead->GetBinContent(i,j));
      zee_sublead_varbins->Fill(zee_sublead->GetXaxis()->GetBinCenter(i),zee_sublead->GetYaxis()->GetBinCenter(j),zee_sublead->GetBinContent(i,j));
      hgg_lead_varbins->Fill(hgg_lead->GetXaxis()->GetBinCenter(i),hgg_lead->GetYaxis()->GetBinCenter(j),hgg_lead->GetBinContent(i,j));
      hgg_sublead_varbins->Fill(hgg_sublead->GetXaxis()->GetBinCenter(i),hgg_sublead->GetYaxis()->GetBinCenter(j),hgg_sublead->GetBinContent(i,j));
    }
  }

  hgg_lead->GetXaxis()->SetTitleSize(0.05);
  hgg_lead->GetYaxis()->SetTitleSize(0.05);
  hgg_lead->GetYaxis()->SetTitle("lead R9");
  hgg_lead->GetXaxis()->SetTitle("lead #eta");
  hgg_sublead->GetXaxis()->SetTitleSize(0.05);
  hgg_sublead->GetYaxis()->SetTitleSize(0.05);
  hgg_sublead->GetYaxis()->SetTitle("sublead R9");
  hgg_sublead->GetXaxis()->SetTitle("sublead #eta");

  zee_lead->GetXaxis()->SetTitleSize(0.05);
  zee_lead->GetYaxis()->SetTitleSize(0.05);
  zee_lead->GetYaxis()->SetTitle("lead R9");
  zee_lead->GetXaxis()->SetTitle("lead #eta");
  zee_sublead->GetXaxis()->SetTitleSize(0.05);
  zee_sublead->GetYaxis()->SetTitleSize(0.05);
  zee_sublead->GetYaxis()->SetTitle("sublead R9");
  zee_sublead->GetXaxis()->SetTitle("sublead #eta"); 
  hgg_lead_varbins->GetXaxis()->SetTitleSize(0.05);
  hgg_lead_varbins->GetYaxis()->SetTitleSize(0.05);
  hgg_lead_varbins->GetYaxis()->SetTitle("lead R9");
  hgg_lead_varbins->GetXaxis()->SetTitle("lead #eta");
  hgg_sublead_varbins->GetXaxis()->SetTitleSize(0.05);
  hgg_sublead_varbins->GetYaxis()->SetTitleSize(0.05);
  hgg_sublead_varbins->GetYaxis()->SetTitle("sublead R9");
  hgg_sublead_varbins->GetXaxis()->SetTitle("sublead #eta");

  zee_both =(TH2*)zee_lead->Clone();
  zee_both->Add(zee_sublead);
  hgg_both =(TH2*)hgg_lead->Clone();
  hgg_both->Add(hgg_sublead);
  zee_both_varbins =(TH2*)zee_lead_varbins->Clone();
  zee_both_varbins->Add(zee_sublead_varbins);
  hgg_both_varbins =(TH2*)hgg_lead_varbins->Clone();
  hgg_both_varbins->Add(hgg_sublead_varbins);

  zee_both->GetYaxis()->SetTitle("R9");
  zee_both->GetXaxis()->SetTitle("#eta");
  hgg_both->GetYaxis()->SetTitle("R9");
  hgg_both->GetXaxis()->SetTitle("#eta");
  zee_both_varbins->GetYaxis()->SetTitle("R9");
  zee_both_varbins->GetXaxis()->SetTitle("#eta");
  hgg_both_varbins->GetYaxis()->SetTitle("R9");
  hgg_both_varbins->GetXaxis()->SetTitle("#eta");

  TCanvas *c_lead = new TCanvas("c_lead","eta vs R9, lead",1200,900);
  c_lead->Divide(2,2);

  c_lead->cd(1);
  gPad->SetLogz();
  zee_lead->Draw("colz");
  gPad->RedrawAxis();

  c_lead->cd(2);
  gPad->SetLogz();
  hgg_lead->Draw("colz");
  gPad->RedrawAxis();

  c_lead->cd(3);
  gPad->SetLogz();
  float sf = hgg_lead->Integral()/zee_lead->Integral();
  zee_lead_varbins->Scale(sf);
  ratio_lead = (TH2*)hgg_lead_varbins->Clone();
  ratio_lead->Divide(zee_lead_varbins);
  //ratio_lead->SetMaximum(1.8);
  //ratio_lead->SetMinimum(0.2);

  for (int i=1; i<26; i++) {
    for (int j=1; j<14; j++) {
      if (ratio_lead->GetBinContent(i,j)==0.) ratio_lead->SetBinContent(i,j,1.);
    }
  }

  ratio_lead->Draw("colz");

  //c_lead->SaveAs("etaR9reweight.png");


  TCanvas *c_sublead = new TCanvas("c_sublead","eta vs R9, sublead",1200,900);
  c_sublead->Divide(2,2);

  c_sublead->cd(1);
  gPad->SetLogz();
  zee_sublead->Draw("colz");
  gPad->RedrawAxis();

  c_sublead->cd(2);
  gPad->SetLogz();
  hgg_sublead->Draw("colz");
  gPad->RedrawAxis();

  c_sublead->cd(3);
  gPad->SetLogz();
  float sf = hgg_sublead->Integral()/zee_sublead->Integral();
  zee_sublead_varbins->Scale(sf);
  ratio_sublead = (TH2*)hgg_sublead_varbins->Clone();
  ratio_sublead->Divide(zee_sublead_varbins);
  //ratio_sublead->SetMaximum(1.8);
  //ratio_sublead->SetMinimum(0.2);

  for (int i=1; i<26; i++) {
    for (int j=1; j<14; j++) {
      if (ratio_sublead->GetBinContent(i,j)==0.) ratio_sublead->SetBinContent(i,j,1.);
    }
  }

  ratio_sublead->Draw("colz");


  //c_sublead->SaveAs("etaR9reweight.png");
 

  TCanvas *c_both = new TCanvas("c_both","eta vs R9, both",1200,900);
  c_both->Divide(2,2);

  c_both->cd(1);
  gPad->SetLogz();
  zee_both->Draw("colz");
  gPad->RedrawAxis();

  c_both->cd(2);
  gPad->SetLogz();
  hgg_both->Draw("colz");
  gPad->RedrawAxis();

  c_both->cd(3);
  gPad->SetLogz();
  float sf = hgg_both->Integral()/zee_both->Integral();
  zee_both_varbins->Scale(sf);
  ratio_both = (TH2*)hgg_both_varbins->Clone();
  ratio_both->Divide(zee_both_varbins);
  //ratio_both->SetMaximum(1.8);
  //ratio_both->SetMinimum(0.2);

  for (int i=1; i<26; i++) {
    for (int j=1; j<14; j++) {
      if (ratio_both->GetBinContent(i,j)==0.) ratio_both->SetBinContent(i,j,1.);
    }
  }

  ratio_both->Draw("colz");

  c_both->cd(4);
  ratio_leadsublead = (TH2*)zee_lead->Clone();
  ratio_leadsublead->GetYaxis()->SetTitle("R9");
  ratio_leadsublead->GetXaxis()->SetTitle("#eta");
  ratio_leadsublead->Divide(zee_sublead);
  ratio_leadsublead->SetMaximum(1.5);
  ratio_leadsublead->SetMinimum(0.5);
  ratio_leadsublead->Draw("colz");

  //c_both->SaveAs("etaR9reweight.png");

  TFile *fout = new TFile("zee_etaR9reweight.root","RECREATE");
  fout->cd();
  ratio_lead->Write("etaR9reweight_lead");
  ratio_sublead->Write("etaR9reweight_sublead");
  fout->Close();

}
