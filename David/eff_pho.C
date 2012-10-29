void eff_pho() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);

  TFile *f = TFile::Open("histograms_CMS-HGG_cic.root");

  TH1* hist_eta_presel[4][2];
  hist_eta_presel[0][0] = (TH1*)pho1_eta_presel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[1][0] = (TH1*)pho1_eta_presel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[2][0] = (TH1*)pho1_eta_presel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[3][0] = (TH1*)pho1_eta_presel_cat4_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[0][1] = (TH1*)pho2_eta_presel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[1][1] = (TH1*)pho2_eta_presel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[2][1] = (TH1*)pho2_eta_presel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_presel[3][1] = (TH1*)pho2_eta_presel_cat4_gluglu_H_gg_124_pu2012->Clone();
  TH1* hist_eta_sel[4][2];
  hist_eta_sel[0][0] = (TH1*)pho1_eta_sel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[1][0] = (TH1*)pho1_eta_sel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[2][0] = (TH1*)pho1_eta_sel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[3][0] = (TH1*)pho1_eta_sel_cat4_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[0][1] = (TH1*)pho2_eta_sel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[1][1] = (TH1*)pho2_eta_sel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[2][1] = (TH1*)pho2_eta_sel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_eta_sel[3][1] = (TH1*)pho2_eta_sel_cat4_gluglu_H_gg_124_pu2012->Clone();
  TH1* hist_pt_presel[4][2];
  hist_pt_presel[0][0] = (TH1*)pho1_pt_presel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[1][0] = (TH1*)pho1_pt_presel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[2][0] = (TH1*)pho1_pt_presel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[3][0] = (TH1*)pho1_pt_presel_cat4_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[0][1] = (TH1*)pho2_pt_presel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[1][1] = (TH1*)pho2_pt_presel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[2][1] = (TH1*)pho2_pt_presel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_presel[3][1] = (TH1*)pho2_pt_presel_cat4_gluglu_H_gg_124_pu2012->Clone();
  TH1* hist_pt_sel[4][2];
  hist_pt_sel[0][0] = (TH1*)pho1_pt_sel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[1][0] = (TH1*)pho1_pt_sel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[2][0] = (TH1*)pho1_pt_sel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[3][0] = (TH1*)pho1_pt_sel_cat4_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[0][1] = (TH1*)pho2_pt_sel_cat1_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[1][1] = (TH1*)pho2_pt_sel_cat2_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[2][1] = (TH1*)pho2_pt_sel_cat3_gluglu_H_gg_124_pu2012->Clone();
  hist_pt_sel[3][1] = (TH1*)pho2_pt_sel_cat4_gluglu_H_gg_124_pu2012->Clone();

  TH1* hist_eta_presel_all[4];
  TH1* hist_eta_sel_all[4];
  TH1* hist_pt_presel_all[4];
  TH1* hist_pt_sel_all[4];
  for (int icat = 0; icat<4; icat++) {
    hist_eta_presel_all[icat] = (TH1*)hist_eta_presel[icat][0]->Clone();
    hist_eta_presel_all[icat]->Add(hist_eta_presel[icat][1]);
    hist_eta_sel_all[icat] = (TH1*)hist_eta_sel[icat][0]->Clone();
    hist_eta_sel_all[icat]->Add(hist_eta_sel[icat][1]);
    hist_pt_presel_all[icat] = (TH1*)hist_pt_presel[icat][0]->Clone();
    hist_pt_presel_all[icat]->Add(hist_pt_presel[icat][1]);
    hist_pt_sel_all[icat] = (TH1*)hist_pt_sel[icat][0]->Clone();
    hist_pt_sel_all[icat]->Add(hist_pt_sel[icat][1]);
  }

  hist_eta_presel_allcat = (TH1F*)hist_eta_presel_all[0]->Clone();
  hist_eta_presel_allcat->Add(hist_eta_presel_all[1]);
  hist_eta_presel_allcat->Add(hist_eta_presel_all[2]);
  hist_eta_presel_allcat->Add(hist_eta_presel_all[3]);

  hist_eta_sel_allcat = (TH1F*)hist_eta_sel_all[0]->Clone();
  hist_eta_sel_allcat->Add(hist_eta_sel_all[1]);
  hist_eta_sel_allcat->Add(hist_eta_sel_all[2]);
  hist_eta_sel_allcat->Add(hist_eta_sel_all[3]);

  hist_pt_presel_allcat = (TH1F*)hist_pt_presel_all[0]->Clone();
  hist_pt_presel_allcat->Add(hist_pt_presel_all[1]);
  hist_pt_presel_allcat->Add(hist_pt_presel_all[2]);
  hist_pt_presel_allcat->Add(hist_pt_presel_all[3]);

  hist_pt_sel_allcat = (TH1F*)hist_pt_sel_all[0]->Clone();
  hist_pt_sel_allcat->Add(hist_pt_sel_all[1]);
  hist_pt_sel_allcat->Add(hist_pt_sel_all[2]);
  hist_pt_sel_allcat->Add(hist_pt_sel_all[3]);

  float N_eta_presel_bin[51][4][2];
  float N_eta_sel_bin[51][4][2];
  float N_pt_presel_bin[45][4][2];
  float N_pt_sel_bin[45][4][2];
  for (int ipho = 0; ipho<2; ipho++) {
    for (int icat = 0; icat<4; icat++) {
      for (int ibin=0; ibin<51; ibin++) {
	N_eta_presel_bin[ibin][icat][ipho] = hist_eta_presel[icat][ipho]->GetEntries()/hist_eta_presel[icat][ipho]->Integral()*hist_eta_presel[icat][ipho]->GetBinContent(ibin+1);
	N_eta_sel_bin[ibin][icat][ipho] = hist_eta_sel[icat][ipho]->GetEntries()/hist_eta_sel[icat][ipho]->Integral()*hist_eta_sel[icat][ipho]->GetBinContent(ibin+1);
      }
      for (int ibin=0; ibin<45; ibin++) {
	N_pt_presel_bin[ibin][icat][ipho] = hist_pt_presel[icat][ipho]->GetEntries()/hist_pt_presel[icat][ipho]->Integral()*hist_pt_presel[icat][ipho]->GetBinContent(ibin+1);
	N_pt_sel_bin[ibin][icat][ipho] = hist_pt_sel[icat][ipho]->GetEntries()/hist_pt_sel[icat][ipho]->Integral()*hist_pt_sel[icat][ipho]->GetBinContent(ibin+1);
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1","Canvas1");
  c1->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_eta_sel[icat][0]->Divide(hist_eta_presel[icat][0]);
    cout << endl;
    for (int ibin=0; ibin<51; ibin++) {
      float err = eff_err(N_eta_sel_bin[ibin][icat][0],N_eta_presel_bin[ibin][icat][0]);
      cout << ibin << " " << N_eta_sel_bin[ibin][icat][0] << " " << N_eta_presel_bin[ibin][icat][0] << " " << err << endl;
      hist_eta_sel[icat][0]->SetBinError(ibin+1,err);
    }
    hist_eta_sel[icat][0]->SetMarkerSize(1.2);
  }
  hist_eta_sel[0][0]->SetMarkerColor(kGreen);
  hist_eta_sel[0][0]->SetMarkerStyle(20);
  hist_eta_sel[0][0]->SetMaximum(1.);
  hist_eta_sel[0][0]->SetMinimum(0.);
  hist_eta_sel[0][0]->GetXaxis()->SetTitle("leading photon SC #eta");
  hist_eta_sel[0][0]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_eta_sel[0][0]->Draw("e");
  hist_eta_sel[1][0]->SetMarkerColor(kCyan);
  hist_eta_sel[1][0]->SetMarkerStyle(21);
  hist_eta_sel[1][0]->Draw("same,e");
  hist_eta_sel[2][0]->SetMarkerColor(kBlue);
  hist_eta_sel[2][0]->SetMarkerStyle(22);
  hist_eta_sel[2][0]->Draw("same,e");
  hist_eta_sel[3][0]->SetMarkerColor(kRed);
  hist_eta_sel[3][0]->SetMarkerStyle(23);
  hist_eta_sel[3][0]->Draw("same,e");

  leg = new TLegend(.35,.15,.7,.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(hist_eta_sel[0][0],"category 1 (EB, R_{9}>0.94)","LP");
  leg->AddEntry(hist_eta_sel[1][0],"category 2 (EB, R_{9}<0.94)","LP");
  leg->AddEntry(hist_eta_sel[2][0],"category 3 (EE, R_{9}>0.94)","LP");
  leg->AddEntry(hist_eta_sel[3][0],"category 4 (EE, R_{9}<0.94)","LP");

  leg->Draw();

  txt = new TText();
  txt->SetNDC();
  txt->SetTextSize(0.04);
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");

  txt1 = new TLatex();
  txt1->SetNDC();
  txt1->SetTextSize(0.04);
  txt1->SetTextAlign(12);
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c2 = new TCanvas("c2","Canvas2");
  c2->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_pt_sel[icat][0]->Divide(hist_pt_presel[icat][0]);
    for (int ibin=0; ibin<45; ibin++) {
      hist_pt_sel[icat][0]->SetBinError(ibin+1,eff_err(N_pt_sel_bin[ibin][icat][0],N_pt_presel_bin[ibin][icat][0]));
    }
    hist_pt_sel[icat][0]->SetMarkerSize(1.2);
  }
  hist_pt_sel[0][0]->SetMarkerColor(kGreen);
  hist_pt_sel[0][0]->SetMarkerStyle(20);
  hist_pt_sel[0][0]->SetMaximum(1.);
  hist_pt_sel[0][0]->SetMinimum(0.);
  hist_pt_sel[0][0]->GetXaxis()->SetRangeUser(40.,100.);
  hist_pt_sel[0][0]->GetXaxis()->SetTitle("leading photon p_{T}");
  hist_pt_sel[0][0]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_pt_sel[0][0]->Draw("e");
  hist_pt_sel[1][0]->SetMarkerColor(kCyan);
  hist_pt_sel[1][0]->SetMarkerStyle(21);
  hist_pt_sel[1][0]->Draw("same,e");
  hist_pt_sel[2][0]->SetMarkerColor(kBlue);
  hist_pt_sel[2][0]->SetMarkerStyle(22);
  hist_pt_sel[2][0]->Draw("same,e");
  hist_pt_sel[3][0]->SetMarkerColor(kRed);
  hist_pt_sel[3][0]->SetMarkerStyle(23);
  hist_pt_sel[3][0]->Draw("same,e");
  leg->Draw();
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c3 = new TCanvas("c3","Canvas3");
  c3->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_eta_sel[icat][1]->Divide(hist_eta_presel[icat][1]);
    for (int ibin=0; ibin<51; ibin++) {
      hist_eta_sel[icat][1]->SetBinError(ibin+1,eff_err(N_eta_sel_bin[ibin][icat][1],N_eta_presel_bin[ibin][icat][1]));
    }
    hist_eta_sel[icat][1]->SetMarkerSize(1.2);
  }
  hist_eta_sel[0][1]->SetMarkerColor(kGreen);
  hist_eta_sel[0][1]->SetMarkerStyle(20);
  hist_eta_sel[0][1]->SetMaximum(1.);
  hist_eta_sel[0][1]->SetMinimum(0.);
  hist_eta_sel[0][1]->GetXaxis()->SetTitle("subleading photon SC #eta");
  hist_eta_sel[0][1]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_eta_sel[0][1]->Draw("e");
  hist_eta_sel[1][1]->SetMarkerColor(kCyan);
  hist_eta_sel[1][1]->SetMarkerStyle(21);
  hist_eta_sel[1][1]->Draw("same,e");
  hist_eta_sel[2][1]->SetMarkerColor(kBlue);
  hist_eta_sel[2][1]->SetMarkerStyle(22);
  hist_eta_sel[2][1]->Draw("same,e");
  hist_eta_sel[3][1]->SetMarkerColor(kRed);
  hist_eta_sel[3][1]->SetMarkerStyle(23);
  hist_eta_sel[3][1]->Draw("same,e");
  leg->Draw();
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c4 = new TCanvas("c4","Canvas4");
  c4->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_pt_sel[icat][1]->Divide(hist_pt_presel[icat][1]);
    for (int ibin=0; ibin<45; ibin++) {
      hist_pt_sel[icat][1]->SetBinError(ibin+1,eff_err(N_pt_sel_bin[ibin][icat][1],N_pt_presel_bin[ibin][icat][1]));
    }
    hist_pt_sel[icat][1]->SetMarkerSize(1.2);
  }
  hist_pt_sel[0][1]->SetMarkerColor(kGreen);
  hist_pt_sel[0][1]->SetMarkerStyle(20);
  hist_pt_sel[0][1]->SetMaximum(1.);
  hist_pt_sel[0][1]->SetMinimum(0.);
  hist_pt_sel[0][1]->GetXaxis()->SetRangeUser(30.,80.);
  hist_pt_sel[0][1]->GetXaxis()->SetTitle("subleading photon p_{T} (GeV/c");
  hist_pt_sel[0][1]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_pt_sel[0][1]->Draw("e");
  hist_pt_sel[1][1]->SetMarkerColor(kCyan);
  hist_pt_sel[1][1]->SetMarkerStyle(21);
  hist_pt_sel[1][1]->Draw("same,e");
  hist_pt_sel[2][1]->SetMarkerColor(kBlue);
  hist_pt_sel[2][1]->SetMarkerStyle(22);
  hist_pt_sel[2][1]->Draw("same,e");
  hist_pt_sel[3][1]->SetMarkerColor(kRed);
  hist_pt_sel[3][1]->SetMarkerStyle(23);
  hist_pt_sel[3][1]->Draw("same,e");
  leg->Draw();
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c5 = new TCanvas("c5","Canvas5");
  c5->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_eta_sel_all[icat]->Divide(hist_eta_presel_all[icat]);
    for (int ibin=0; ibin<51; ibin++) {
      float err = eff_err(N_eta_sel_bin[ibin][icat][0]+N_eta_sel_bin[ibin][icat][1],N_eta_presel_bin[ibin][icat][0]+N_eta_presel_bin[ibin][icat][1]);
      hist_eta_sel_all[icat]->SetBinError(ibin+1,err);
    }
    hist_eta_sel_all[icat]->SetMarkerSize(1.2);
  }
  hist_eta_sel_all[0]->SetMarkerColor(kGreen);
  hist_eta_sel_all[0]->SetMarkerStyle(20);
  hist_eta_sel_all[0]->SetMaximum(1.);
  hist_eta_sel_all[0]->SetMinimum(0.);
  hist_eta_sel_all[0]->GetXaxis()->SetTitle("photon SC #eta");
  hist_eta_sel_all[0]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_eta_sel_all[0]->Draw("e");
  hist_eta_sel_all[1]->SetMarkerColor(kCyan);
  hist_eta_sel_all[1]->SetMarkerStyle(21);
  hist_eta_sel_all[1]->Draw("same,e");
  hist_eta_sel_all[2]->SetMarkerColor(kBlue);
  hist_eta_sel_all[2]->SetMarkerStyle(22);
  hist_eta_sel_all[2]->Draw("same,e");
  hist_eta_sel_all[3]->SetMarkerColor(kRed);
  hist_eta_sel_all[3]->SetMarkerStyle(23);
  hist_eta_sel_all[3]->Draw("same,e");
  leg->Draw();
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c6 = new TCanvas("c6","Canvas6");
  c6->SetGrid();

  for (int icat = 0; icat<4; icat++) {
    hist_pt_sel_all[icat]->Divide(hist_pt_presel_all[icat]);
    for (int ibin=0; ibin<45; ibin++) {
      hist_pt_sel_all[icat]->SetBinError(ibin+1,eff_err(N_pt_sel_bin[ibin][icat][0]+N_pt_sel_bin[ibin][icat][1],N_pt_presel_bin[ibin][icat][0]+N_pt_presel_bin[ibin][icat][1]));
    }
    hist_pt_sel_all[icat]->SetMarkerSize(1.2);
  }
  hist_pt_sel_all[0]->SetMarkerColor(kGreen);
  hist_pt_sel_all[0]->SetMarkerStyle(20);
  hist_pt_sel_all[0]->SetMaximum(1.);
  hist_pt_sel_all[0]->SetMinimum(0.);
  hist_pt_sel_all[0]->GetXaxis()->SetRangeUser(30.,100.);
  hist_pt_sel_all[0]->GetXaxis()->SetTitle("photon p_{T} (GeV/c)");
  hist_pt_sel_all[0]->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_pt_sel_all[0]->Draw("e");
  hist_pt_sel_all[1]->SetMarkerColor(kCyan);
  hist_pt_sel_all[1]->SetMarkerStyle(21);
  hist_pt_sel_all[1]->Draw("same,e");
  hist_pt_sel_all[2]->SetMarkerColor(kBlue);
  hist_pt_sel_all[2]->SetMarkerStyle(22);
  hist_pt_sel_all[2]->Draw("same,e");
  hist_pt_sel_all[3]->SetMarkerColor(kRed);
  hist_pt_sel_all[3]->SetMarkerStyle(23);
  hist_pt_sel_all[3]->Draw("same,e");
  leg->Draw();
  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.45,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c7 = new TCanvas("c7","Canvas7");
  c7->SetGrid();

  hist_eta_sel_allcat->Divide(hist_eta_presel_allcat);

  for (int ibin=0; ibin<51; ibin++) {
    float ntot_sel=0.;
    float ntot_presel=0.;
    for (int ipho = 0; ipho<2; ipho++) {
      for (int icat = 0; icat<4; icat++) {
	ntot_sel+= N_eta_sel_bin[ibin][icat][ipho];
	ntot_presel+= N_eta_presel_bin[ibin][icat][ipho];
      }
    }
    float err = eff_err(ntot_sel,ntot_presel);
    hist_eta_sel_allcat->SetBinError(ibin+1,err);
  }
  hist_eta_sel_allcat->SetMarkerSize(1.);
  hist_eta_sel_allcat->SetMarkerSize(1.);

  hist_eta_sel_allcat->SetMarkerStyle(20);
  hist_eta_sel_allcat->SetMaximum(1.);
  hist_eta_sel_allcat->SetMinimum(0.);
  hist_eta_sel_allcat->GetXaxis()->SetTitle("photon SC #eta");
  hist_eta_sel_allcat->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_eta_sel_allcat->Draw("e");
  hist_eta_sel_allcat->SetMarkerStyle(20);
  hist_eta_sel_allcat->Draw("same,e");

  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.38,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  TCanvas *c8 = new TCanvas("c8","Canvas8");
  c8->SetGrid();

  hist_pt_sel_allcat->Divide(hist_pt_presel_allcat);

  for (int ibin=0; ibin<45; ibin++) {
    float ntot_sel=0.;
    float ntot_presel=0.;
    for (int ipho = 0; ipho<2; ipho++) {
      for (int icat = 0; icat<4; icat++) {
	ntot_sel+= N_pt_sel_bin[ibin][icat][ipho];
	ntot_presel+= N_pt_presel_bin[ibin][icat][ipho];
      }
    }
    float err = eff_err(ntot_sel,ntot_presel);
    hist_pt_sel_allcat->SetBinError(ibin+1,err);
  }
  hist_pt_sel_allcat->SetMarkerSize(1.);
  hist_pt_sel_allcat->SetMarkerSize(1.);

  hist_pt_sel_allcat->SetMarkerStyle(20);
  hist_pt_sel_allcat->SetMaximum(1.);
  hist_pt_sel_allcat->SetMinimum(0.);
  hist_pt_sel_allcat->GetXaxis()->SetTitle("photon p_{T} (GeV/c)");
  hist_pt_sel_allcat->GetYaxis()->SetTitle("Photon ID Efficiency");
  hist_pt_sel_allcat->Draw("e");
  hist_pt_sel_allcat->SetMarkerStyle(20);
  hist_pt_sel_allcat->Draw("same,e");

  txt->DrawText(0.55,0.92,"CMS Simulation Preliminary");
  txt1->DrawLatex(0.38,0.45,"gg #rightarrow H, m_{H} = 124 GeV");

  c1->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_eta_lead.png");
  c2->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_pt_lead.png");
  c3->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_eta_sublead.png");
  c4->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_pt_sublead.png");
  c5->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_eta.png");
  c6->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_pt.png");
  c7->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_eta_all.png");
  c8->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/unblinding/CiC_Jun14/phoEff/eff_pt_all.png");

}

float eff_err(float a, float b) {
  if (b==0.) return 0.;

  float eff=a/b;
  return sqrt(eff*(1-eff)/b);

}
