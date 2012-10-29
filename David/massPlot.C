void massPlot(TString category="0", TString var="all_mass", bool blind=true, bool breakdown=false) {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetTitleFont(42,"XYZ");

  gStyle->SetMarkerSize(0.8);
  gStyle->SetTitleYOffset(1.5);

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);

  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");

  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.24);

  TFile *f = TFile::Open("histograms_CMS-HGG_cic.root");

  TH1 *hist_mass_sig;
  TH1 *hist_mass_data;
  TH1 *hist_mass_born;
  TH1 *hist_mass_box;
  TH1 *hist_mass_gjet_pf;
  TH1 *hist_mass_gjet20_pf;
  TH1 *hist_mass_qcd_pf;
  TH1 *hist_mass_qcd30_pf;
  TH1 *hist_mass_qcd_ff;
  TH1 *hist_mass_qcd30_ff;
  TH1 *hist_mass_dy;
  TH1 *hist_mass_bkg;
  THStack* hist_mass_bkg_stack;

  hist_mass_data = (TH1*)(f->Get(var+"_cat"+category+"_Data"))->Clone();
  hist_mass_sig = (TH1*)(f->Get(var+"_cat"+category+"_ggh_m125_8TeV"))->Clone();
  hist_mass_sig->Add((TH1*)(f->Get(var+"_cat"+category+"_vbf_m125_8TeV")));
  hist_mass_sig->Add((TH1*)(f->Get(var+"_cat"+category+"_wzh_m125_8TeV")));
  hist_mass_sig->Add((TH1*)(f->Get(var+"_cat"+category+"_tth_m125_8TeV")));
  hist_mass_born = (TH1*)(f->Get(var+"_cat"+category+"_DiPhotonJets"))->Clone();
  hist_mass_box = (TH1*)(f->Get(var+"_cat"+category+"_DiPhotonBox25"))->Clone();

  hist_mass_gjet_pf = (TH1*)(f->Get(var+"_cat"+category+"_GJet40PF"))->Clone();
  hist_mass_gjet20_pf = (TH1*)(f->Get(var+"_cat"+category+"_GJet20PF"))->Clone();
  hist_mass_gjet_pf->Add(hist_mass_gjet20_pf);

  hist_mass_qcd_ff = (TH1*)(f->Get(var+"_cat"+category+"_QCD40FF"))->Clone();
  hist_mass_qcd30_ff = (TH1*)(f->Get(var+"_cat"+category+"_QCD30FF"))->Clone();
  hist_mass_qcd_ff->Add(hist_mass_qcd30_ff);

  hist_mass_qcd_pf = (TH1*)(f->Get(var+"_cat"+category+"_QCD40PF"))->Clone();
  hist_mass_qcd30_pf = (TH1*)(f->Get(var+"_cat"+category+"_QCD30PF"))->Clone();
  hist_mass_qcd_pf->Add(hist_mass_qcd30_pf);

  hist_mass_dy = (TH1*)(f->Get(var+"_cat"+category+"_DYJetsToLL"))->Clone();

  if (var=="all_mass" || var=="pt") {
    hist_mass_data->Rebin(2);
    hist_mass_sig->Rebin(2);
    hist_mass_born->Rebin(2);
    hist_mass_box->Rebin(2);
    hist_mass_gjet_pf->Rebin(2);
    hist_mass_qcd_pf->Rebin(2);
    hist_mass_qcd_ff->Rebin(2);
    hist_mass_dy->Rebin(2);
  }

  float lumiSF=1.;
  hist_mass_born->Scale(lumiSF);
  hist_mass_box->Scale(lumiSF*18/17);
  hist_mass_gjet_pf->Scale(lumiSF);
  hist_mass_qcd_pf->Scale(lumiSF);
  hist_mass_qcd_ff->Scale(lumiSF);
  hist_mass_dy->Scale(lumiSF);
  hist_mass_sig->Scale(lumiSF);

  //fix DY normalization
  //hist_mass_dy->Scale(2.089*3532.8/(2950.*1.15));
  //hist_mass_dy->Scale((2950.*1.15)/(3532.8));
  //hist_mass_dy->Scale(1./2.089);

//   cout << hist_mass_born->Integral() << ", ";
//   cout << hist_mass_box->Integral() << ", ";
//   cout << hist_mass_gjet_pf->Integral() << ", ";
//   cout << hist_mass_qcd_pf->Integral() << ", ";
//   cout << hist_mass_qcd_ff->Integral() << ", ";
//   cout << hist_mass_dy->Integral() << endl;

  hist_mass_pp = (TH1*)hist_mass_born->Clone();
  hist_mass_pp->Add(hist_mass_box);
  hist_mass_pf = (TH1*)hist_mass_gjet_pf->Clone();
  hist_mass_pf->Add(hist_mass_qcd_pf);
  hist_mass_ff = (TH1*)hist_mass_qcd_ff->Clone();

  if (var=="all_mass") {
    float n_pp = hist_mass_pp->Integral(31,70);
    float n_pf = hist_mass_pf->Integral(31,70);
    float n_ff = hist_mass_ff->Integral(31,70);
    float frac_pp = n_pp/(n_pp+n_pf+n_ff);
    cout << "frac_pp = " << frac_pp << " " << hist_mass_pp->GetXaxis()->GetBinLowEdge(31) << " " << hist_mass_pp->GetXaxis()->GetBinLowEdge(71) << endl;

    float n_bkg=0.;
    n_bkg += hist_mass_pp->Integral(31,50);
    n_bkg += hist_mass_pf->Integral(31,50);
    n_bkg += hist_mass_ff->Integral(31,50);
    n_bkg += hist_mass_dy->Integral(31,50);
    cout << "Total background (110-130): " << n_bkg << endl;
  }

  hist_mass_sig->SetLineColor(2);
  hist_mass_sig->SetLineWidth(2.5);
  hist_mass_data->SetMarkerStyle(20);
  hist_mass_data->SetMarkerSize(0.8);
  hist_mass_pp->SetFillColor(kBlue-8);
  hist_mass_pf->SetFillColor(kGreen-10);
  hist_mass_ff->SetFillColor(kCyan-10);
  hist_mass_dy->SetFillColor(kMagenta-10);
  //hist_mass_pp->SetFillColor(591);
  //hist_mass_pf->SetFillColor(406);
  //hist_mass_ff->SetFillColor(422);
  //hist_mass_dy->SetFillColor(606);
  hist_mass_pp->SetLineWidth(2.5);
  hist_mass_pf->SetLineWidth(2.5);
  hist_mass_ff->SetLineWidth(2.5);
  hist_mass_dy->SetLineWidth(2.5);

  hist_mass_born->SetFillColor(kGreen-2);
  hist_mass_box->SetFillColor(kGreen-1);
  hist_mass_gjet_pf->SetFillColor(kOrange-2);
  hist_mass_qcd_pf->SetFillColor(kOrange-3);
  hist_mass_qcd_ff->SetFillColor(kOrange+2);

  float bkg_tot=0.;
  bkg_tot += hist_mass_pp->Integral();
  bkg_tot += hist_mass_pf->Integral();
  bkg_tot += hist_mass_ff->Integral();
  bkg_tot += hist_mass_dy->Integral();

  //float normFac = hist_mass_data->Integral()/bkg_tot;
  //float normFac = hist_mass_dy_s4->Integral()/hist_mass_dy->Integral();
  //cout << normFac << endl;

  //float normFac = 1.065;
  //float normFac = 1.09569;
  //hist_mass_pp->Scale(normFac);
  //hist_mass_pf->Scale(normFac);
  //hist_mass_ff->Scale(normFac);

  //hist_mass_dy->Scale(normFac);

  hist_mass_bkg_stack = new THStack("hist_mass_bkg_stack","Background");
  if (breakdown) {
    hist_mass_bkg_stack->Add(hist_mass_box);
    hist_mass_bkg_stack->Add(hist_mass_born);
    hist_mass_bkg_stack->Add(hist_mass_gjet_pf);
    hist_mass_bkg_stack->Add(hist_mass_qcd_pf);
    hist_mass_bkg_stack->Add(hist_mass_qcd_ff);
  } else {
    hist_mass_bkg_stack->Add(hist_mass_pp);
    hist_mass_bkg_stack->Add(hist_mass_pf);
    hist_mass_bkg_stack->Add(hist_mass_ff);
  }
  hist_mass_bkg_stack->Add(hist_mass_dy);

  hist_mass_bkg = (TH1*)hist_mass_box->Clone();
  hist_mass_bkg->Add(hist_mass_born);
  hist_mass_bkg->Add(hist_mass_gjet_pf);
  hist_mass_bkg->Add(hist_mass_qcd_pf);
  hist_mass_bkg->Add(hist_mass_qcd_ff);
  hist_mass_bkg->Add(hist_mass_dy);
  hist_mass_sig->Scale(5.);

  TCanvas *c_mgg = new TCanvas("c_mgg","Mgg",1000,700);
  c_mgg->SetFillColor(0);

  float max = hist_mass_bkg_stack->GetMaximum();
  if (hist_mass_sig->GetMaximum()>max) max=hist_mass_sig->GetMaximum();
  if (hist_mass_data->GetMaximum()>max) max=hist_mass_data->GetMaximum();
  if (var=="all_mass") {
    //hist_mass_bkg_stack->SetMaximum(1600.);
    //hist_mass_bkg_stack->SetMaximum(max*0.85);
    //hist_mass_bkg_stack->SetMaximum(max*0.65);
    //hist_mass_bkg_stack->SetMaximum(max*1.2);
  } else if (var=="pho1_eta" || var=="pho2_eta") {
    hist_mass_bkg_stack->SetMaximum(max*1.85);
  } else if (var=="pho1_pt" || var=="pho2_pt") {
    hist_mass_bkg_stack->SetMaximum(max*1.3);
  } else if (var=="pho1_r9" || var=="pho2_r9") {
    hist_mass_bkg_stack->SetMaximum(max*1.15);
  }

  hist_mass_bkg_stack->Draw("hist");
  //if (var=="all_mass") hist_mass_bkg_stack->GetXaxis()->SetRangeUser(100.,179.);
  if (var=="all_mass") hist_mass_bkg_stack->GetXaxis()->SetRangeUser(90.,190.);
  //hist_mass_bkg_stack->GetXaxis()->SetRangeUser(80.,120.);
  hist_mass_bkg_stack_err = (TH1*)(hist_mass_bkg_stack->GetStack()->Last())->Clone();
  for (int ibin=0; ibin<hist_mass_dy->GetNbinsX(); ibin++) {
    float err_pp = 0.15*hist_mass_pp->GetBinContent(ibin+1);
    float err_pf = 0.20*hist_mass_pf->GetBinContent(ibin+1);
    float err_ff = 0.50*hist_mass_ff->GetBinContent(ibin+1);
    float err_dy = 0.04*hist_mass_dy->GetBinContent(ibin+1);
    float err_tot = sqrt(err_pp*err_pp + err_pf*err_pf + err_ff*err_ff + err_dy*err_dy);
    hist_mass_bkg_stack_err->SetBinError(ibin+1,err_tot);
  }
  hist_mass_bkg_stack_err->SetFillStyle(3004);
  hist_mass_bkg_stack_err->SetFillColor(1);
  hist_mass_bkg_stack_err->Draw("same,e2");
  //hist_mass_bkg_stack->GetYaxis()->SetTitleSize(0.05);
  hist_mass_bkg_stack->GetYaxis()->SetTitle("Events / 2 GeV");
  //hist_mass_bkg_stack->GetXaxis()->SetTitleOffset(0.9);
  //hist_mass_bkg_stack->GetXaxis()->SetTitleSize(0.05);
  if (var=="all_mass") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  } else if (var=="pho1_pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon p_{T} (GeV)");
  } else if (var=="pho2_pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon p_{T} (GeV)");
  } else if (var=="pho1_r9") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon R_{9} (GeV)");
  } else if (var=="pho2_r9") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon R_{9} (GeV)");
  } else if (var=="pho1_eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon #eta (GeV)");
  } else if (var=="pho2_eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("sublead photon #eta (GeV)");
  } else if (var=="eta") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("di-photon #eta (GeV)");
  } else if (var=="pt") {
    hist_mass_bkg_stack->GetXaxis()->SetTitle("di-photon p_{T} (GeV)");
    hist_mass_bkg_stack->GetXaxis()->SetRangeUser(0.,120.);
  } else {
    cout << "unknown variable: " << var << endl;
  }
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon #pt");
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon p_{T}");
  //hist_mass_bkg_stack->GetXaxis()->SetTitle("lead photon R_{9}");
  hist_mass_sig->Draw("hist,same");
  //if (var!="all_mass") hist_mass_data->Draw("same,e");
  if (var=="all_mass"&& blind) {
    for (int ibin=16; ibin<36; ibin++) {
      hist_mass_data->SetBinContent(ibin,-999.);
    }
  }
  hist_mass_data->Draw("same,e");

  TLegend *leg_mass;
  if (var=="pho1_r9" || var=="pho2_r9") {
    leg_mass = new TLegend(.3,.48,.72,.78);
  }else {
    leg_mass = new TLegend(.6,.6,.92,.9);
  }
  leg_mass->SetBorderSize(0);
  leg_mass->SetFillColor(10);
  leg_mass->SetTextSize(.03);
  if (breakdown) {
    leg_mass->AddEntry(hist_mass_data,"Data","P");
    leg_mass->AddEntry(hist_mass_sig,"H#rightarrow#gamma#gamma (125 GeV) #times5","F");
    leg_mass->AddEntry(hist_mass_dy,"DYee+Z","F");
    leg_mass->AddEntry(hist_mass_qcd_ff,"QCD fake-fake","F");
    leg_mass->AddEntry(hist_mass_qcd_pf,"QCD prompt-fake","F");
    leg_mass->AddEntry(hist_mass_gjet_pf,"GJet prompt-fake","F");
    leg_mass->AddEntry(hist_mass_born,"Madgraph DiPhotonJets","F");
    leg_mass->AddEntry(hist_mass_box,"Pythia Box","F");
  } else {
    leg_mass->AddEntry(hist_mass_data,"Data","P");
    leg_mass->AddEntry(hist_mass_pp,"2 prompt #gamma","F");
    leg_mass->AddEntry(hist_mass_pf,"1 prompt #gamma 1 fake #gamma","F");
    leg_mass->AddEntry(hist_mass_ff,"2 fake #gamma","F");
    leg_mass->AddEntry(hist_mass_dy,"Drell-Yan","F");
    leg_mass->AddEntry(hist_mass_sig,"H#rightarrow#gamma#gamma (125 GeV) #times5","F");
  }
  if (var!="eta") leg_mass->Draw();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.05);
  txt->SetTextAlign(12);
  txt->DrawLatex(0.31,0.85,"#scale[0.8]{#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 12.2 fb^{-1}}}");
  //txt->DrawLatex(0.26,0.82,"CMS preliminary");
  //txt->DrawLatex(0.26,0.75,"#sqrt{s} = 7 TeV L = 1.66 fb^{-1}");

  //c_mgg->SaveAs("mass.gif");
  //c_mgg->SaveAs("mass.png");

  TString outcat="";
  TString outvar=var;
  if (category=="1") outcat = "_cat0";
  if (category=="2") outcat = "_cat1";
  if (category=="3") outcat = "_cat2";
  if (category=="4") outcat = "_cat3";
  if (category=="5") outcat = "_cat4";
  if (category=="6") outcat = "_cat5";
  if (category=="7") outcat = "_cat6";
  if (category=="8") outcat = "_cat7";
  if (category=="9") outcat = "_cat8";
  if (var=="all_mass") outvar = "mass";
  //c_mgg->SaveAs("/afs/cern.ch/user/f/futyand/www/hgg/topup/dataMC/massfacmva/"+outvar+"/"+outvar+outcat+".png");
  c_mgg->SaveAs(outvar+outcat+".png");

}
