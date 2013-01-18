void BDT_Zee_bdtin_cat() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  
  gStyle->SetPalette(1);

  TFile *file = TFile::Open("histograms_CMS-HGG_zeevalidation.root");
  file->cd();

  txt = new TLatex();
  txt->SetNDC();
  txt->SetTextSize(0.08);

  mass_cat0_Data->GetXaxis()->SetTitleSize(0.05);
  mass_cat0_DYJetsToLL->SetFillColor(38);
  mass_cat0_DYJetsToLL->SetLineColor(1);
  mass_cat0_Data->SetMarkerStyle(20);
  mass_cat0_Data->SetMarkerSize(0.4);
  mass_cat0_Data->SetLineColor(1);

  TLegend *leg;
  leg = new TLegend(.6,.65,.87,.87);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->AddEntry(mass_cat0_Data,"Data (19.6fb^{-1})");
  leg->AddEntry(mass_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLegend *leg2;
  leg2 = new TLegend(.2,.65,.52,.87);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(10);
  leg2->SetTextSize(.035);
  leg2->AddEntry(mass_cat0_Data,"Data (19.6fb^{-1})");
  leg2->AddEntry(mass_cat0_DYJetsToLL,"DYJetsToLL MC","F");

  TLine *line[13];
  for (int i=0; i<13; i++) {
    line[i] = new TLine(-1.,1.,1.,1.);
    line[i]->SetLineColor(4);
    line[i]->SetLineWidth(2);
    //if (i==4) line[i]->SetX1(0.3);
    if (i==0 || i==1) {
      line[i]->SetX1(-0.4);
      line[i]->SetX2(0.4);
    }
    if (i==2 || i==3) {
      line[i]->SetX1(0.);
      line[i]->SetX2(0.06);
    }
    if (i==5) {
      line[i]->SetX1(0.2);
      line[i]->SetX2(1.2);
    }
    if (i==6) {
      line[i]->SetX1(0.2);
      line[i]->SetX2(0.9);
    }
    if (i==7 || i==8) {
      line[i]->SetX1(-2.5);
      line[i]->SetX2(2.5);
    }
    if (i==11) {
      line[i]->SetX1(0.);
      line[i]->SetX2(200.);
    }
  }

  TString var[13] = {"pho1_phoidMva","pho2_phoidMva","sigmaMOverM","sigmaMOverM_wrongVtx","vtxProb","pho1_ptOverM","pho2_ptOverM","pho1_eta","pho2_eta","cosDeltaPhi","mass","pt","nvtx"};
  TString var_name[13] = {"lead photon ID MVA output","sublead photon ID MVA output","#sigma_{M}/M_{#gamma#gamma} (right vertex)","#sigma_{M}/M_{#gamma#gamma} (wrong vertex)","Vertex probabilty","lead p_{T}/M_{#gamma#gamma}","sublead p_{T}/M_{#gamma#gamma}","lead #eta","sublead #eta","cos(#Delta#phi)","m_{e^{+}e^{-}}","di-photon p_{T}","no. of primary vertices"};

  TCanvas *c_invar[13];
  TH1* hist_cat1_Data[13];
  TH1* hist_cat2_Data[13];
  TH1* hist_cat3_Data[13];
  TH1* hist_cat4_Data[13];
  TH1* hist_cat5_Data[13];
  TH1* hist_cat6_Data[13];
  TH1* hist_cat1_DYJetsToLL[13];
  TH1* hist_cat2_DYJetsToLL[13];
  TH1* hist_cat3_DYJetsToLL[13];
  TH1* hist_cat4_DYJetsToLL[13];
  TH1* hist_cat5_DYJetsToLL[13];
  TH1* hist_cat6_DYJetsToLL[13];

  for (int i=0; i<13; i++) {

    bool logy = false;
    if (i==4 || i==9) logy=true;

    float xmin,xmax;
    bool setRange=true;
    if (i==0 || i==1) {
      xmin=-0.4;
      xmax=0.4;
    } else if (i==4) {
      xmin=0.;
      xmax=1.;
    } else if (i==5) {
      xmin=0.2;
      xmax=1.2;
    } else if (i==6) {
      xmin=0.2;
      xmax=0.9;
    } else {
      setRange=false;
    }

    hist_cat1_Data[i] = (TH1*)(file->Get(var[i]+"_cat1_Data"))->Clone();
    hist_cat2_Data[i] = (TH1*)(file->Get(var[i]+"_cat2_Data"))->Clone();
    hist_cat3_Data[i] = (TH1*)(file->Get(var[i]+"_cat3_Data"))->Clone();
    hist_cat4_Data[i] = (TH1*)(file->Get(var[i]+"_cat4_Data"))->Clone();
    hist_cat5_Data[i] = (TH1*)(file->Get(var[i]+"_cat5_Data"))->Clone();
    hist_cat6_Data[i] = (TH1*)(file->Get(var[i]+"_cat6_Data"))->Clone();
    hist_cat1_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat1_DYJetsToLL"))->Clone();
    hist_cat2_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat2_DYJetsToLL"))->Clone();
    hist_cat3_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat3_DYJetsToLL"))->Clone();
    hist_cat4_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat4_DYJetsToLL"))->Clone();
    hist_cat5_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat5_DYJetsToLL"))->Clone();
    hist_cat6_DYJetsToLL[i] = (TH1*)(file->Get(var[i]+"_cat6_DYJetsToLL"))->Clone();

    hist_cat1_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (category 0)");
    hist_cat2_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (category 1)");
    hist_cat3_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (category 2)");
    hist_cat4_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (category 3)");
    hist_cat5_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (dijet cat1)");
    hist_cat6_Data[i]->GetXaxis()->SetTitle(var_name[i]+" (dijet cat2)");
    hist_cat1_Data[i]->GetYaxis()->SetTitle("");
    hist_cat2_Data[i]->GetYaxis()->SetTitle("");
    hist_cat3_Data[i]->GetYaxis()->SetTitle("");
    hist_cat4_Data[i]->GetYaxis()->SetTitle("");
    hist_cat5_Data[i]->GetYaxis()->SetTitle("");
    hist_cat6_Data[i]->GetYaxis()->SetTitle("");
    if (setRange) {
      hist_cat2_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      hist_cat3_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      hist_cat4_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      hist_cat5_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      hist_cat6_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    }
    hist_cat1_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat2_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat3_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat4_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat5_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat6_Data[i]->GetXaxis()->SetTitleSize(0.05);
    hist_cat1_DYJetsToLL[i]->SetFillColor(38);
    hist_cat2_DYJetsToLL[i]->SetFillColor(38);
    hist_cat3_DYJetsToLL[i]->SetFillColor(38);
    hist_cat4_DYJetsToLL[i]->SetFillColor(38);
    hist_cat5_DYJetsToLL[i]->SetFillColor(38);
    hist_cat6_DYJetsToLL[i]->SetFillColor(38);
    hist_cat1_DYJetsToLL[i]->SetLineColor(1);
    hist_cat2_DYJetsToLL[i]->SetLineColor(1);
    hist_cat3_DYJetsToLL[i]->SetLineColor(1);
    hist_cat4_DYJetsToLL[i]->SetLineColor(1);
    hist_cat5_DYJetsToLL[i]->SetLineColor(1);
    hist_cat6_DYJetsToLL[i]->SetLineColor(1);
    hist_cat1_Data[i]->SetMarkerStyle(20);
    hist_cat2_Data[i]->SetMarkerStyle(20);
    hist_cat3_Data[i]->SetMarkerStyle(20);
    hist_cat4_Data[i]->SetMarkerStyle(20);
    hist_cat5_Data[i]->SetMarkerStyle(20);
    hist_cat6_Data[i]->SetMarkerStyle(20);
    hist_cat1_Data[i]->SetMarkerSize(0.4);
    hist_cat2_Data[i]->SetMarkerSize(0.4);
    hist_cat3_Data[i]->SetMarkerSize(0.4);
    hist_cat4_Data[i]->SetMarkerSize(0.4);
    hist_cat5_Data[i]->SetMarkerSize(0.4);
    hist_cat6_Data[i]->SetMarkerSize(0.4);
    hist_cat1_Data[i]->SetLineColor(1);
    hist_cat2_Data[i]->SetLineColor(1);
    hist_cat3_Data[i]->SetLineColor(1);
    hist_cat4_Data[i]->SetLineColor(1);
    hist_cat5_Data[i]->SetLineColor(1);
    hist_cat6_Data[i]->SetLineColor(1);

    c_invar[i] = new TCanvas("c_"+var[i],var[i]+" in di-photon BDT categories",1600,600);
    c_invar[i]->Divide(6,2);

    c_invar[i]->cd(1);
    if (logy) gPad->SetLogy();
    float sf = hist_cat1_Data[i]->Integral()/hist_cat1_DYJetsToLL[i]->Integral();
    hist_cat1_DYJetsToLL[i]->Scale(sf);
    //if (i==11 || i==9) {
      hist_cat1_Data[i]->Rebin(2);
      hist_cat1_DYJetsToLL[i]->Rebin(2);
      if (setRange) hist_cat1_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    //}
    hist_cat1_Data[i]->Draw("e");
    hist_cat1_DYJetsToLL[i]->Draw("hist,same");
    hist_cat1_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    if (i==0 || i==1 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8 && i!=9 && i!=5) {
      leg->Draw();
    }

    c_invar[i]->cd(2);
    if (logy) gPad->SetLogy();
    float sf = hist_cat2_Data[i]->Integral()/hist_cat2_DYJetsToLL[i]->Integral();
    hist_cat2_DYJetsToLL[i]->Scale(sf);
    hist_cat2_Data[i]->Draw("e");
    hist_cat2_DYJetsToLL[i]->Draw("hist,same");
    hist_cat2_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8) {
      leg->Draw();
    }

    c_invar[i]->cd(3);
    if (logy) gPad->SetLogy();
    float sf = hist_cat3_Data[i]->Integral()/hist_cat3_DYJetsToLL[i]->Integral();
    hist_cat3_DYJetsToLL[i]->Scale(sf);
    hist_cat3_Data[i]->Draw("e");
    hist_cat3_DYJetsToLL[i]->Draw("hist,same");
    hist_cat3_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8) {
      leg->Draw();
    }

    c_invar[i]->cd(4);
    if (logy) gPad->SetLogy();
    float sf = hist_cat4_Data[i]->Integral()/hist_cat4_DYJetsToLL[i]->Integral();
    hist_cat4_DYJetsToLL[i]->Scale(sf);
    hist_cat4_Data[i]->Draw("e");
    hist_cat4_DYJetsToLL[i]->Draw("hist,same");
    hist_cat4_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8 && i!=1) {
      leg->Draw();
    }

    c_invar[i]->cd(5);
    if (logy) gPad->SetLogy();
    float sf = hist_cat5_Data[i]->Integral()/hist_cat5_DYJetsToLL[i]->Integral();
    hist_cat5_DYJetsToLL[i]->Scale(sf);
    //if (i==11 || i==9) {
    if (i==7 || i==8) {
      hist_cat5_Data[i]->Rebin(5);
      hist_cat5_DYJetsToLL[i]->Rebin(5);
    } else {
      hist_cat5_Data[i]->Rebin(4);
      hist_cat5_DYJetsToLL[i]->Rebin(4);
    }
    if (setRange) hist_cat5_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    hist_cat5_Data[i]->Draw("e");
    hist_cat5_DYJetsToLL[i]->Draw("hist,same");
    hist_cat5_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8 && i!=1) {
      leg->Draw();
    }

    c_invar[i]->cd(6);
    if (logy) gPad->SetLogy();
    float sf = hist_cat6_Data[i]->Integral()/hist_cat6_DYJetsToLL[i]->Integral();
    hist_cat6_DYJetsToLL[i]->Scale(sf);
    //if (i==11 || i==9) {
    if (i==7 || i==8) {
      hist_cat6_Data[i]->Rebin(5);
      hist_cat6_DYJetsToLL[i]->Rebin(5);
    } else {
      hist_cat6_Data[i]->Rebin(4);
      hist_cat6_DYJetsToLL[i]->Rebin(4);
    }
    if (setRange) hist_cat6_Data[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    hist_cat6_Data[i]->Draw("e");
    hist_cat6_DYJetsToLL[i]->Draw("hist,same");
    hist_cat6_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 || i==4) {
      leg2->Draw();
    } else if (i!=7 && i!=8 && i!=1) {
      leg->Draw();
    }

    /*
    c_invar[i]->cd(6);
    gPad->SetLogy();
    hist_cat1_Data[i]->Draw("e");
    hist_cat1_DYJetsToLL[i]->Draw("hist,same");
    hist_cat1_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 ) {
      leg2->Draw();
    } else {
      leg->Draw();
    }

    c_invar[i]->cd(7);
    gPad->SetLogy();
    hist_cat2_Data[i]->Draw("e");
    hist_cat2_DYJetsToLL[i]->Draw("hist,same");
    hist_cat2_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 ) {
      leg2->Draw();
    } else {
      leg->Draw();
    }

    c_invar[i]->cd(8);
    gPad->SetLogy();
    hist_cat3_Data[i]->Draw("e");
    hist_cat3_DYJetsToLL[i]->Draw("hist,same");
    hist_cat3_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 ) {
      leg2->Draw();
    } else {
      leg->Draw();
    }

    c_invar[i]->cd(9);
    gPad->SetLogy();
    hist_cat4_Data[i]->Draw("e");
    hist_cat4_DYJetsToLL[i]->Draw("hist,same");
    hist_cat4_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 ) {
      leg2->Draw();
    } else {
      leg->Draw();
    }

    c_invar[i]->cd(10);
    gPad->SetLogy();
    hist_cat5_Data[i]->Draw("e");
    hist_cat5_DYJetsToLL[i]->Draw("hist,same");
    hist_cat5_Data[i]->Draw("e,same");
    gPad->RedrawAxis();
    
    if (i==0 || i==1 || i==3 ) {
      leg2->Draw();
    } else {
      leg->Draw();
    }
    */

    c_invar[i]->cd(7);
    gPad->SetGrid();
    ratio_cat1 = (TH1*)hist_cat1_Data[i]->Clone();
    ratio_cat1->Sumw2();
    ratio_cat1->Divide(hist_cat1_DYJetsToLL[i]);
    ratio_cat1->SetMaximum(1.8);
    ratio_cat1->SetMinimum(0.2);
    ratio_cat1->Draw("e");
    line[i]->Draw();

    c_invar[i]->cd(8);
    gPad->SetGrid();
    ratio_cat2 = (TH1*)hist_cat2_Data[i]->Clone();
    ratio_cat2->Sumw2();
    ratio_cat2->Divide(hist_cat2_DYJetsToLL[i]);
    ratio_cat2->SetMaximum(1.8);
    ratio_cat2->SetMinimum(0.2);
    ratio_cat2->Draw("e");
    line[i]->Draw();

    c_invar[i]->cd(9);
    gPad->SetGrid();
    ratio_cat3 = (TH1*)hist_cat3_Data[i]->Clone();
    ratio_cat3->Sumw2();
    ratio_cat3->Divide(hist_cat3_DYJetsToLL[i]);
    ratio_cat3->SetMaximum(1.8);
    ratio_cat3->SetMinimum(0.2);
    ratio_cat3->Draw("e");
    line[i]->Draw();

    c_invar[i]->cd(10);
    gPad->SetGrid();
    ratio_cat4 = (TH1*)hist_cat4_Data[i]->Clone();
    ratio_cat4->Sumw2();
    ratio_cat4->Divide(hist_cat4_DYJetsToLL[i]);
    ratio_cat4->SetMaximum(1.8);
    ratio_cat4->SetMinimum(0.2);
    ratio_cat4->Draw("e");
    line[i]->Draw();

    c_invar[i]->cd(11);
    gPad->SetGrid();
    ratio_cat5 = (TH1*)hist_cat5_Data[i]->Clone();
    ratio_cat5->Sumw2();
    ratio_cat5->Divide(hist_cat5_DYJetsToLL[i]);
    ratio_cat5->SetMaximum(1.8);
    ratio_cat5->SetMinimum(0.2);
    ratio_cat5->Draw("e");
    line[i]->Draw();

    c_invar[i]->cd(12);
    gPad->SetGrid();
    ratio_cat6 = (TH1*)hist_cat6_Data[i]->Clone();
    ratio_cat6->Sumw2();
    ratio_cat6->Divide(hist_cat6_DYJetsToLL[i]);
    ratio_cat6->SetMaximum(1.8);
    ratio_cat6->SetMinimum(0.2);
    ratio_cat6->Draw("e");
    line[i]->Draw();

    c_invar[i]->SaveAs(var[i]+"_cat.png");

  }

}
