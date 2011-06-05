
void compare(){

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  TFile * fin = new TFile("CMS-HGG_105_140_0_5.root","READONLY");

  const int nmassMC = 4;
  Float_t massMC[nmassMC] = {110, 115, 120, 130} ;

  const int ncat = 8;
  
  //   for (int i=0; i<nmassMC; i++) {
    
  //     c = new TCanvas("c","c",1200,600);
  //     c->Divide(4,2);

  //     for (int j=0; j<ncat; j++) {

  //       c->cd(j+1);

  //       char hname[100];
  //       sprintf(hname,"th1f_sig_mass_m%i_cat%i",massMC[i],j);
  //       TH1F * hist = (TH1F*)fin->Get(hname);
  //       sprintf(hname,"th1f_sig_mass_m%i_cat%i_MC",massMC[i],j);
  //       TH1F * histActual = (TH1F*)fin->Get(hname);

  //       hist->SetLineColor(kBlue);
  //       hist->SetLineWidth(2);

  //       histActual->SetLineColor(2);
  //       histActual->SetLineStyle(kDashed);
  //       histActual->SetLineWidth(2);

  //       sprintf(hname,"Interpolated Mass at %iGeV cat%i",massMC[i],j);
  //       hist->SetTitle(hname);

  //       hist->SetMaximum(1.1*TMath::Max(hist->GetMaximum(),histActual->GetMaximum()));
  //       hist->Draw();
  //       histActual->Draw("same");

  //       TLegend *leg = new TLegend(0.6522989,0.6800847,0.8534483,0.8580508,NULL,"brNDC");
  //       leg->SetBorderSize(1);
  //       leg->SetTextFont(62);
  //       leg->SetTextSize(0.03225806);
  //       leg->SetLineColor(0);
  //       leg->SetFillColor(0);
  //       leg->AddEntry(histActual,"expected");
  //       leg->AddEntry(hist,"interpolated");
  //       leg->Draw();
      
  //     } 
  //     sprintf(hname,"interpolation_m%i.png",massMC[i]);
  //     c->Print(hname);   
  //   }  

  const int nmass = 70;
  int massBegin = 105;

  float yields[ncat][nmass];
  float masses[ncat][nmass];

  for (int i=0; i<nmass; i++) {

    for (int j=0; j<ncat; j++) {
      
      float mass = massBegin+i*0.5;
      
      char hname[100];
      TH1F* hist;
      if (i%2) {
        sprintf(hname,"th1f_sig_mass_m%i_5_cat%i",mass-0.5,j);
        hist = (TH1F*)fin->Get(hname);
        cout << " hello1 " << i << " " << j << " " << hname << " " << hist << " " << hist->Integral() << endl;
      } else {
        sprintf(hname,"th1f_sig_mass_m%i_cat%i",mass,j);
        hist = (TH1F*)fin->Get(hname);
        cout << " hello2 " << i << " " << j << " " << hname << " " << hist << " " << hist->Integral() << endl;
      }

      yields[j][i] = hist->Integral();
      masses[j][i] = mass;

      cout << i << " " << mass << " " << yields[j][i] << endl;

    }

  }

  float mcyields[ncat][nmassMC];
  
  for (int i=0; i<nmassMC; i++) {

    for (int j=0; j<ncat; j++) {

      char hname[100];
      TH1F* hist;
      sprintf(hname,"th1f_sig_mass_m%i_cat%i_MC",massMC[i],j);
      hist = (TH1F*)fin->Get(hname);
      cout << " hello3 " << i << " " << j << " " << hname << " " << hist << " " << hist->Integral() << endl;
      mcyields[j][i] = hist->Integral();

      cout << i << " " << massMC[i] << mcyields[j][i] << endl;
    }

  }

  TGraph* gr[ncat];
  TGraph* grMC[ncat];
  
  c1 = new TCanvas("c1","c1",600,600);
  c1->cd();

  TLegend *leg = new TLegend(0.7332215,0.6730769,0.8909396,0.8776224,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03225806);
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  for (int j=0; j<ncat; j++) {

    gr[j] = new TGraph(nmass,masses[j],yields[j]);
    gr[j]->SetMarkerStyle(20);
    gr[j]->SetMarkerSize(1);
    gr[j]->SetMarkerColor(2+j);

    grMC[j] = new TGraph(nmassMC,massMC,mcyields[j]);
    grMC[j]->SetMarkerStyle(20);
    grMC[j]->SetMarkerSize(1);
    grMC[j]->SetMarkerColor(1);
    
    if (j==0) {
      gr[j]->SetTitle("Interpolated Yields");
      gr[j]->SetMaximum(0.9);
      gr[j]->SetMinimum(0.0);
      gr[j]->Draw("ap");
      grMC[j]->Draw("p");
    } else {
      gr[j]->Draw("p");
      grMC[j]->Draw("p");
    }
    leg->AddEntry(gr[j],Form("cat%i",j));
  }

  leg->Draw();

  c1->Print("yields_v_mass_all.png");


}
