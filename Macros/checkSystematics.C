#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

using namespace std;

void checkSystematics(string filename){

  gROOT->SetBatch();

  TFile *tFile = TFile::Open(filename.c_str());
   
  const int ncats=6;

  vector<string> processes;
  processes.push_back("ggh");
  processes.push_back("vbf");
  processes.push_back("wzh");
  processes.push_back("tth");

  vector<string> systematics;
  systematics.push_back("E_res");
  systematics.push_back("E_scale");
  systematics.push_back("idEff");
  systematics.push_back("pdfWeight");
  systematics.push_back("phoIdMva");
  systematics.push_back("regSig");
  systematics.push_back("triggerEff");
  systematics.push_back("vtxEff");
  
  system("mkdir systematicPlots");

  for (int m=110; m<=150; m++){
    for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++) {
      for (int cat=0; cat<ncats; cat++){

        TH1F *nomSig = (TH1F*)tFile->Get(Form("th1f_sig_%s_mass_m%3d_cat%d",proc->c_str(),m,cat)); 
        TH1F *sigRV = (TH1F*)tFile->Get(Form("th1f_sig_%s_mass_m%3d_rv_cat%d",proc->c_str(),m,cat)); 
        TH1F *sigWV = (TH1F*)tFile->Get(Form("th1f_sig_%s_mass_m%3d_wv_cat%d",proc->c_str(),m,cat)); 
        TCanvas *canv = new TCanvas();
        TLegend *leg = new TLegend(0.7,0.7,0.89,0.89);
        nomSig->SetLineWidth(2);
        nomSig->SetLineColor(kBlack);
        sigRV->SetLineColor(kBlue);
        sigRV->SetFillColor(kBlue);
        sigRV->SetFillStyle(3003);
        sigWV->SetLineColor(kRed);
        sigWV->SetFillColor(kRed);
        sigWV->SetFillStyle(3004);
        leg->AddEntry(nomSig,"Nominal Signal","f");
        leg->AddEntry(sigRV,"Right vertex","f");
        leg->AddEntry(sigWV,"Wrong vertex","f");
        nomSig->GetYaxis()->SetRangeUser(0.,1.3*TMath::Max(nomSig->GetMaximum(),TMath::Max(sigRV->GetMaximum(),sigWV->GetMaximum())));
        nomSig->Draw();
        sigRV->Draw("same");
        sigWV->Draw("same");
        leg->Draw("same");
        canv->Print(Form("systematicPlots/rvwv_%s_m%d_cat%d.pdf",proc->c_str(),m,cat));
        canv->Print(Form("systematicPlots/rvwv_%s_m%d_cat%d.png",proc->c_str(),m,cat));
        delete leg;

        for (vector<string>::iterator syst=systematics.begin(); syst!=systematics.end(); syst++){
          
          TLegend *systLeg = new TLegend(0.7,0.7,0.89,0.89);
          TH1F *up = (TH1F*)tFile->Get(Form("th1f_sig_%s_mass_m%3d_cat%d_%sUp01_sigma",proc->c_str(),m,cat,syst->c_str()));
          TH1F *down = (TH1F*)tFile->Get(Form("th1f_sig_%s_mass_m%3d_cat%d_%sDown01_sigma",proc->c_str(),m,cat,syst->c_str()));
          up->SetLineColor(kRed);
          down->SetLineColor(kBlue);
          systLeg->AddEntry(nomSig,"Nominal Signal","f");
          systLeg->AddEntry(up,Form("%s Up",syst->c_str()),"f");
          systLeg->AddEntry(down,Form("%s Down",syst->c_str()),"f");
          nomSig->GetYaxis()->SetRangeUser(0.,1.3*TMath::Max(nomSig->GetMaximum(),TMath::Max(up->GetMaximum(),down->GetMaximum())));
          nomSig->Draw();
          up->Draw("same");
          down->Draw("same");
          systLeg->Draw("same");
          canv->Print(Form("systematicPlots/%s_%s_m%d_cat%d.pdf",syst->c_str(),proc->c_str(),m,cat));
          canv->Print(Form("systematicPlots/%s_%s_m%d_cat%d.png",syst->c_str(),proc->c_str(),m,cat));
          delete systLeg;
        
        }
        delete canv;
      }
    }
  }

  system("python make_syst_html.py systematicPlots/");

  cout << "Done. Please cp -r systematicPlots to some public web space to view files" << endl;

}
