#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

using namespace std;

double getCSangle(TLorentzVector p1, TLorentzVector p2){
  
  TLorentzVector higgs=p1+p2;
  return 2.*(p2.E()*p1.Pz()-p1.E()*p2.Pz())/(higgs.M()*TMath::Sqrt(higgs.M()*higgs.M()+higgs.Pt()*higgs.Pt()));

}

void fillDists(vector<TTree*> trees, TH1F *nom, TH1F *up, TH1F *down, double pt_shift){
 
  float evweight;
  double lead_E;
  double lead_px;
  double lead_py;
  double lead_pz;
  double sublead_E;
  double sublead_px;
  double sublead_py;
  double sublead_pz;

  for (vector<TTree*>::iterator it=trees.begin(); it!=trees.end(); it++){
    (*it)->SetBranchAddress("evweight",&evweight);
    (*it)->SetBranchAddress("lead_E",&lead_E);
    (*it)->SetBranchAddress("lead_px",&lead_px);
    (*it)->SetBranchAddress("lead_py",&lead_py);
    (*it)->SetBranchAddress("lead_pz",&lead_pz);
    (*it)->SetBranchAddress("sublead_E",&sublead_E);
    (*it)->SetBranchAddress("sublead_px",&sublead_px);
    (*it)->SetBranchAddress("sublead_py",&sublead_py);
    (*it)->SetBranchAddress("sublead_pz",&sublead_pz);

    for (int ev=0; ev<(*it)->GetEntries(); ev++){
      (*it)->GetEntry(ev);
      TLorentzVector p1(lead_px,lead_py,lead_pz,lead_E);
      TLorentzVector p2(sublead_px,sublead_py,sublead_pz,sublead_E);
      TLorentzVector higgs=p1+p2;
      //nom->Fill(higgs.Pt(),evweight);
      //up->Fill(higgs.Pt()*(1.+pt_shift),evweight);
      //down->Fill(higgs.Pt()*(1.-pt_shift),evweight);
      nom->Fill(higgs.Pt());
      up->Fill(higgs.Pt()*(1.+pt_shift));
      down->Fill(higgs.Pt()*(1.-pt_shift));
    }
  }


}

void reweightPt(string infile, int sqrts=8, double pt_shift=0.1){

  TFile *inFile = TFile::Open(infile.c_str());

  TFile *outFile = new TFile(Form("ptweights_%dTeV.root",sqrts),"RECREATE");
  vector<int> masses;
  masses.push_back(125);
	masses.push_back(126);

  TCanvas *canv = new TCanvas();
  for (vector<int>::iterator mass=masses.begin(); mass!=masses.end(); mass++){
    
		int fmass=*mass;
		if (sqrts==7) fmass=125;
		TTree *gghTree = (TTree*)inFile->Get(Form("spin_trees/ggh_m%d_%dTeV",fmass,sqrts));
    TTree *vbfTree = (TTree*)inFile->Get(Form("spin_trees/vbf_m%d_%dTeV",fmass,sqrts));
    TTree *wzhTree = (TTree*)inFile->Get(Form("spin_trees/wzh_m%d_%dTeV",fmass,sqrts));
    TTree *tthTree = (TTree*)inFile->Get(Form("spin_trees/tth_m%d_%dTeV",fmass,sqrts));
    TTree *ggh_gravTree = (TTree*)inFile->Get(Form("spin_trees/gg_grav_m%d_%dTeV",fmass,sqrts));
    TTree *vbf_gravTree = (TTree*)inFile->Get(Form("spin_trees/qq_grav_m%d_%dTeV",fmass,sqrts));

    vector<TTree*> smTrees;
    smTrees.push_back(gghTree);
    smTrees.push_back(vbfTree);
    smTrees.push_back(wzhTree);
    smTrees.push_back(tthTree);
    
    vector<TTree*> gg_gravTrees;
    gg_gravTrees.push_back(ggh_gravTree);

    vector<TTree*> qq_gravTrees;
    qq_gravTrees.push_back(vbf_gravTree);


    TH1F *smNom = new TH1F(Form("smNom%d",*mass),"",20,0.,100.);
    TH1F *smUp = new TH1F(Form("smUp%d",*mass),"",20,0.,100.);
    TH1F *smDown = new TH1F(Form("smDown%d",*mass),"",20,0.,100.);
    TH1F *ggNom = new TH1F(Form("ggNom%d",*mass),"",20,0.,100.);
    TH1F *ggUp = new TH1F(Form("ggUp%d",*mass),"",20,0.,100.);
    TH1F *ggDown = new TH1F(Form("ggDown%d",*mass),"",20,0.,100.);
    TH1F *qqNom = new TH1F(Form("qqNom%d",*mass),"",20,0.,100.);
    TH1F *qqUp = new TH1F(Form("qqUp%d",*mass),"",20,0.,100.);
    TH1F *qqDown = new TH1F(Form("qqDown%d",*mass),"",20,0.,100.);

    smNom->Sumw2();
    smUp->Sumw2();
    smDown->Sumw2();
    ggNom->Sumw2();
    ggUp->Sumw2();
    ggDown->Sumw2();
    qqNom->Sumw2();
    qqUp->Sumw2();
    qqDown->Sumw2();

    fillDists(smTrees,smNom,smUp,smDown,pt_shift);
    fillDists(gg_gravTrees,ggNom,ggUp,ggDown,pt_shift);
    fillDists(qq_gravTrees,qqNom,qqUp,qqDown,pt_shift);

    TH1F *smNomRat = (TH1F*)smNom->Clone(Form("smNomRat%d",*mass));
    smNomRat->Divide(smNom);
    TH1F *smUpRat = (TH1F*)smUp->Clone(Form("smUpRat%d",*mass));
    smUpRat->Divide(smNom);
    TH1F *smDownRat = (TH1F*)smDown->Clone(Form("smDownRat%d",*mass));
    smDownRat->Divide(smNom);

    TH1F *ggNomRat = (TH1F*)ggNom->Clone(Form("ggNomRat%d",*mass));
    ggNomRat->Divide(ggNom);
    TH1F *ggUpRat = (TH1F*)ggUp->Clone(Form("ggUpRat%d",*mass));
    ggUpRat->Divide(ggNom);
    TH1F *ggDownRat = (TH1F*)ggDown->Clone(Form("ggDownRat%d",*mass));
    ggDownRat->Divide(ggNom);
    
    TH1F *qqNomRat = (TH1F*)qqNom->Clone(Form("qqNomRat%d",*mass));
    qqNomRat->Divide(qqNom);
    TH1F *qqUpRat = (TH1F*)qqUp->Clone(Form("qqUpRat%d",*mass));
    qqUpRat->Divide(qqNom);
    TH1F *qqDownRat = (TH1F*)qqDown->Clone(Form("qqDownRat%d",*mass));
    qqDownRat->Divide(qqNom);

    smNom->Write();
    smUp->Write();
    smDown->Write();
    ggNom->Write();
    ggUp->Write();
    ggDown->Write();
    qqNom->Write();
    qqUp->Write();
    qqDown->Write();

    smNomRat->Write();
    smUpRat->Write();
    smDownRat->Write();
    ggNomRat->Write();
    ggUpRat->Write();
    ggDownRat->Write();
    qqNomRat->Write();
    qqUpRat->Write();
    qqDownRat->Write();

    smNomRat->SetLineColor(1);
    smUpRat->SetLineColor(4);
    smDownRat->SetLineColor(2);
    smNomRat->Draw("HIST");
    smUpRat->Draw("HISTsame");
    smDownRat->Draw("HISTsame");
    canv->Print(Form("ptweights_sm_m%d.pdf",*mass));
    ggNomRat->SetLineColor(1);
    ggUpRat->SetLineColor(4);
    ggDownRat->SetLineColor(2);
    ggNomRat->Draw("HIST");
    ggUpRat->Draw("HISTsame");
    ggDownRat->Draw("HISTsame");
    canv->Print(Form("ptweights_gg_m%d.pdf",*mass));
    qqNomRat->SetLineColor(1);
    qqUpRat->SetLineColor(4);
    qqDownRat->SetLineColor(2);
    qqNomRat->Draw("HIST");
    qqUpRat->Draw("HISTsame");
    qqDownRat->Draw("HISTsame");
    canv->Print(Form("ptweights_qq_m%d.pdf",*mass));
  }

  outFile->Close();
}
