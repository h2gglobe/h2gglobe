#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

int gp_n_;
short int gp_status_[10000];
short int gp_pdgid_[10000];
short int gp_mother_[10000];
TClonesArray *gp_p4_ = new TClonesArray("TLorentzVector",10000);
int pho_n_;
TClonesArray *pho_p4_ = new TClonesArray("TLorentzVector",10000);
TClonesArray *gh_pho1 = new TClonesArray("TLorentzVector",10000);
TClonesArray *gh_pho2 = new TClonesArray("TLorentzVector",10000);


double getCosTheta(TLorentzVector b1, TLorentzVector b2, TLorentzVector g1, TLorentzVector g2){
    
    TLorentzVector diphoton=g1+g2;
    TVector3 boostToDiphotonFrame = -diphoton.BoostVector();

    // Boost to higgs frame
    TLorentzVector refDIPHO_g1 = g1; refDIPHO_g1.Boost(boostToDiphotonFrame);
    TLorentzVector refDIPHO_b1 = b1; refDIPHO_b1.Boost(boostToDiphotonFrame);
    TLorentzVector refDIPHO_b2 = b2; refDIPHO_b2.Boost(boostToDiphotonFrame);

    // Getting beam 3-vector from 4-vectors
    TVector3 refDIPHO_vb1_direction = refDIPHO_b1.Vect().Unit();
    TVector3 refDIPHO_vb2_direction = refDIPHO_b2.Vect().Unit();

    // Definition of zz directions
    TVector3 direction_cs = (refDIPHO_vb1_direction - refDIPHO_vb2_direction).Unit(); // CS direction

    return TMath::Cos(direction_cs.Angle(refDIPHO_g1.Vect()));
}

map<string,vector<string> > getListOfFiles(string datfile){
  
  map<string,vector<string> > fileMap;
  ifstream file;
  file.open(datfile.c_str());
  if (file.fail()) exit(1);
  while (file.good()){
    string line;
    getline(file,line);
    if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
    string process = line.substr(0,line.find_first_of(" "));
    string filename = line.substr(line.find_first_of(" ")+1,string::npos);
    if (fileMap.find(process)==fileMap.end()){
      vector<string> temp;
      temp.push_back(filename);
      fileMap.insert(pair<string,vector<string> >(process,temp));
    }
    else {
      fileMap[process].push_back(filename);
    }
  }
  file.close();
  return fileMap;

}

void makePlot(string datfile, string outfile, double sqrtS=8.){
 
 	double beamE = 500.*sqrtS;
  TLorentzVector b1(0.,0.,beamE,beamE);
  TLorentzVector b2(0.,0.,-beamE,beamE);
	
	cout << "Using beam energy " << beamE << endl;

  TFile *outf = new TFile(outfile.c_str(),"RECREATE");
  int type=-1;
  int category=-1;
  float cosTheta=-999.;
  float cosThetaCS=-999.;

  TTree *outTree = new TTree("CosThetaTree","CosThetaTree");
  outTree->Branch("type",&type,"type/I");
  outTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
  outTree->Branch("cosThetaCS",&cosThetaCS,"cosThetaCS/F");
  outTree->Branch("category",&category,"category/I");
  

  map<string,vector<string> > fileMap = getListOfFiles(datfile);
  int procCount=0;
  for (map<string,vector<string> >::iterator mapIt=fileMap.begin(); mapIt!=fileMap.end(); mapIt++){
    string name = mapIt->first.substr(0,mapIt->first.find_first_of(":"));
    string temp = mapIt->first.substr(mapIt->first.find_first_of(":")+1,string::npos);
    string title = temp.substr(0,temp.find_first_of(":"));
    int color = atoi(temp.substr(temp.find_first_of(":")+1,string::npos).c_str());
    type=color;
    int fileCount=0;
    cout << "Running process " << name << " " << title << " " << color << endl;
    for (vector<string>::iterator file=mapIt->second.begin(); file!=mapIt->second.end(); file++){
     
      cout << "\t File: " << fileCount << "/" << mapIt->second.size() << endl;
      cout << "\t      - " << *file << endl;
      TFile *inFile = TFile::Open(file->c_str());
      TTree *tree = (TTree*)inFile->Get("event");

      tree->SetBranchAddress("gp_n",&gp_n_);
      tree->SetBranchAddress("gp_status",&gp_status_);
      tree->SetBranchAddress("gp_pdgid",&gp_pdgid_);
      tree->SetBranchAddress("gp_mother",&gp_mother_);
      tree->SetBranchAddress("gp_p4",&gp_p4_);

			//tree->SetBranchAddress("gh_pho1_p4",&gh_pho1);
			//tree->SetBranchAddress("gh_pho2_p4",&gh_pho2);
			tree->SetBranchAddress("pho_n",&pho_n_);
			tree->SetBranchAddress("pho_p4",&pho_p4_);
     
      for (int ev=0; ev<tree->GetEntries(); ev++){
        if (ev%100==0) cout << ev << "/" << tree->GetEntries() << "\r" << flush;
        tree->GetEntry(ev);

				/*
				//cout << "pho_n " << pho_n_ << endl;
				if (pho_n_<2) continue;

				int lead_ind=-1;
				int sublead_ind=-1;
				double maxpt=0.;
				double submaxpt=0.;
				for (int pi=0; pi<pho_n_; pi++){
					TLorentzVector *pho_p4 = (TLorentzVector*)(*pho_p4_)[pi];
					//cout << "\t pi: " << pi << " pt: " << pho_p4->Pt() << endl;
					if (pho_p4->Pt()>=maxpt) {
						// switch lead to sublead
						sublead_ind = lead_ind;
						lead_ind = pi;
						submaxpt = maxpt;
						maxpt = pho_p4->Pt();
					}
					else {
						if (pho_p4->Pt()>submaxpt && pho_p4->Pt()<maxpt){
							sublead_ind=pi;
							submaxpt = pho_p4->Pt();
						}
					}
				}
				assert(lead_ind!=-1);
				assert(sublead_ind!=-1);
				assert(lead_ind!=sublead_ind);
				//cout << "Chosen: " << lead_ind << " " << sublead_ind << endl;
				TLorentzVector *lead_p4 = (TLorentzVector*)(*pho_p4_)[lead_ind];
				TLorentzVector *sublead_p4 = (TLorentzVector*)(*pho_p4_)[sublead_ind];
       
        cosThetaCS = getCosTheta(b1,b2,*lead_p4,*sublead_p4);
				//cout << cosThetaCS << endl;
				*/
        vector<int> protonInd;
        vector<int> partonInd;
        vector<int> photonInd;
        int higgsIndex;
        
        for (int gi=0; gi<gp_n_; gi++){
          if (gp_status_[gi]!=3) continue;
          if (gp_pdgid_[gi]==2212) protonInd.push_back(gi);
          if (gp_pdgid_[gi]==25 || gp_pdgid_[gi]==39) higgsIndex = gi;
        }
        assert(protonInd.size()==2);
        for (int gi=0; gi<gp_n_; gi++){
          if (gp_status_[gi]!=3) continue;
          if (gp_mother_[gi]==higgsIndex) photonInd.push_back(gi);
          if (gp_mother_[gi]==protonInd[0] || gp_mother_[gi]==protonInd[1]) partonInd.push_back(gi); 
        }

        //cout << "protons at.." << protonInd[0] << "," << protonInd[1] << endl;
        //cout << "partons at.." << partonInd[0] << "," << partonInd[1] << endl;
        //cout << "photons at.." << photonInd[0] << "," << photonInd[1] << endl;
        //cout << "higgs   at.." << higgsIndex << endl;
        TLorentzVector *gen_higgs = (TLorentzVector*)(*gp_p4_)[higgsIndex];
        TLorentzVector *gen_part1 = (TLorentzVector*)(*gp_p4_)[partonInd[0]];
        TLorentzVector *gen_part2 = (TLorentzVector*)(*gp_p4_)[partonInd[1]];
        TLorentzVector *gen_pho1 = (TLorentzVector*)(*gp_p4_)[photonInd[0]];
        TLorentzVector *gen_pho2 = (TLorentzVector*)(*gp_p4_)[photonInd[1]];
        //cout << gen_higgs->Px() << " " << gen_higgs->Py() << " " << gen_higgs->Pz() << " " << gen_higgs->E() << endl;
        //category = getCategory(lead_eta,sublead_eta,lead_r9,sublead_r9);

        cosThetaCS = getCosTheta(b1,b2,*gen_pho1,*gen_pho2);
				//cout << cosThetaCS << endl;
        //cosTheta   = getCosTheta(*gen_part1,*gen_part2,*gen_pho1,*gen_pho2);
        
        outTree->Fill();
      } // loop events
      cout << endl;
      fileCount++;
    } // loop files
    procCount++;
  } // loop procs
  outf->cd();
  outTree->Write();
  outf->Close();

}
