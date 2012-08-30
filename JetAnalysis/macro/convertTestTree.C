#include "TROOT.h"

struct MultiBDT { 
	float Signal, bg0, bg1; 
};


void convertTestTree(TString fname="naiveOptimization.root",TTree *tr=0, int maxEntries=-1)
{
	if( tr == 0 ) { tr = (TTree*)gROOT->FindObject("TestTree"); }
	
	TFile * fout = TFile::Open(fname,"recreate");
	fout->cd();
	
	TTree * sig  = new TTree("sig","sig");
	TTree * bkg1 = new TTree("bkg1","bkg1");
	TTree * bkg2 = new TTree("bkg2","bkg2");
	
	TTree * tout[] = { sig, bkg1, bkg2 };
	int nclass = 3; // sizeof(tuot) / sizeof(TTree *);

	int classID;
	MultiBDT mvas;
	tr->SetBranchAddress("BDTG", &mvas );
	tr->SetBranchAddress("classID", &classID );
	for(int ii=0; ii<nclass; ++ii) {
		tout[ii]->Branch("mva0",&mvas.Signal);
		tout[ii]->Branch("mva1",&mvas.bg0);
		tout[ii]->Branch("mva2",&mvas.bg1);
	}
	
	for(int ii=0;ii<tr->GetEntries(); ++ii) {
		tr->GetEntry(ii);
		assert( classID < nclass );
		// convert bg0 and bg1 such that signal peaks at 1.
		mvas.bg0 = 1. - mvas.bg0;
		mvas.bg1 = 1. - mvas.bg1;
		if( maxEntries < 0 || tout[classID]->GetEntries() < maxEntries ) {
			tout[classID]->Fill();
		}
	}
		
	for(int ii=0; ii<nclass; ++ii) {
		tout[ii]->Write();
	}
	
	fout->Close();
}
