/// READ HERE ----
/// to run this script you need to load some things in root first:
/// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include -I$CMSSW_BASE/src");
/// .L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so
/// .L ../../libLoopAll.so
/// .L collapseSpinCategories.C+g

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"
#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"

using namespace std;
using namespace RooFit;

void collapseSpinCategories(string infile, int nOldCats=20, int nNewCats=4, int sqrtS=8) {

gROOT->SetBatch();

TFile *inf = TFile::Open(infile.c_str());
RooWorkspace *work=(RooWorkspace*)inf->Get("cms_hgg_workspace");
RooRealVar *mass = (RooRealVar*)work->var("CMS_hgg_mass");

double checkSum=0.;

vector <RooDataSet*> newData;
vector<RooDataHist*> newSig_ggh;
vector<RooDataHist*> newSig_vbf;
vector<RooDataHist*> newSig_zh;
vector<RooDataHist*> newSig_wh;
vector<RooDataHist*> newSig_tth;
vector<RooDataHist*> newSig_gg_grav;
vector<RooDataHist*> newSig_qq_grav;
int newCatCounter=-1;
for (int cat=0; cat<nOldCats; cat++) {

	RooDataSet *data = (RooDataSet*)work->data(Form("data_mass_cat%d",cat));
	RooDataHist *ggh = (RooDataHist*)work->data(Form("roohist_sig_ggh_mass_m125_cat%d",cat));
	RooDataHist *vbf = (RooDataHist*)work->data(Form("roohist_sig_vbf_mass_m125_cat%d",cat));
	RooDataHist *wh = (RooDataHist*)work->data(Form("roohist_sig_wh_mass_m125_cat%d",cat));
	RooDataHist *zh = (RooDataHist*)work->data(Form("roohist_sig_zh_mass_m125_cat%d",cat));
	RooDataHist *tth = (RooDataHist*)work->data(Form("roohist_sig_tth_mass_m125_cat%d",cat));
	RooDataHist *gg_grav = (RooDataHist*)work->data(Form("roohist_sig_gg_grav_mass_m125_cat%d",cat));
	RooDataHist *qq_grav = (RooDataHist*)work->data(Form("roohist_sig_qq_grav_mass_m125_cat%d",cat));
	if (cat%(nOldCats/nNewCats)==0) 
	{ 
		cout << "check -- " << checkSum << endl;
		checkSum=data->sumEntries();
		newCatCounter++;
		RooDataSet *newDataTemp = (RooDataSet*)data->Clone(Form("data_mass_cat%d",newCatCounter));
		RooDataHist *new_ggh = (RooDataHist*)ggh->Clone(Form("roohist_sig_ggh_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_vbf = (RooDataHist*)vbf->Clone(Form("roohist_sig_vbf_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_wh = (RooDataHist*)wh->Clone(Form("roohist_sig_wh_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_zh = (RooDataHist*)zh->Clone(Form("roohist_sig_zh_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_tth = (RooDataHist*)tth->Clone(Form("roohist_sig_tth_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_gg_grav = (RooDataHist*)gg_grav->Clone(Form("roohist_sig_gg_grav_mass_m125_cat%d",newCatCounter));
		RooDataHist *new_qq_grav = (RooDataHist*)qq_grav->Clone(Form("roohist_sig_qq_grav_mass_m125_cat%d",newCatCounter));
		
		newData.push_back(newDataTemp);
		newSig_ggh.push_back(new_ggh);
		newSig_vbf.push_back(new_vbf);
		newSig_wh.push_back(new_wh);
		newSig_zh.push_back(new_zh);
		newSig_tth.push_back(new_tth);
		newSig_gg_grav.push_back(new_gg_grav);
		newSig_qq_grav.push_back(new_qq_grav);
	}
	else
	{
		checkSum += data->sumEntries();
		newData[newCatCounter]->append(*data);
		newSig_ggh[newCatCounter]->add(*ggh);
		newSig_vbf[newCatCounter]->add(*vbf);
		newSig_wh[newCatCounter]->add(*wh);
		newSig_zh[newCatCounter]->add(*zh);
		newSig_tth[newCatCounter]->add(*tth);
		newSig_gg_grav[newCatCounter]->add(*gg_grav);
		newSig_qq_grav[newCatCounter]->add(*qq_grav);
	}
cout<< cat<< endl;
}
cout << "check -- " << checkSum << endl;

	TCanvas *canv = new TCanvas("c","c",1200,800);

	RooWorkspace *newWS = new RooWorkspace("cms_hgg_workspace");
	for (unsigned int i=0; i<newData.size(); i++){
		newWS->import(*newData[i]);	
		newWS->import(*newSig_ggh[i]);	
		newWS->import(*newSig_vbf[i]);	
		newWS->import(*newSig_wh[i]);	
		newWS->import(*newSig_zh[i]);	
		newWS->import(*newSig_tth[i]);	
		newWS->import(*newSig_gg_grav[i]);	
		newWS->import(*newSig_qq_grav[i]);
		cout << "cr check -- " << newData[i]->sumEntries() << endl;
		
		// make bkg
		RooArgList *params = new RooArgList();
		for (int p=0; p<5; p++){
			RooRealVar *var = new RooRealVar(Form("CMS_hgg_quintic%d_%dTeV_cat%d",p,sqrtS,i),"",0.5,-1.,1.);
			params->add(*var);
		}
		RooBernsteinFast<5> *bkg = new RooBernsteinFast<5>(Form("pdf_data_pol_model_%dTeV_cat%d",sqrtS,i),"",*mass,*params);
		RooRealVar *norm = new RooRealVar(Form("pdf_data_pol_model_%dTeV_cat%d_norm",sqrtS,i),"",1000.,0.,1.e7);
		RooExtendPdf *pdf = new RooExtendPdf(Form("data_pol_model_%dTeV_cat%d",sqrtS,i),"",*bkg,*norm);
		pdf->fitTo(*newData[i]);
		newWS->import(*pdf);

		RooPlot *plot = mass->frame();
		newData[i]->plotOn(plot);
		bkg->plotOn(plot);
		plot->Draw();
		canv->Print(Form("fits_%dTeV_cat%d.pdf",sqrtS,i));
	}
	string outfile = infile+"_collapsed.root";
	TFile *hfile = new TFile(outfile.c_str(),"RECREATE");
	newWS->Write();
	newWS->Print();
	hfile->Close();
}
