/// READ HERE ----
/// to run this script you need to load some things in root first:
/// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include -I$CMSSW_BASE/src");
/// .L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so
/// .L ../../libLoopAll.so
/// .L collapseSpinCategories.C+g

#include <iostream>
#include <vector>
#include <string>
#include <map>

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

void collapseSpinCategories(string infile, int nOldCats=20, int nNewCats=4, int sqrtS=8, bool justSignal=false) {

	gROOT->SetBatch();

	TFile *inf = TFile::Open(infile.c_str());
	RooWorkspace *work=(RooWorkspace*)inf->Get("cms_hgg_workspace");
	RooRealVar *mass = (RooRealVar*)work->var("CMS_hgg_mass");
	mass->setBins(160);
	RooRealVar *intL = (RooRealVar*)work->var("IntLumi");

	vector<string> systList;
	systList.push_back("E_res");
	systList.push_back("E_scale");
	systList.push_back("idEff");
	systList.push_back("triggerEff");
	systList.push_back("vtxEff");
	systList.push_back("r9Eff");
	systList.push_back("ptSpin");
	vector<string> procList;
	procList.push_back("ggh");
	procList.push_back("vbf");
	procList.push_back("wh");
	procList.push_back("zh");
	procList.push_back("tth");
	procList.push_back("gg_grav");
	procList.push_back("qq_grav");

	vector <RooDataSet*> newData;
	map<string,vector<RooDataHist*> > newSig;
	map<string,vector<RooDataSet*> > newSigUnbinned;
	vector<pair<double,double> > checkDataSum;
	map<string,vector<pair<double,double> > > checkSigSum;
	
	// initiase check sums etc.
	for (int cat=0; cat<nNewCats; cat++){
		checkDataSum.push_back(make_pair(0.,0.));
	}
	for (vector<string>::iterator it=procList.begin(); it!=procList.end(); it++){
		checkSigSum.insert(make_pair(*it,checkDataSum));
		vector<RooDataHist*> dummy;
		newSig.insert(make_pair(*it,dummy));
		vector<RooDataSet*> dummyUnb;
		newSigUnbinned.insert(make_pair(*it,dummyUnb));
		newSigUnbinned.insert(make_pair(*it+"rv",dummyUnb));
		newSigUnbinned.insert(make_pair(*it+"wv",dummyUnb));
		for (vector<string>::iterator sys=systList.begin(); sys!=systList.end(); sys++){
			checkSigSum.insert(make_pair(*it+*sys+"Up",checkDataSum));
			checkSigSum.insert(make_pair(*it+*sys+"Down",checkDataSum));
			newSig.insert(make_pair(*it+*sys+"Up",dummy));
			newSig.insert(make_pair(*it+*sys+"Down",dummy));
			newSigUnbinned.insert(make_pair(*it+*sys+"Up",dummyUnb));
			newSigUnbinned.insert(make_pair(*it+*sys+"Down",dummyUnb));
		}
	}

	int newCatCounter=-1;
	// loop old cats
	for (int cat=0; cat<nOldCats; cat++) {

		// if it's a new cat then clone
		RooDataSet *data = (RooDataSet*)work->data(Form("data_mass_cat%d",cat));
		if (cat%(nOldCats/nNewCats)==0) 
		{ 
			newCatCounter++;
			if (!justSignal) {
				// data
				checkDataSum[newCatCounter].first = data->sumEntries();
				RooDataSet *newDataTemp = (RooDataSet*)data->Clone(Form("data_mass_cat%d",newCatCounter));
				newData.push_back(newDataTemp);
			}
			// sig
			for (unsigned int iProc=0; iProc<procList.size(); iProc++) {
				RooDataHist *sig = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d",procList[iProc].c_str(),cat));
				checkSigSum[procList[iProc]][newCatCounter].first = sig->sumEntries();
				RooDataHist *newSigTemp = (RooDataHist*)sig->Clone(Form("roohist_sig_%s_mass_m125_cat%d",procList[iProc].c_str(),newCatCounter));
				newSig[procList[iProc]].push_back(newSigTemp);
				RooDataSet *sigUn = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d",procList[iProc].c_str(),cat));
				RooDataSet *newSigUnTemp = (RooDataSet*)sigUn->Clone(Form("sig_%s_mass_m125_cat%d",procList[iProc].c_str(),newCatCounter));
				newSigUnbinned[procList[iProc]].push_back(newSigUnTemp);
				RooDataSet *sigUnRV = (RooDataSet*)work->data(Form("sig_%s_mass_m125_rv_cat%d",procList[iProc].c_str(),cat));
				RooDataSet *newSigUnRVTemp = (RooDataSet*)sigUnRV->Clone(Form("sig_%s_mass_m125_rv_cat%d",procList[iProc].c_str(),newCatCounter));
				newSigUnbinned[procList[iProc]+"rv"].push_back(newSigUnRVTemp);
				RooDataSet *sigUnWV = (RooDataSet*)work->data(Form("sig_%s_mass_m125_wv_cat%d",procList[iProc].c_str(),cat));
				RooDataSet *newSigUnWVTemp = (RooDataSet*)sigUnWV->Clone(Form("sig_%s_mass_m125_wv_cat%d",procList[iProc].c_str(),newCatCounter));
				newSigUnbinned[procList[iProc]+"wv"].push_back(newSigUnWVTemp);
				
				for (vector<string>::iterator sys=systList.begin(); sys!=systList.end(); sys++){
					RooDataHist *sigUp = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					RooDataHist *sigDown = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					checkSigSum[procList[iProc]+*sys+"Up"][newCatCounter].first = sigUp->sumEntries();
					checkSigSum[procList[iProc]+*sys+"Down"][newCatCounter].first = sigDown->sumEntries();
					RooDataHist *newSigUpTemp = (RooDataHist*)sigUp->Clone(Form("roohist_sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),newCatCounter,sys->c_str()));
					RooDataHist *newSigDownTemp = (RooDataHist*)sigDown->Clone(Form("roohist_sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),newCatCounter,sys->c_str()));
					newSig[procList[iProc]+*sys+"Up"].push_back(newSigUpTemp);
					newSig[procList[iProc]+*sys+"Down"].push_back(newSigDownTemp);
					RooDataSet *sigUpUnb = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					RooDataSet *sigDownUnb = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					RooDataSet *newSigUpUnbTemp = (RooDataSet*)sigUpUnb->Clone(Form("sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),newCatCounter,sys->c_str()));
					RooDataSet *newSigDownUnbTemp = (RooDataSet*)sigDownUnb->Clone(Form("sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),newCatCounter,sys->c_str()));
					newSigUnbinned[procList[iProc]+*sys+"Up"].push_back(newSigUpUnbTemp);
					newSigUnbinned[procList[iProc]+*sys+"Down"].push_back(newSigDownUnbTemp);
				}
			}
		}
		// else append
		else
		{
			if (!justSignal) {
				// data
				checkDataSum[newCatCounter].first += data->sumEntries();
				newData[newCatCounter]->append(*data);
			}
			// sig
			for (unsigned int iProc=0; iProc<procList.size(); iProc++) {
				RooDataHist *sig = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d",procList[iProc].c_str(),cat));
				checkSigSum[procList[iProc]][newCatCounter].first += sig->sumEntries();
				newSig[procList[iProc]][newCatCounter]->add(*sig);
				RooDataSet *sigUnb = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d",procList[iProc].c_str(),cat));
				newSigUnbinned[procList[iProc]][newCatCounter]->append(*sigUnb);
				RooDataSet *sigUnbRV = (RooDataSet*)work->data(Form("sig_%s_mass_m125_rv_cat%d",procList[iProc].c_str(),cat));
				newSigUnbinned[procList[iProc]+"rv"][newCatCounter]->append(*sigUnbRV);
				RooDataSet *sigUnbWV = (RooDataSet*)work->data(Form("sig_%s_mass_m125_wv_cat%d",procList[iProc].c_str(),cat));
				newSigUnbinned[procList[iProc]+"wv"][newCatCounter]->append(*sigUnbWV);
				for (vector<string>::iterator sys=systList.begin(); sys!=systList.end(); sys++){
					RooDataHist *sigUp = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					RooDataHist *sigDown = (RooDataHist*)work->data(Form("roohist_sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					checkSigSum[procList[iProc]+*sys+"Up"][newCatCounter].first += sigUp->sumEntries();
					checkSigSum[procList[iProc]+*sys+"Down"][newCatCounter].first += sigDown->sumEntries();
					newSig[procList[iProc]+*sys+"Up"][newCatCounter]->add(*sigUp);
					newSig[procList[iProc]+*sys+"Down"][newCatCounter]->add(*sigDown);
					RooDataSet *sigUpUnb = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d_%sUp01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					RooDataSet *sigDownUnb = (RooDataSet*)work->data(Form("sig_%s_mass_m125_cat%d_%sDown01_sigma",procList[iProc].c_str(),cat,sys->c_str()));
					newSigUnbinned[procList[iProc]+*sys+"Up"][newCatCounter]->append(*sigUpUnb);
					newSigUnbinned[procList[iProc]+*sys+"Down"][newCatCounter]->append(*sigDownUnb);
				}
			}
		}
	}

	TCanvas *canv = new TCanvas("c","c",1200,800);

	RooWorkspace *newWS = new RooWorkspace("cms_hgg_workspace");
	newWS->import(*intL);
	for (int i=0; i<nNewCats; i++){
		if (!justSignal) {
			// import data
			newWS->import(*newData[i]);
			checkDataSum[i].second = newData[i]->sumEntries();
			cout << "cat" << i << Form("  %7s","data") << " -- " << checkDataSum[i].first << " : " << checkDataSum[i].second << endl;
		// import data binned clone
		RooDataHist *binnedData = (RooDataHist*)newData[i]->binnedClone(Form("roohist_data_mass_cat%d",i));
		newWS->import(*binnedData);
		}
		for (unsigned int iProc=0; iProc<procList.size(); iProc++) {
			newWS->import(*newSig[procList[iProc]][i]);
			newWS->import(*newSigUnbinned[procList[iProc]][i]);
			newWS->import(*newSigUnbinned[procList[iProc]+"rv"][i]);
			newWS->import(*newSigUnbinned[procList[iProc]+"wv"][i]);
			checkSigSum[procList[iProc]][i].second = newSig[procList[iProc]][i]->sumEntries();
			cout << "cat" << i << Form("  %7s",procList[iProc].c_str()) << " -- " << checkSigSum[procList[iProc]][i].first << " : " << checkSigSum[procList[iProc]][i].second << endl;
			for (vector<string>::iterator sys=systList.begin(); sys!=systList.end(); sys++){
				newWS->import(*newSig[procList[iProc]+*sys+"Up"][i]);
				newWS->import(*newSig[procList[iProc]+*sys+"Down"][i]);
				newWS->import(*newSigUnbinned[procList[iProc]+*sys+"Up"][i]);
				newWS->import(*newSigUnbinned[procList[iProc]+*sys+"Down"][i]);
				checkSigSum[procList[iProc]+*sys+"Up"][i].second = newSig[procList[iProc]+*sys+"Up"][i]->sumEntries();
				checkSigSum[procList[iProc]+*sys+"Down"][i].second = newSig[procList[iProc]+*sys+"Down"][i]->sumEntries();
				cout << "cat" << i << "  " << Form("%10s",sys->c_str()) << "Up   -- " << checkSigSum[procList[iProc]+*sys+"Up"][i].first << " : " << checkSigSum[procList[iProc]+*sys+"Up"][i].second << endl;
				cout << "cat" << i << "  " << Form("%10s",sys->c_str()) << "Down -- " << checkSigSum[procList[iProc]+*sys+"Down"][i].first << " : " << checkSigSum[procList[iProc]+*sys+"Down"][i].second << endl;
			}
		}
	
		// make bkg
		if (!justSignal) {
			RooArgList *params = new RooArgList();
			for (int p=0; p<5; p++){
				RooRealVar *var = new RooRealVar(Form("CMS_hgg_quintic%d_%dTeV_cat%d",p,sqrtS,i),"",0.5,-1.,1.);
				RooFormulaVar *var2 = new RooFormulaVar(Form("CMS_hgg_modquintic%d_%dTeV_cat%d",p,sqrtS,i),"","@0*@0",RooArgList(*var));
				params->add(*var2);
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
	}
	string outfile = infile+"_collapsed.root";
	TFile *hfile = new TFile(outfile.c_str(),"RECREATE");
	newWS->Write();
	newWS->Print();
	hfile->Close();

	cout << "File written. Check sums: " << endl;
	for (int i=0; i<nNewCats; i++){
		cout << "cat" << i << Form("  %7s","data") << " -- " << checkDataSum[i].first << " : " << checkDataSum[i].second << endl;
		if (TMath::Abs(checkDataSum[i].first-checkDataSum[i].second)>1.e-6) {
			cout << "ERROR -- Check sums do not match!" << endl;
			exit(1);
		}
		for (unsigned int iProc=0; iProc<procList.size(); iProc++) {
			cout << "cat" << i << Form("  %7s",procList[iProc].c_str()) << " -- " << checkSigSum[procList[iProc]][i].first << " : " << checkSigSum[procList[iProc]][i].second << endl;
			if (TMath::Abs(checkSigSum[procList[iProc]][i].first-checkSigSum[procList[iProc]][i].second)>1.e-6) {
				cout << "ERROR -- Check sums do not match!" << endl;
				exit(1);
			}
			for (vector<string>::iterator sys=systList.begin(); sys!=systList.end(); sys++){
				cout << "cat" << i << "  " << Form("%10s",sys->c_str()) << "Up   -- " << checkSigSum[procList[iProc]+*sys+"Up"][i].first << " : " << checkSigSum[procList[iProc]+*sys+"Up"][i].second << endl;
				if (TMath::Abs(checkSigSum[procList[iProc]+*sys+"Up"][i].first-checkSigSum[procList[iProc]+*sys+"Up"][i].second)>1.e-6) {
					cout << "ERROR -- Check sums do not match!" << endl;
					exit(1);
				}
				cout << "cat" << i << "  " << Form("%10s",sys->c_str()) << "Down -- " << checkSigSum[procList[iProc]+*sys+"Down"][i].first << " : " << checkSigSum[procList[iProc]+*sys+"Down"][i].second << endl;
				if (TMath::Abs(checkSigSum[procList[iProc]+*sys+"Down"][i].first-checkSigSum[procList[iProc]+*sys+"Down"][i].second)>1.e-6) {
					cout << "ERROR -- Check sums do not match!" << endl;
					exit(1);
				}
			}
		}
	}
}
