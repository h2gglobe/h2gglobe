#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TF1.h"
#include "RooPlot.h"
#include "RooVoigtian.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "boost/algorithm/string/erase.hpp"

#include "../interface/FinalModelConstruction.h"

using namespace std;
using namespace RooFit;
using namespace boost;

FinalModelConstruction::FinalModelConstruction(RooRealVar *massVar, RooRealVar *MHvar, RooRealVar *intL, int mhLow, int mhHigh, string proc, int cat, bool doSecMods, string systematicsFileName, int verbosity, bool isCB, bool is2011):
  mass(massVar),
  MH(MHvar),
  intLumi(intL),
  mhLow_(mhLow),
  mhHigh_(mhHigh),
  proc_(proc),
  cat_(cat),
  doSecondaryModels(doSecMods),
  isCutBased_(isCB),
	is2011_(is2011),
  verbosity_(verbosity),
  systematicsSet_(false),
  rvFractionSet_(false)
{
  allMH_ = getAllMH();
	if (is2011_) sqrts_ = 7;
	else sqrts_ = 8;
  // load xs and br info from Normalization_8TeV
  norm = new Normalization_8TeV();
  norm->Init(sqrts_);
  TGraph *brGraph = norm->GetBrGraph();
	brSpline = graphToSpline(Form("fbr_%dTeV",sqrts_),brGraph);
  
  string procs[8] = {"ggh","vbf","wzh","wh","zh","tth","gg_grav","qq_grav"};
  for (int i=0; i<8; i++){
    TGraph *xsGraph = norm->GetSigmaGraph(procs[i].c_str());
    RooSpline1D *xsSpline = graphToSpline(Form("fxs_%s_%dTeV",procs[i].c_str(),sqrts_),xsGraph);
    xsSplines.insert(pair<string,RooSpline1D*>(procs[i],xsSpline));
  }

  // const smearing terms for resolution systematic hard-coded here!!
  //constSmearVals.push_back(0.0080838); // cat0
  //constSmearVals.push_back(0.0083661); // cat1
  //constSmearVals.push_back(0.0086261); // cat2
  //constSmearVals.push_back(0.0134325); // cat3
  //constSmearVals.push_back(0.0119792); // cat4 (all exc cats)
  loadSignalSystematics(systematicsFileName);
	if (verbosity_) printSignalSystematics();
	// should pick up different values for the cut-based but not sure what they are!
}

FinalModelConstruction::~FinalModelConstruction(){}

void FinalModelConstruction::printSignalSystematics(){
	
	// nuisance parameters
	cout << "Implementing the following floating nuisance parameters" << endl;
	for (map<string,RooRealVar*>::iterator sys=photonSystematics.begin(); sys!=photonSystematics.end(); sys++){
		cout << "\t " << sys->first << " -- "; sys->second->Print();
	}
	// const parameters
	cout << "Implementing the following constant parameters" << endl;
	for (map<string,RooRealVar*>::iterator sys=photonSystematicConsts.begin(); sys!=photonSystematicConsts.end(); sys++){
		cout << "\t " << sys->first << " -- "; sys->second->Print();
	}
}

void FinalModelConstruction::loadSignalSystematics(string filename){
	
	int diphotonCat=-1;
	string proc;
	ifstream datfile;
	datfile.open(filename.c_str());
	if (datfile.fail()) exit(1);
	while (datfile.good()){
		
		string line;
		getline(datfile,line);
		
		// The input file needs correct ordering
		if (line=="\n" || line.substr(0,1)=="#" || line==" " || line.empty()) continue;
		// First read the photon categories
		if (starts_with(line,"photonCats")){
			line = line.substr(line.find("=")+1,string::npos);
			split(photonCats,line,boost::is_any_of(","));
			if (verbosity_){
				cout << "PhotonVec: ";
				if (verbosity_) printVec(photonCats);
			}
		}
		// Then read diphoton cat
		else if (starts_with(line,"diphotonCat")){
			diphotonCat = boost::lexical_cast<int>(line.substr(line.find("=")+1,string::npos));
		}
		// And the process
		else if (starts_with(line,"proc")){
			proc = line.substr(line.find('=')+1,string::npos);
			if (verbosity_) cout << "Process:  " << proc << "  DiphoCat: " << diphotonCat << endl;
		}
		// Then read values
		else {
			stripSpace(line);
			vector<string> els;
			split(els,line,boost::is_any_of(" "));
			if (verbosity_) {cout << "\t"; printVec(els);}
			if (els.size()!=4) {
				cout << "I cant read this datfile " << endl;
				exit(1);
			}
			string phoSystName = els[0];
			double meanCh = lexical_cast<double>(els[1]);
			if( meanCh != 0. ) { 
				addToSystMap(meanSysts,proc,diphotonCat,phoSystName,meanCh);
				RooRealVar *meanChVar = new RooRealVar(Form("const_%s_cat%d_mean_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),Form("const_%s_cat%d_mean_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),meanCh);
				meanChVar->setConstant(true);
				photonSystematicConsts.insert(make_pair(meanChVar->GetName(),meanChVar));
			}
			double sigmaCh = lexical_cast<double>(els[2]);
			if( sigmaCh != 0. ) {
				addToSystMap(sigmaSysts,proc,diphotonCat,phoSystName,sigmaCh);
				RooRealVar *sigmaChVar = new RooRealVar(Form("const_%s_cat%d_sigma_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),Form("const_%s_cat%d_sigma_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),sigmaCh);
				sigmaChVar->setConstant(true);
				photonSystematicConsts.insert(make_pair(sigmaChVar->GetName(),sigmaChVar));
			}
			double rateCh = lexical_cast<double>(els[3]);
			if( rateCh != 0. ) {
				addToSystMap(rateSysts,proc,diphotonCat,phoSystName,rateCh);
				RooRealVar *rateChVar = new RooRealVar(Form("const_%s_cat%d_rate_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),Form("const_%s_cat%d_rate_%s",proc.c_str(),diphotonCat,phoSystName.c_str()),rateCh);
				rateChVar->setConstant(true);
				photonSystematicConsts.insert(make_pair(rateChVar->GetName(),rateChVar));
			}
		}
	}
	datfile.close();
	for (vector<string>::iterator it=photonCats.begin(); it!=photonCats.end(); it++){
		RooRealVar *varScale = new RooRealVar(Form("CMS_hgg_nuisance%sscale",it->c_str()),Form("CMS_hgg_nuisance%sscale",it->c_str()),0.,-5.,5.);
		RooRealVar *varSmear = new RooRealVar(Form("CMS_hgg_nuisance%ssmear",it->c_str()),Form("CMS_hgg_nuisance%ssmear",it->c_str()),0.,-5.,5.);
		photonSystematics.insert(make_pair(varScale->GetName(),varScale));
		photonSystematics.insert(make_pair(varSmear->GetName(),varSmear));
	}
}

void FinalModelConstruction::stripSpace(string &line){
	stringstream lineSt(line);
	line="";
	string word;
	while (lineSt >> word) {
		line.append(word).append(" ");
	}
	line = line.substr(0,line.size()-1);
}

void FinalModelConstruction::printVec(vector<string> vec){

	cout << "[";
	for (unsigned i=0; i<vec.size()-1;i++){
		cout << vec[i] << ",";
	}
	cout << vec[vec.size()-1] << "]" << endl;
}

void FinalModelConstruction::printSystMap(map<string,map<int,map<string,double> > > &theMap){

	for (map<string,map<int,map<string,double> > >::iterator p = theMap.begin(); p != theMap.end(); p++) {
		for (map<int,map<string,double> >::iterator c = p->second.begin(); c != p->second.end(); c++){
			cout << "Proc = " << p->first << "  DiphotonCat: " << c->first << endl;
			for (map<string,double>::iterator m = c->second.begin(); m != c->second.end(); m++){
				cout << "\t " << m->first << " : " << m->second << endl;
			}
		}
	}
}

void FinalModelConstruction::addToSystMap(map<string,map<int,map<string,double> > > &theMap, string proc, int diphotonCat, string phoSystName, double var){
	// does proc map exist?
	if (theMap.find(proc)!=theMap.end()) {
		// does the cat map exist?
		if (theMap[proc].find(diphotonCat)!=theMap[proc].end()){
			theMap[proc][diphotonCat].insert(make_pair(phoSystName,var));
		}
		else{
			map<string,double> temp;
			temp.insert(make_pair(phoSystName,var));
			theMap[proc].insert(make_pair(diphotonCat,temp));
		}
	}
	else {
		map<string,double> temp;
		map<int,map<string,double> > cTemp;
		temp.insert(make_pair(phoSystName,var));
		cTemp.insert(make_pair(diphotonCat,temp));
		theMap.insert(make_pair(proc,cTemp));
	}
}

vector<int> FinalModelConstruction::getAllMH(){
  vector<int> result;
  for (int m=mhLow_; m<=mhHigh_; m+=5){
    if (verbosity_>=1) cout << "FinalModelConstruction - Adding mass: " << m << endl;
    result.push_back(m);
  }
  return result;
}

void FinalModelConstruction::setSecondaryModelVars(RooRealVar *mh_sm, RooRealVar *deltam, RooAddition *mh_2, RooRealVar *width){
  MH_SM = mh_sm;
  DeltaM = deltam;
  MH_2 = mh_2;
  higgsDecayWidth = width;
  
  TGraph *brGraph = norm->GetBrGraph();
	brSpline_SM = graphToSpline(Form("fbr_%dTeV_SM",sqrts_),brGraph,MH_SM);
	brSpline_2 = graphToSpline(Form("fbr_%dTeV_2",sqrts_),brGraph,MH_2);
	brSpline_NW = graphToSpline(Form("fbr_%dTeV_NW",sqrts_),brGraph,MH);
  
	string procs[8] = {"ggh","vbf","wzh","wh","zh","tth","gg_grav","qq_grav"};
  for (int i=0; i<8; i++){
    TGraph *xsGraph = norm->GetSigmaGraph(procs[i].c_str());
    RooSpline1D *xsSpline_SM = graphToSpline(Form("fxs_%s_%dTeV_SM",procs[i].c_str(),sqrts_),xsGraph,MH_SM);
    RooSpline1D *xsSpline_2 = graphToSpline(Form("fxs_%s_%dTeV_2",procs[i].c_str(),sqrts_),xsGraph,MH_2); 
    RooSpline1D *xsSpline_NW = graphToSpline(Form("fxs_%s_%dTeV_NW",procs[i].c_str(),sqrts_),xsGraph,MH);
    xsSplines_SM.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_SM));
    xsSplines_2.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_2));
    xsSplines_NW.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_NW));
  }
  secondaryModelVarsSet=true;
}

void FinalModelConstruction::getRvFractionFunc(string name){
  
  assert(allMH_.size()==rvDatasets.size());
  assert(allMH_.size()==wvDatasets.size());
  vector<double> mhValues, rvFracValues;
  for (unsigned int i=0; i<allMH_.size(); i++){
    int mh = allMH_[i];
    mhValues.push_back(mh);
    double rvN = rvDatasets[mh]->sumEntries();
    double wvN = wvDatasets[mh]->sumEntries();
		double rvF = rvN/(rvN+wvN);
		if (rvF != rvF) rvF=1.; // incase nan when no entries
    rvFracValues.push_back(rvF);
  }
  rvFracFunc = new RooSpline1D(name.c_str(),name.c_str(),*MH,mhValues.size(),&(mhValues[0]),&(rvFracValues[0]));
  rvFractionSet_=true;
}

RooAbsReal* FinalModelConstruction::getMeanWithPhotonSyst(RooAbsReal *dm, string name){
	
	string formula="(@0+@1)*(1.+@2";
	RooArgList *dependents = new RooArgList();
	dependents->add(*MH); // MH sits at @0
	dependents->add(*dm); // dm sits at @1
	dependents->add(*globalScale); // sits at @2

	for (unsigned int i=0; i<photonCats.size(); i++){
		string phoCat = photonCats[i];
		int formPlace = dependents->getSize();		
		bool hasEffect = false;
		// want scale effect on mean
		if( photonSystematicConsts.find(Form("const_%s_cat%d_mean_%sscale",proc_.c_str(),cat_,phoCat.c_str())) 
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvScale = photonSystematicConsts[Form("const_%s_cat%d_mean_%sscale",proc_.c_str(),cat_,phoCat.c_str())]; 
			RooRealVar *nuisScale = photonSystematics[Form("CMS_hgg_nuisance%sscale",phoCat.c_str())];
			
			// and res effect on mean
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			dependents->add(*cvScale);
			dependents->add(*nuisScale);
			formPlace += 2;
			hasEffect = true;
		}
		// check if smearing affects mean
		if( photonSystematicConsts.find(Form("const_%s_cat%d_mean_%ssmear",proc_.c_str(),cat_,phoCat.c_str())) 
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvSmear = photonSystematicConsts[Form("const_%s_cat%d_mean_%ssmear",proc_.c_str(),cat_,phoCat.c_str())];
			RooRealVar *nuisSmear = photonSystematics[Form("CMS_hgg_nuisance%ssmear",phoCat.c_str())];
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			dependents->add(*cvSmear);
			dependents->add(*nuisSmear);
			hasEffect = true;
		}
		if( ! hasEffect ) {
			std::cerr << "WARNING: Photon category systematic " << phoCat 
				  << " doesn't affect the signal model scale." << std::endl; 
		}
	}
	formula+=")";
	RooFormulaVar *formVar = new RooFormulaVar(name.c_str(),name.c_str(),formula.c_str(),*dependents);
	return formVar;
}

RooAbsReal* FinalModelConstruction::getSigmaWithPhotonSyst(RooAbsReal *sig_fit, string name){

	string formula="@0*(1.";
	RooArgList *dependents = new RooArgList();
	dependents->add(*sig_fit); // sig_fit sits at @0
	
	for (unsigned int i=0; i<photonCats.size(); i++){
		string phoCat = photonCats[i];
		int formPlace = dependents->getSize();
		bool hasEffect = false;
		// want scale effect on sigma
		if( photonSystematicConsts.find(Form("const_%s_cat%d_sigma_%sscale",proc_.c_str(),cat_,phoCat.c_str()))
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvScale = photonSystematicConsts[Form("const_%s_cat%d_sigma_%sscale",proc_.c_str(),cat_,phoCat.c_str())]; 
			RooRealVar *nuisScale = photonSystematics[Form("CMS_hgg_nuisance%sscale",phoCat.c_str())];
			dependents->add(*cvScale);
			dependents->add(*nuisScale);
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			formPlace += 2;
			hasEffect = true;
		}
		
		// and res effect on sigma
		if( photonSystematicConsts.find(Form("const_%s_cat%d_sigma_%ssmear",proc_.c_str(),cat_,phoCat.c_str()))
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvSmear = photonSystematicConsts[Form("const_%s_cat%d_sigma_%ssmear",proc_.c_str(),cat_,phoCat.c_str())];
			RooRealVar *nuisSmear = photonSystematics[Form("CMS_hgg_nuisance%ssmear",phoCat.c_str())];
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			dependents->add(*cvSmear);
			dependents->add(*nuisSmear);
			hasEffect = true;
		}
		if( ! hasEffect ) {
			std::cerr << "WARNING: Photon category systematic " << 
				phoCat << " doesn't affect the signal model width." << std::endl; 
		}
	}
	formula+=")";
	formula = Form("TMath::Max(%s,0.)",formula.c_str()); // consider smooth cutoff ? 
	RooFormulaVar *formVar = new RooFormulaVar(name.c_str(),name.c_str(),formula.c_str(),*dependents);
	return formVar;
}

RooAbsReal* FinalModelConstruction::getRateWithPhotonSyst(string name){
	
	string formula="(1.";
	RooArgList *dependents = new RooArgList();

	for (unsigned int i=0; i<photonCats.size(); i++){
		string phoCat = photonCats[i];
		int formPlace = dependents->getSize();
		bool hasEffect = false;
		// want scale effect on rate
		if( photonSystematicConsts.find(Form("const_%s_cat%d_rate_%sscale",proc_.c_str(),cat_,phoCat.c_str()))
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvScale = photonSystematicConsts[Form("const_%s_cat%d_rate_%sscale",proc_.c_str(),cat_,phoCat.c_str())]; 
			RooRealVar *nuisScale = photonSystematics[Form("CMS_hgg_nuisance%sscale",phoCat.c_str())];
			dependents->add(*cvScale);
			dependents->add(*nuisScale);
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			formPlace += 2;
			hasEffect = true;
		}
		// and res effect on rate
		if( photonSystematicConsts.find(Form("const_%s_cat%d_rate_%ssmear",proc_.c_str(),cat_,phoCat.c_str()))
		    != photonSystematicConsts.end() ) {
			RooRealVar *cvSmear = photonSystematicConsts[Form("const_%s_cat%d_rate_%ssmear",proc_.c_str(),cat_,phoCat.c_str())];
			RooRealVar *nuisSmear = photonSystematics[Form("CMS_hgg_nuisance%ssmear",phoCat.c_str())];
			//formula += Form("@%d*@%d+@%d*@%d",formPlace,formPlace+1,formPlace+2,formPlace+3);
			dependents->add(*cvSmear);
			dependents->add(*nuisSmear);
			formula += Form("+@%d*@%d",formPlace,formPlace+1);
			hasEffect = true;
		}
		if( ! hasEffect ) {
			std::cerr << "WARNING: Photon category systematic " << 
				phoCat << " don't affect the signal normalization." << std::endl; 
		}
	}
	formula+=")";	
	if (isCutBased_){
		// category hard code should be tidied up
		// should be actually taken out and be put in the datfile
		if (cat_==0 || cat_==1) {
			formula += Form("*(1.+@%d)",dependents->getSize());
			dependents->add(*r9barrelNuisance);
		}
		if (cat_==2 || cat_==3) {
			formula += Form("*(1.+@%d)",dependents->getSize());
			dependents->add(*r9mixedNuisance);
		}
	}
	RooFormulaVar *formVar = new RooFormulaVar(name.c_str(),name.c_str(),formula.c_str(),*dependents);
	return formVar;
}

void FinalModelConstruction::setupSystematics(){
  
  //int nuisCat = cat_;
  // this is a hack for correlation with previous 2011 ws
  // for legacy paper - this MUST BE UPDATED
  //if (cat_>=nIncCats_) nuisCat = nIncCats_;
  //categoryScale = new RooRealVar(Form("CMS_hgg_nuissancedeltamcat%d",nuisCat),Form("CMS_hgg_nuissancedeltamcat%d",nuisCat),0.,-5.,5.);
  //categoryScale->setConstant(true);
  //categorySmear = new RooConstVar(Form("CMS_hgg_constsmearcat%d",nuisCat),Form("CMS_hgg_constsmearcat%d",nuisCat),constSmearVals[nuisCat]);
  //categoryResolution = new RooRealVar(Form("CMS_hgg_nuissancedeltasmearcat%d",nuisCat),Form("CMS_hgg_nuissancedeltasmearcat%d",nuisCat),0.0,-0.2,0.2);
	
	vertexNuisance = new RooRealVar(Form("CMS_hgg_nuisancedeltafracright_%dTeV",sqrts_),Form("CMS_hgg_nuisancedeltafracright_%dTeV",sqrts_),0.,-1.,1.);
	vertexNuisance->setConstant(true);
	globalScale = new RooRealVar("CMS_hgg_globalscale","CMS_hgg_globalscale",0.,-5.,5.);
	globalScale->setConstant(true);
	if (isCutBased_) {
		r9barrelNuisance = new RooRealVar(Form("CMS_hgg_nuisancedeltar9barrel_%dTeV",sqrts_),Form("CMS_hgg_nuisancedeltar9barrel_%dTeV",sqrts_),0.,-1.,1.);
		r9mixedNuisance = new RooRealVar(Form("CMS_hgg_nuisancedeltar9mixed_%dTeV",sqrts_),Form("CMS_hgg_nuisancedeltar9mixed_%dTeV",sqrts_),0.,-1.,1.);
		r9barrelNuisance->setConstant(true);
		r9mixedNuisance->setConstant(true);
	}
  systematicsSet_=true;
}

void FinalModelConstruction::buildStdPdf(string name, int nGaussians, bool recursive){

  if (!systematicsSet_) setupSystematics();
  vector<RooAddPdf*> pdfs;
  pdfs = buildPdf(name.c_str(),nGaussians,recursive,stdSplines);
  finalPdf = pdfs[0];
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    finalPdf_SM = pdfs[1];
    finalPdf_2 = pdfs[2];
    finalPdf_NW = pdfs[3];
  }
}

void FinalModelConstruction::buildRvWvPdf(string name, int nGrv, int nGwv, bool recursive){

  if (!rvFractionSet_) getRvFractionFunc(Form("%s_%s_cat%d_rvFracFunc",name.c_str(),proc_.c_str(),cat_));
  if (!systematicsSet_) setupSystematics();
  RooFormulaVar *rvFraction = new RooFormulaVar(Form("%s_%s_cat%d_rvFrac",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_rvFrac",name.c_str(),proc_.c_str(),cat_),"TMath::Min(@0+@1,1.0)",RooArgList(*vertexNuisance,*rvFracFunc));
  vector<RooAddPdf*> rvPdfs = buildPdf(name,nGrv,recursive,rvSplines,Form("_rv_%dTeV",sqrts_)); 
  vector<RooAddPdf*> wvPdfs = buildPdf(name,nGwv,recursive,wvSplines,Form("_wv_%dTeV",sqrts_)); 
  finalPdf = new RooAddPdf(Form("%s_%s_cat%d",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[0],*wvPdfs[0]),RooArgList(*rvFraction));
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    finalPdf_SM = new RooAddPdf(Form("%s_%s_cat%d_SM",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_SM",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[1],*wvPdfs[1]),RooArgList(*rvFraction));
    finalPdf_2 = new RooAddPdf(Form("%s_%s_cat%d_2",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_2",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[2],*wvPdfs[2]),RooArgList(*rvFraction));
    finalPdf_NW = new RooAddPdf(Form("%s_%s_cat%d_NW",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_NW",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[3],*wvPdfs[3]),RooArgList(*rvFraction));
  }
}

vector<RooAddPdf*> FinalModelConstruction::buildPdf(string name, int nGaussians, bool recursive, map<string,RooSpline1D*> splines, string add){
  
  vector<RooAddPdf*> result;

  RooArgList *gaussians = new RooArgList();
  RooArgList *coeffs = new RooArgList();
  string ext = Form("%s_cat%d%s",proc_.c_str(),cat_,add.c_str());
  
  // for SM Higgs as Background
  RooArgList *gaussians_SM = new RooArgList();
  RooArgList *coeffs_SM = new RooArgList();
  
  // for Second Higgs
  RooArgList *gaussians_2 = new RooArgList();
  RooArgList *coeffs_2 = new RooArgList();

  // for Natural Width
  RooArgList *voigtians_NW = new RooArgList();
  RooArgList *coeffs_NW = new RooArgList();

  for (int g=0; g<nGaussians; g++){
    RooAbsReal *dm = splines[Form("dm_g%d",g)];
    dm->SetName(Form("dm_g%d_%s",g,ext.c_str()));
    RooAbsReal *mean = getMeanWithPhotonSyst(dm,Form("mean_g%d_%s",g,ext.c_str()));
		//RooAbsReal *mean = new RooFormulaVar(Form("mean_g%d_%s",g,ext.c_str()),Form("mean_g%d_%s",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH,*dm,*globalScale,*categoryScale));
    RooAbsReal *sig_fit = splines[Form("sigma_g%d",g)];
    sig_fit->SetName(Form("sigma_g%d_%s",g,ext.c_str()));
    RooAbsReal *sigma = getSigmaWithPhotonSyst(sig_fit,Form("sig_g%d_%s",g,ext.c_str()));
    //RooAbsReal *sigma = new RooFormulaVar(Form("sig_g%d_%s",g,ext.c_str()),Form("sig_g%d_%s",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fit,*MH,*categorySmear,*categoryResolution));
		RooGaussian *gaus = new RooGaussian(Form("gaus_g%d_%s",g,ext.c_str()),Form("gaus_g%d_%s",g,ext.c_str()),*mass,*mean,*sigma);
    gaussians->add(*gaus);
    // add secondary models as well
    if (doSecondaryModels){
      assert(secondaryModelVarsSet);
      // sm higgs as background
      RooAbsReal *dmSM = splines[Form("dm_g%d_SM",g)];
      dmSM->SetName(Form("dm_g%d_%s_SM",g,ext.c_str()));
			RooAbsReal *meanSM = getMeanWithPhotonSyst(dmSM,Form("mean_g%d_%s_SM",g,ext.c_str()));
      //RooAbsReal *meanSM = new RooFormulaVar(Form("mean_g%d_%s_SM",g,ext.c_str()),Form("mean_g%d_%s_SM",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH_SM,*dmSM,*globalScale,*categoryScale));
      RooAbsReal *sig_fitSM = splines[Form("sigma_g%d_SM",g)];
      sig_fitSM->SetName(Form("sigma_g%d_%s_SM",g,ext.c_str()));
			RooAbsReal *sigmaSM = getSigmaWithPhotonSyst(sig_fitSM,Form("sig_g%d_%s_SM",g,ext.c_str()));
      //RooAbsReal *sigmaSM = new RooFormulaVar(Form("sig_g%d_%s_SM",g,ext.c_str()),Form("sig_g%d_%s_SM",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fitSM,*MH_SM,*categorySmear,*categoryResolution));
      RooGaussian *gausSM = new RooGaussian(Form("gaus_g%d_%s_SM",g,ext.c_str()),Form("gaus_g%d_%s_SM",g,ext.c_str()),*mass,*meanSM,*sigmaSM);
      gaussians_SM->add(*gausSM);
      // second degen higgs
      RooAbsReal *dm2 = splines[Form("dm_g%d_2",g)];
      dm2->SetName(Form("dm_g%d_%s_2",g,ext.c_str()));
			RooAbsReal *mean2 = getMeanWithPhotonSyst(dm2,Form("mean_g%d_%s_2",g,ext.c_str()));
      //RooAbsReal *mean2 = new RooFormulaVar(Form("mean_g%d_%s_2",g,ext.c_str()),Form("mean_g%d_%s_2",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH_2,*dm2,*globalScale,*categoryScale));
      RooAbsReal *sig_fit2 = splines[Form("sigma_g%d_2",g)];
      sig_fit2->SetName(Form("sigma_g%d_%s_2",g,ext.c_str()));
			RooAbsReal *sigma2 = getSigmaWithPhotonSyst(sig_fit2,Form("sig_g%d_%s_2",g,ext.c_str()));
      //RooAbsReal *sigma2 = new RooFormulaVar(Form("sig_g%d_%s_2",g,ext.c_str()),Form("sig_g%d_%s_2",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fit2,*MH_2,*categorySmear,*categoryResolution));
      RooGaussian *gaus2 = new RooGaussian(Form("gaus_g%d_%s_2",g,ext.c_str()),Form("gaus_g%d_%s_2",g,ext.c_str()),*mass,*mean2,*sigma2);
      gaussians_2->add(*gaus2);
      // natural width
      RooVoigtian *voigNW = new RooVoigtian(Form("voig_g%d_%s_NW",g,ext.c_str()),Form("voig_g%d_%s_NW",g,ext.c_str()),*mass,*mean,*higgsDecayWidth,*sigma);
      voigtians_NW->add(*voigNW);
    }
    if (g<nGaussians-1) {
      RooAbsReal *frac = splines[Form("frac_g%d",g)];
      frac->SetName(Form("frac_g%d_%s",g,ext.c_str()));
      coeffs->add(*frac);
      // add secondary models as well
      if (doSecondaryModels){
        assert(secondaryModelVarsSet);
        // sm higgs as background
        RooAbsReal *fracSM = splines[Form("frac_g%d_SM",g)];
        fracSM->SetName(Form("frac_g%d_%s_SM",g,ext.c_str()));
        coeffs_SM->add(*fracSM);
        // second degen higgs
        RooAbsReal *frac2 = splines[Form("frac_g%d_2",g)];
        frac2->SetName(Form("frac_g%d_%s_2",g,ext.c_str()));
        coeffs_2->add(*frac2);
        // natural width
        coeffs_NW->add(*frac);
      }
    }
  }
  assert(gaussians->getSize()==nGaussians && coeffs->getSize()==nGaussians-1);
  RooAddPdf *pdf = new RooAddPdf(Form("%s_%s",name.c_str(),ext.c_str()),Form("%s_%s",name.c_str(),ext.c_str()),*gaussians,*coeffs,recursive);
  result.push_back(pdf);
  
  // add secondary models as well
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    // sm higgs as background
    RooAddPdf *pdf_SM = new RooAddPdf(Form("%s_%s_SM",name.c_str(),ext.c_str()),Form("%s_%s_SM",name.c_str(),ext.c_str()),*gaussians_SM,*coeffs_SM,recursive);
    result.push_back(pdf_SM);
    // second degen higgs
    RooAddPdf *pdf_2 = new RooAddPdf(Form("%s_%s_2",name.c_str(),ext.c_str()),Form("%s_%s_2",name.c_str(),ext.c_str()),*gaussians_2,*coeffs_2,recursive);
    result.push_back(pdf_2);
    // natural width
    RooAddPdf *pdf_NW = new RooAddPdf(Form("%s_%s_NW",name.c_str(),ext.c_str()),Form("%s_%s_NW",name.c_str(),ext.c_str()),*voigtians_NW,*coeffs_NW,recursive);
    result.push_back(pdf_NW);
  }

  return result;
}

void FinalModelConstruction::setRVsplines(map<string,RooSpline1D*> splines){
  rvSplines = splines;
}

void FinalModelConstruction::setWVsplines(map<string,RooSpline1D*> splines){
  wvSplines = splines;
}

void FinalModelConstruction::setSTDsplines(map<string,RooSpline1D*> splines){
  stdSplines = splines;
}

void FinalModelConstruction::setRVdatasets(map<int,RooDataSet*> data){
  rvDatasets = data;
}

void FinalModelConstruction::setWVdatasets(map<int,RooDataSet*> data){
  wvDatasets = data;
}

void FinalModelConstruction::setSTDdatasets(map<int,RooDataSet*> data){
  stdDatasets = data;
}

void FinalModelConstruction::makeSTDdatasets(){
  for (unsigned int i=0; i<allMH_.size(); i++){
    int mh=allMH_[i];
		RooDataSet *data = (RooDataSet*)rvDatasets[mh]->Clone(Form("sig_%s_mass_m%d_cat%d",proc_.c_str(),mh,cat_));
		data->append(*wvDatasets[mh]);
		stdDatasets.insert(pair<int,RooDataSet*>(mh,data));
	}	
}

void FinalModelConstruction::plotPdf(string outDir){
  system(Form("mkdir -p %s",outDir.c_str()));
  
  TCanvas *canv = new TCanvas();
  RooPlot *dataPlot = mass->frame(Title(Form("%s_cat%d",proc_.c_str(),cat_)),Range(100,160));
  for (unsigned int i=0; i<allMH_.size(); i++){
    int mh=allMH_[i];
    stdDatasets[mh]->plotOn(dataPlot,Binning(160));
    MH->setVal(mh);
    extendPdf->plotOn(dataPlot);
  }
  dataPlot->Draw();
  canv->Print(Form("%s/%s_cat%d_fits.pdf",outDir.c_str(),proc_.c_str(),cat_));
  canv->Print(Form("%s/%s_cat%d_fits.png",outDir.c_str(),proc_.c_str(),cat_));
  
  RooPlot *pdfPlot = mass->frame(Title(Form("%s_cat%d",proc_.c_str(),cat_)),Range(100,160));
	pdfPlot->GetYaxis()->SetTitle(Form("Pdf projection / %2.1f GeV",(mass->getMax()-mass->getMin())/160.));
  for (int mh=mhLow_; mh<=mhHigh_; mh++){
    MH->setVal(mh);
		// to get correct normlization need to manipulate with bins and range
    extendPdf->plotOn(pdfPlot,Normalization(mass->getBins()/160.*(mass->getMax()-mass->getMin())/60.,RooAbsReal::RelativeExpected));
  }
  pdfPlot->Draw();
  canv->Print(Form("%s/%s_cat%d_interp.pdf",outDir.c_str(),proc_.c_str(),cat_));
  canv->Print(Form("%s/%s_cat%d_interp.png",outDir.c_str(),proc_.c_str(),cat_));
  delete canv;

}

RooSpline1D* FinalModelConstruction::graphToSpline(string name, TGraph *graph){
  
  vector<double> xValues, yValues;
  for (double mh=mhLow_; mh<(mhHigh_+0.25); mh+=0.5){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*MH,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

RooSpline1D* FinalModelConstruction::graphToSpline(string name, TGraph *graph, RooAbsReal *var){
  
  vector<double> xValues, yValues;
  for (double mh=mhLow_; mh<(mhHigh_+0.25); mh+=0.5){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*var,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

void FinalModelConstruction::getNormalization(){
 
  TGraph *temp = new TGraph();
  TF1 *pol2 = new TF1("pol","pol2",110,150);
  for (unsigned int i=0; i<allMH_.size(); i++){
    double mh = double(allMH_[i]);
    RooDataSet *data = stdDatasets[mh];
    // calcu eA as sumEntries / totalxs * totalbr * intL
    double effAcc = (data->sumEntries()/(intLumi->getVal()*norm->GetXsection(mh,proc_)*norm->GetBR(mh)));
    temp->SetPoint(i,mh,effAcc);
  }
  verbosity_ >=2 ?
    temp->Fit(pol2,"EMFEX0") :
    temp->Fit(pol2,"QEMFEX0")
  ;
  TGraph *eaGraph = new TGraph(pol2);
  RooSpline1D *eaSpline = graphToSpline(Form("fea_%s_cat%d_%dTeV",proc_.c_str(),cat_,sqrts_),eaGraph);
  RooSpline1D *xs = xsSplines[proc_];
	RooAbsReal *rateNuisTerm = getRateWithPhotonSyst(Form("rate_%s_cat%d_%dTeV",proc_.c_str(),cat_,sqrts_));
  
	finalNorm = new RooFormulaVar(Form("%s_norm",finalPdf->GetName()),Form("%s_norm",finalPdf->GetName()),"@0*@1*@2*@3",RooArgList(*xs,*brSpline,*eaSpline,*rateNuisTerm));
  // these are for plotting
  finalNormThisLum = new RooFormulaVar(Form("%s_normThisLumi",finalPdf->GetName()),Form("%s_normThisLumi",finalPdf->GetName()),"@0*@1*@2*@3*@4",RooArgList(*xs,*brSpline,*eaSpline,*rateNuisTerm,*intLumi));
  
	extendPdfRel = new RooExtendPdf(Form("extend%s",finalPdf->GetName()),Form("extend%s",finalPdf->GetName()),*finalPdf,*finalNorm);
  extendPdf = new RooExtendPdf(Form("extend%sThisLumi",finalPdf->GetName()),Form("extend%sThisLumi",finalPdf->GetName()),*finalPdf,*finalNormThisLum);
  // do secondary models
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    // sm higgs as bkg
    RooSpline1D *eaSpline_SM = graphToSpline(Form("fea_%s_cat%d_%dTeV_SM",proc_.c_str(),cat_,sqrts_),eaGraph,MH_SM);
    RooSpline1D *xs_SM = xsSplines_SM[proc_];
    finalNorm_SM = new RooFormulaVar(Form("%s_norm",finalPdf_SM->GetName()),Form("%s_norm",finalPdf_SM->GetName()),"@0*@1*@2*@3",RooArgList(*xs_SM,*brSpline_SM,*eaSpline_SM,*rateNuisTerm));
    // second degen higgs
    RooSpline1D *eaSpline_2 = graphToSpline(Form("fea_%s_cat%d_%dTeV_2",proc_.c_str(),cat_,sqrts_),eaGraph,MH_2);
    RooSpline1D *xs_2 = xsSplines_2[proc_];
    finalNorm_2 = new RooFormulaVar(Form("%s_norm",finalPdf_2->GetName()),Form("%s_norm",finalPdf_2->GetName()),"@0*@1*@2*@3",RooArgList(*xs_2,*brSpline_2,*eaSpline_2,*rateNuisTerm));
    // natural width
    RooSpline1D *eaSpline_NW = graphToSpline(Form("fea_%s_cat%d_%dTeV_NW",proc_.c_str(),cat_,sqrts_),eaGraph,MH);
    RooSpline1D *xs_NW = xsSplines_NW[proc_];
    finalNorm_NW = new RooFormulaVar(Form("%s_norm",finalPdf_NW->GetName()),Form("%s_norm",finalPdf_NW->GetName()),"@0*@1*@2*@3",RooArgList(*xs_NW,*brSpline_NW,*eaSpline_NW,*rateNuisTerm));
  }

}

void FinalModelConstruction::save(RooWorkspace *work){
  work->import(*finalPdf,RecycleConflictNodes());
  work->import(*finalNorm,RecycleConflictNodes());
  work->import(*finalNormThisLum,RecycleConflictNodes());
  work->import(*extendPdf,RecycleConflictNodes());
  work->import(*extendPdfRel,RecycleConflictNodes());
  for (map<int,RooDataSet*>::iterator it=stdDatasets.begin(); it!=stdDatasets.end(); it++){
    work->import(*(it->second));
  }

  // do secondary models
	if (doSecondaryModels){
		work->import(*finalPdf_SM,RecycleConflictNodes());
		work->import(*finalPdf_2,RecycleConflictNodes());
		work->import(*finalPdf_NW,RecycleConflictNodes());
		work->import(*finalNorm_SM,RecycleConflictNodes());
		work->import(*finalNorm_2,RecycleConflictNodes());
		work->import(*finalNorm_NW,RecycleConflictNodes());
	}
}
