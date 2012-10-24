#include "boost/lexical_cast.hpp"

#include "TMath.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "../interface/FMTTree.h"

using namespace std;

FMTTree::FMTTree(string infilename, string outfilename, string bdtname, string weightsFile, double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, int nIncCategories, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose):
	bdtname_(bdtname),
	crossCheck_(true),
 FMTBase(intLumi, is2011, mHMinimum, mHMaximum, mHStep, massMin, massMax, nDataBins, signalRegionWidth, sidebandWidth, numberOfSidebands, numberOfSidebandsForAlgos, numberOfSidebandGaps, massSidebandMin, massSidebandMax, nIncCategories, includeVBF, nVBFCategories, includeLEP, nLEPCategories, systematics, rederiveOptimizedBinEdges, AllBinEdges, verbose)
  {
    // open files and workspaces etc.
		inFile_ = TFile::Open(infilename.c_str());
		outFile_ = new TFile(outfilename.c_str(),"RECREATE");
    if (verbose_) cout << "Number of cats.. " << getNcats() << endl;
    if (verbose_) cout << "Init datasets..." << endl;
		ws_ = new RooWorkspace("cms_hgg_workspace");
		massVar_ = new RooRealVar("CMS_hgg_mass","CMS_hgg_mass",100.,180.);
		dataSet_ = new RooDataSet("data_mass","data_mass",RooArgList(*massVar_));
    for (int cat=0; cat<getNcats(); cat++){
      if (verbose_) cout << "\t\r" << cat << flush;
		  catDataSets_.push_back(new RooDataSet(Form("data_mass_cat%d",cat),Form("data_mass_cat%d",cat),RooArgList(*massVar_)));
    }
    if (verbose_) cout << endl;
    if (verbose_) cout << "Init tmva..." << endl;
    // tmva
		tmvaReader_ = new TMVA::Reader();
		tmvaReader_->AddVariable("bdtoutput",&diphotonBDT_);
		tmvaReader_->AddVariable("deltaMoverM",&deltaMOverM_);
		tmvaReader_->BookMVA(bdtname_.c_str(),weightsFile.c_str());

    if (verbose_) cout << "Init tree..." << endl;
		// set up vars
		initVariables();
		// initialize histograms
		vector<double> masses = getAllMH();
		vector<int> mcMasses = getMCMasses();
		vector<string> processes = getProcesses();
		
    if (verbose_) cout << "Init th1fs..." << endl;
		// data + bkg
		for (vector<double>::iterator mIt=masses.begin(); mIt!=masses.end(); mIt++){
			for (int cat=0; cat<getNcats(); cat++){
				th1fs_.insert(pair<string,TH1F*>(Form("th1f_data_BDT_grad_%3.1f_cat%d",*mIt,cat),new TH1F(Form("th1f_data_BDT_grad_%3.1f_cat%d",*mIt,cat),"BDT",5000,-1.,1.)));
				th1fs_.insert(pair<string,TH1F*>(Form("th1f_bkg_BDT_grad_%3.1f_cat%d",*mIt,cat),new TH1F(Form("th1f_bkg_BDT_grad_%3.1f_cat%d",*mIt,cat),"BDT",5000,-1.,1.)));
				for (int i=1; i<=getnumberOfSidebands(); i++){
					th1fs_.insert(pair<string,TH1F*>(Form("th1f_data_%dhigh_BDT_grad_%3.1f_cat%d",i,*mIt,cat),new TH1F(Form("th1f_data_%dhigh_BDT_grad_%3.1f_cat%d",i,*mIt,cat),"BDT",5000,-1.,1.)));
					th1fs_.insert(pair<string,TH1F*>(Form("th1f_data_%dlow_BDT_grad_%3.1f_cat%d",i,*mIt,cat),new TH1F(Form("th1f_data_%dlow_BDT_grad_%3.1f_cat%d",i,*mIt,cat),"BDT",5000,-1.,1.)));
					th1fs_.insert(pair<string,TH1F*>(Form("th1f_bkg_%dhigh_BDT_grad_%3.1f_cat%d",i,*mIt,cat),new TH1F(Form("th1f_bkg_%dhigh_BDT_grad_%3.1f_cat%d",i,*mIt,cat),"BDT",5000,-1.,1.)));
					th1fs_.insert(pair<string,TH1F*>(Form("th1f_bkg_%dlow_BDT_grad_%3.1f_cat%d",i,*mIt,cat),new TH1F(Form("th1f_bkg_%dlow_BDT_grad_%3.1f_cat%d",i,*mIt,cat),"BDT",5000,-1.,1.)));
				}
			}
		}
		// sig
		for (vector<int>::iterator mIt=mcMasses.begin(); mIt!=mcMasses.end(); mIt++){
			for (int cat=0; cat<getNcats(); cat++){
				for (vector<string>::iterator pIt=processes.begin(); pIt!=processes.end(); pIt++){
					th1fs_.insert(pair<string,TH1F*>(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d",pIt->c_str(),*mIt,cat), new TH1F(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d",pIt->c_str(),*mIt,cat),"BDT",5000,-1.,1.)));
					for (vector<string>::iterator sysIt=systematics.begin(); sysIt!=systematics.end(); sysIt++){
						th1fs_.insert(pair<string,TH1F*>(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d_%sDown01_sigma",pIt->c_str(),*mIt,cat,sysIt->c_str()), new TH1F(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d_%sDown01_sigma",pIt->c_str(),*mIt,cat,sysIt->c_str()),"BDT",5000,-1.,1.)));
						th1fs_.insert(pair<string,TH1F*>(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d_%sUp01_sigma",pIt->c_str(),*mIt,cat,sysIt->c_str()), new TH1F(Form("th1f_sig_BDT_grad_%s_%d.0_cat%d_%sUp01_sigma",pIt->c_str(),*mIt,cat,sysIt->c_str()),"BDT",5000,-1.,1.)));
					}
				}
			}
		}
}

FMTTree::~FMTTree(){
	
	// save histograms to file
	outFile_->cd();
  if (dataSet_->numEntries()>0){
    ws_->import(*dataSet_);
    cout << "All Data: --- " << dataSet_->numEntries() << endl;
    for (vector<RooDataSet*>::iterator dat=catDataSets_.begin(); dat!=catDataSets_.end(); dat++){
      cout << (*dat)->GetName() << " --- " << (*dat)->numEntries() << endl;
      ws_->import(**dat);
    }
    outFile_->cd();
    ws_->Write();
  }
  
  cout << "Writing histograms into file. Please be patient....." << endl;
  int i=0;
  int size=th1fs_.size();
	for (map<string,TH1F*>::iterator mapIt=th1fs_.begin(); mapIt!=th1fs_.end(); mapIt++){
    if (i%1000==0) cout << "\r" << i << "/" << size << flush;
		mapIt->second->Write();
		delete mapIt->second;
    i++;
	}
  cout << endl;
	outFile_->Close();
	inFile_->Close();

	delete inFile_;
	delete outFile_;
	delete tmvaReader_;
	delete massVar_;
	delete dataSet_;
	delete ws_;

}

void FMTTree::setdirname(string dir){
	dirname_=dir;
}

void FMTTree::addTreeToMap(map<string,TTree*>& theMap, string name, string label) {

	if (label=="0") label=name;
	
	TTree *temp = (TTree*)inFile_->Get(name.c_str());
	if (!temp) cerr << "WARNING -- tree " << name << " not found " << endl;
	else theMap.insert(pair<string,TTree*>(label,temp));

}

map<string,TTree*> FMTTree::getSignalTrees(string option){
  map<string,TTree*> result;
  vector<int> mcMasses = getMCMasses();
  vector<string> processes = getProcesses();
  for (vector<int>::iterator mSig=mcMasses.begin(); mSig!=mcMasses.end(); mSig++){
    for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++){
      if (option=="all") addTreeToMap(result,Form("%s/%s_m%d_8TeV",dirname_.c_str(),proc->c_str(),*mSig));
      else if (option==*proc) addTreeToMap(result,Form("%s/%s_m%d_8TeV",dirname_.c_str(),proc->c_str(),*mSig));
    }
  }
	//addTreeToMap(result,"full_mva_trees/ggh_m124_pu2012");
  return result;
}

map<string,TTree*> FMTTree::getDataTrees(){
  map<string,TTree*> result;
	addTreeToMap(result,Form("%s/Data",dirname_.c_str()));//,"data");
  return result;
}

map<string,TTree*> FMTTree::getBackgroundTrees(){

  map<string,TTree*> result;
	//addTreeToMap(result,Form("%s/qcd_30_8TeV_ff",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/qcd_30_8TeV_pf",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/qcd_30_8TeV_pp",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/qcd_40_8TeV_ff",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/qcd_40_8TeV_pf",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/qcd_40_8TeV_pp",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_20_8TeV_ff",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_20_8TeV_pf",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_20_8TeV_pp",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_40_8TeV_ff",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_40_8TeV_pf",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/gjet_40_8TeV_pp",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/diphojet_8TeV",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/dipho_Box_10_8TeV",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/dipho_Box_25_8TeV",dirname_.c_str()));//,"bkg");
	addTreeToMap(result,Form("%s/dipho_Box_250_8TeV",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/dipho_Born_10_8TeV",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/dipho_Born_25_8TeV",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/dipho_Born_250_8TeV",dirname_.c_str()));//,"bkg");
	//addTreeToMap(result,Form("%s/dyjetsll_50_8TeV",dirname_.c_str()));//,"bkg");
  return result;
}

void FMTTree::initVariables(){
  // initialise vector vars
  vector<string> systematics=getsystematics();
  for (unsigned int s=0; s<systematics.size(); s++){
    massSyst_.push_back(pair<double,double>(0.,0.));
    bdtoutputSyst_.push_back(pair<double,double>(0.,0.));
    weightSyst_.push_back(pair<double,double>(0.,0.));
    categorySyst_.push_back(pair<int,int>(0.,0.));
  }

}

void FMTTree::setBranchVariables(TTree *tree){
  
  tree->SetBranchAddress("mass",&mass_);
  tree->SetBranchAddress("bdtoutput",&bdtoutput_);
  tree->SetBranchAddress("weight",&weight_);
  tree->SetBranchAddress("category",&category_);
  vector<string> systs = getsystematics();
  for (unsigned int s=0; s<systs.size(); s++){
    tree->SetBranchAddress(Form("mass_%s_Down",systs[s].c_str()),&massSyst_[s].first);
    tree->SetBranchAddress(Form("mass_%s_Up",systs[s].c_str()),&massSyst_[s].second);
    tree->SetBranchAddress(Form("bdtoutput_%s_Down",systs[s].c_str()),&bdtoutputSyst_[s].first);
    tree->SetBranchAddress(Form("bdtoutput_%s_Up",systs[s].c_str()),&bdtoutputSyst_[s].second);
    tree->SetBranchAddress(Form("weight_%s_Down",systs[s].c_str()),&weightSyst_[s].first);
    tree->SetBranchAddress(Form("weight_%s_Up",systs[s].c_str()),&weightSyst_[s].second);
    tree->SetBranchAddress(Form("category_%s_Down",systs[s].c_str()),&categorySyst_[s].first);
    tree->SetBranchAddress(Form("category_%s_Up",systs[s].c_str()),&categorySyst_[s].second);
  }

}

float FMTTree::tmvaGetVal(float dMoM, float bdt){
	deltaMOverM_ = dMoM;
	diphotonBDT_ = bdt;
	return tmvaReader_->EvaluateMVA(bdtname_.c_str());
}

int FMTTree::icCat(int cat){
  if (cat<0) return -1;
  else if (cat>=0 && cat<4) return 0;
  else if (cat==4) return 2;
  else if (cat==5) return 1;
  else return cat-3;
}

void FMTTree::FillHist(string type, int sideband, double mh){
  int cat = icCat(category_);
  if (cat<0) return;
  if (sideband==0) {
    float val = tmvaGetVal((mass_-mh)/mh,bdtoutput_);
    th1fs_[Form("th1f_%s_BDT_grad_%5.1f_cat%d",type.c_str(),mh,cat)]->Fill(val,weight_);
  }
  else if (sideband<0) {
    double hypothesisModifier = (1.-sidebandWidth_)/(1.+sidebandWidth_);
    double evalMH = (mh*(1.-signalRegionWidth_)/(1+sidebandWidth_))*(TMath::Power(hypothesisModifier,(-1*sideband)-1));
    float val = tmvaGetVal((mass_-evalMH)/evalMH,bdtoutput_);
    th1fs_[Form("th1f_%s_%dlow_BDT_grad_%5.1f_cat%d",type.c_str(),sideband*-1,mh,cat)]->Fill(val,weight_);
  }
  else if (sideband>0) {
    double hypothesisModifier = (1.+sidebandWidth_)/(1.-sidebandWidth_);
    double evalMH = (mh*(1.+signalRegionWidth_)/(1-sidebandWidth_))*(TMath::Power(hypothesisModifier,sideband-1));
    float val = tmvaGetVal((mass_-evalMH)/evalMH,bdtoutput_);
    th1fs_[Form("th1f_%s_%dhigh_BDT_grad_%5.1f_cat%d",type.c_str(),sideband,mh,cat)]->Fill(val,weight_);
  }
}

void FMTTree::FillSigHist(string proc, double mh){
  int cat = icCat(category_);
  if (cat<0) return;
	if (mass_>=(1.-sidebandWidth_)*mh && mass_<=(1.+sidebandWidth_)*mh){
		float val = tmvaGetVal((mass_-mh)/mh,bdtoutput_);
		th1fs_[Form("th1f_sig_BDT_grad_%s_%5.1f_cat%d",proc.c_str(),mh,cat)]->Fill(val,weight_);
	}
}

void FMTTree::FillSystHist(string proc, double mh){
	double cutLow = (1.-sidebandWidth_)*mh;
	double cutHigh = (1.+sidebandWidth_)*mh;
	
	for (unsigned int s=0; s<systematics_.size(); s++){
		int cat;
		float mass, val, weight;
		string shift[2] = {"Down","Up"};
		for (int t=0; t<2; t++){
			if (t==0){
				cat = icCat(categorySyst_[s].first);
				mass = massSyst_[s].first;
				val = tmvaGetVal((mass-mh)/mh,bdtoutputSyst_[s].first);
				weight = weightSyst_[s].first;
			}
			else {
				cat = icCat(categorySyst_[s].second);
				mass = massSyst_[s].second;
				val = tmvaGetVal((mass-mh)/mh,bdtoutputSyst_[s].second);
				weight = weightSyst_[s].second;
			}
      if (cat<0) return;
			if (mass>=cutLow && mass<=cutHigh) {
				th1fs_[Form("th1f_sig_BDT_grad_%s_%5.1f_cat%d_%s%s01_sigma",proc.c_str(),mh,cat,systematics_[s].c_str(),shift[t].c_str())]->Fill(val,weight);
			}
		}
	}
}

string FMTTree::getProc(string name){
  vector<string> processes=getProcesses();
  for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++){
    if (name.find(*proc)!=string::npos) return *proc;
    else continue;
  }
}

int FMTTree::getMH(string name){
	int m = boost::lexical_cast<int>(name.substr(name.find_last_of("m")+1,3));
  return m; 
}

void FMTTree::FillMassDatasets(){
  massVar_->setVal(mass_);
	dataSet_->add(RooArgList(*massVar_));
  catDataSets_[icCat(category_)]->add(RooArgList(*massVar_));
}

void FMTTree::doCrossCheck(vector<pair<int,map<string,TTree*> > > allTrees, int mH){
	
	TH1F *sig = new TH1F("sig125","sig125",100,-1,1);
	TH1F *bkg = new TH1F("bkg","bkg",100,-1,1);
	TH1F *data = new TH1F("data","data",100,-1,1);
	vector<pair<TH1F*,TH1F*> > syst;
	for (unsigned int s=0; s<systematics_.size(); s++){
		TH1F *tempDown = new TH1F(Form("%sDown",systematics_[s].c_str()),Form("%sDown",systematics_[s].c_str()),100,-1,1);
		TH1F *tempUp = new TH1F(Form("%sUp",systematics_[s].c_str()),Form("%sUp",systematics_[s].c_str()),100,-1,1);
		syst.push_back(pair<TH1F*,TH1F*>(tempDown,tempUp));
	}

	gROOT->SetBatch();
	for (unsigned int it=0; it<allTrees.size(); it++) {
    int type = allTrees[it].first;
    map<string,TTree*> treeMap = allTrees[it].second;
		for (map<string,TTree*>::iterator treeIt=treeMap.begin(); treeIt!=treeMap.end(); treeIt++){
			if (type==0){
				treeIt->second->Draw("bdtoutput>>+data","weight");
			}
			else if (type>0){
				treeIt->second->Draw("bdtoutput>>+bkg","weight");
			}
			else if (type<0){
				if (treeIt->first.find("m125")!=string::npos) {
					treeIt->second->Draw("bdtoutput>>+sig125","weight");
					for (vector<string>::iterator sIt=systematics_.begin(); sIt!=systematics_.end(); sIt++){
						treeIt->second->Draw(Form("bdtoutput_%s_Down>>+%sDown",sIt->c_str(),sIt->c_str()),"weight");
						treeIt->second->Draw(Form("bdtoutput_%s_Up>>+%sUp",sIt->c_str(),sIt->c_str()),"weight");
					}
				}
			}
		}
	}
	sig->SetLineColor(kRed);
	sig->SetFillColor(kRed);
	sig->SetFillStyle(3004);
	bkg->SetLineColor(kBlue-9);
	bkg->SetFillColor(kBlue-9);
	data->SetMarkerStyle(20);
	gStyle->SetOptStat(0);
	TCanvas *canv = new TCanvas();
	bkg->GetYaxis()->SetRangeUser(0,bkg->GetMaximum()*1.2);
	bkg->Draw();
	sig->Draw("same");
	data->Draw("lep same");
	canv->Print("plots/pdf/crossCheck.pdf");
	canv->Print("plots/png/crossCheck.png");
	for (unsigned int s=0; s<systematics_.size(); s++){
		sig->SetLineColor(kBlack);
		sig->SetLineWidth(2);
		syst[s].first->SetLineColor(kRed);
		syst[s].second->SetLineColor(kBlue);
		sig->GetYaxis()->SetRangeUser(0,sig->GetMaximum()*1.2);
		sig->Draw();
		syst[s].first->Draw("same");
		syst[s].second->Draw("same");
		canv->Print(Form("plots/pdf/ccheck_%s.pdf",systematics_[s].c_str()));
		canv->Print(Form("plots/png/ccheck_%s.pdf",systematics_[s].c_str()));
	}
	delete canv;
}

void FMTTree::printTrees(vector<pair<int,map<string,TTree*> > > allTrees){
	cout << "Available trees: " << endl;
	for (unsigned int it=0; it<allTrees.size(); it++) {
    int type = allTrees[it].first;
    map<string,TTree*> treeMap = allTrees[it].second;
		for (map<string,TTree*>::iterator mapIt = treeMap.begin(); mapIt!=treeMap.end(); mapIt++){
			cout << Form("%20s --  %20s",mapIt->first.c_str(),mapIt->second->GetName()) << endl;
		}
	}
}

void FMTTree::run(string option){
 	
	setdirname("full_mva_trees");
  map<string,TTree*> sigTrees = getSignalTrees(option);
  map<string,TTree*> bkgTrees = getBackgroundTrees();
  map<string,TTree*> dataTrees = getDataTrees();

  vector<pair<int,map<string,TTree*> > > allTrees;
  if (option=="data" || option=="all") allTrees.push_back(pair<int,map<string,TTree*> >(0,dataTrees));
  if (option=="bkg" || option=="all") allTrees.push_back(pair<int,map<string,TTree*> >(1,bkgTrees));
  if (option=="ggh" || option=="vbf" || option=="wzh" || option=="tth" || option=="all") allTrees.push_back(pair<int,map<string,TTree*> >(-1,sigTrees));
  
  

	printTrees(allTrees);
	if (crossCheck_) doCrossCheck(allTrees);
  
  vector<double> masses = getAllMH();
  vector<string> systematics = getsystematics();

	// loop over trees
  int number = sigTrees.size() + bkgTrees.size() + dataTrees.size();
  cout << number << " trees in total" << endl;
	for (unsigned int it=0; it<allTrees.size(); it++) {
  //for (vector<pair<int,map<string,TTree*> > >::iterator typeIt=allTrees.begin(); typeIt!=allTrees.end(); typeIt++){
    int type = allTrees[it].first;
    map<string,TTree*> treeMap = allTrees[it].second;
    for (map<string,TTree*>::iterator mapIt=treeMap.begin(); mapIt!=treeMap.end(); mapIt++){
      setBranchVariables(mapIt->second);
      string thisProc;
      int thisMH;
      if (type<0) {
        thisProc = getProc(mapIt->first);
        thisMH = boost::lexical_cast<double>(getMH(mapIt->first));
      }
      cout << "Tree - " << mapIt->first << " -- " << mapIt->second->GetEntries() << endl;
      TStopwatch sw;
      sw.Start();
      for (int entry=0; entry<(mapIt->second)->GetEntries(); entry++){
        (mapIt->second)->GetEntry(entry);
        if (entry%1000==0) cout << "\r" << entry << "/" << mapIt->second->GetEntries() << flush;
        // Data and Bkg
				if (type==0){
					FillMassDatasets();
				}
        if (type>=0){
          for (vector<double>::iterator mIt=masses.begin(); mIt!=masses.end(); mIt++){
            //cout << "MH: " << *mIt << endl;
            // in signal region
            //cout << "si: " << (1.-sidebandWidth_)*(*mIt) << "  " << (1.+sidebandWidth_)*(*mIt) << endl;
            if (mass_>=(1.-sidebandWidth_)*(*mIt) && mass_<=(1.+sidebandWidth_)*(*mIt)){
              if (type==0) FillHist("data",0,*mIt);
              if (type>0) FillHist("bkg",0,*mIt);
              //cout << " \t in signal region " << *mIt << endl;
            }
            // in sidebands
            int nL = (getNsidebandsUandD(*mIt)).first;
            int nH = (getNsidebandsUandD(*mIt)).second;
            vector<double> lEdge = getLowerSidebandEdges(*mIt);
            vector<double> hEdge = getUpperSidebandEdges(*mIt);
            
            for (int l=0; l<nL; l++){
              //cout << "l"  << l << ": " << lEdge[l+1] << "  " << lEdge[l] << endl;
              if (mass_>=lEdge[l+1] && mass_<=lEdge[l]){
                if (type==0) FillHist("data",-1*(l+1+getnumberOfSidebandGaps()),*mIt);
                if (type>0) FillHist("bkg",-1*(l+1+getnumberOfSidebandGaps()),*mIt);
                //cout << " \t in sideband " << -1*(l+1+getnumberOfSidebandGaps()) << " " << *mIt << endl;
              }
            }
            for (int h=0; h<nH; h++){
              //cout << "h" << h << ": " << hEdge[h] << "  " << hEdge[h+1] << endl;
              if (mass_>=hEdge[h] && mass_<=hEdge[h+1]){
                if (type==0) FillHist("data",(h+1+getnumberOfSidebandGaps()),*mIt);
                if (type>0) FillHist("bkg",(h+1+getnumberOfSidebandGaps()),*mIt);
                //cout << " \t in sideband " << h+1+getnumberOfSidebandGaps() << " " << *mIt << endl;
              }
            }
          }
        }
        // Signal
        else {
					FillSigHist(thisProc,thisMH);
					FillSystHist(thisProc,thisMH);
        }
      }
      sw.Stop();
      cout << endl;
      cout << "\t Took -- ";
      sw.Print();
    }
  }


}
