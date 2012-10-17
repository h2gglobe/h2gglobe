#include "../interface/FMTTree.h"

using namespace std;

FMTTree::FMTTree(double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, int nIncCategories, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose):
  FMTBase(double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, int nIncCategories, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose)
  {
    initVariables();
		outFile_ = new TFile("CMS-HGG_fmt_trees.root","RECREATE");
		ws_ = new RooWorkspace("cms_hgg_workspace");
		massVar_ = new RooRealVar("CMS_hgg_mass","CMS_hgg_mass",100.,180.);
		dataSet_ = new RooDataSet("data_mass","data_mass",RooArgList(*massVar_));

		// initialize histograms
		vector<double> masses = getMHMasses();
		vector<int> mcMasses = getMCMasses();
		vector<string> systematics = getsystematics();
		vector<string> processes = getProcesses();
		
		// data + bkg
		for (vector<double>::iterator mIt=masses.begin(); mIt!=masses.end(); mIt++){
			for (int cat=0; cat<getNcats(); cat++){
				th1fs_.insert(pair<string,TH1F*>(Form("th1f_data_BDT_grad_%3.1f_cat%d",*mIt,cat),new TH1F(Form("th1f_data_BDT_grad_%3.1f_cat%d",*mIt,cat),"BDT",5000,-1.,1.)));
				th1fs_.insert(pair<string,TH1F*>(Form("th1f_bkg_BDT_grad_%3.1f_cat%d",*mIt,cat),new TH1F(Form("th1f_bkg_BDT_grad_%3.1f_cat%d",*mIt,cat),"BDT",5000,-1.,1.)));
				for (int i=1; i<=getnumberOfSidebandsForAlgos(); i++){
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
	ws_->Import(*dataSet_);
	ws_->Write();
	for (map<string,TH1F*>::iterator mapIt=th1fs_.begin(); mapIt!=th1fs_.end(); mapIt++){
		mapIt->second->Write();
		if (verbose_) cout << "Writing " << mapIt->first << " to file" << endl;
		delete mapIt->second;
	}
	outFile_->Close();
	inFile_->Close();

	delete inFile;
	delete outFile;
	delete massVar_;
	delete dataSet_;
	delete ws_;

}

map<string,TTree*> FMTTree::getSignalTrees(TFile *file, string dir){
  map<string,TTree*> result;
  vector<int> mcMasses = getMCMasses();
  vector<string> processes = getProcesses();
  for (vector<int>::iterator mSig=mcMasses.begin(); mSig!=mcMasses.end(); mSig++){
    for (vector<string>::iterator proc=processes.begin(); proc!=processes.end(); proc++){
      result.insert(pair<string,TTree*>(Form("sig_%s_m%d",proc->c_str(),*mSig),(TTree*)file->Get(Form("%s/%s_m%d_pu2012",dir.c_str(),proc->c_str(),*mSig))));
    }
  }
  return result;
}

map<string,TTree*> FMTTree::getDataTrees(TFile *file, string dir){
  map<string,TTree*> result;
  result.insert(pair<string,TTree*>("data",(TTree*)file->Get(Form("%s/Data",dir.c_str()))));
  return result;
}

map<string,TTree*> FMTTree::getBackgroundTrees(TFile *file, string dir){

  map<string,TTree*> result;
  
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_30_8TeV_ff",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_30_8TeV_pf",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_30_8TeV_pp",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_40_8TeV_ff",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_40_8TeV_pf",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/qcd_40_8TeV_pp",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_20_8TeV_ff",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_20_8TeV_pf",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_20_8TeV_pp",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_40_8TeV_ff",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_40_8TeV_pf",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/gjet_40_8TeV_pp",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/diphojet_8TeV",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Box_10_8TeV",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Box_25_8TeV",dir.c_str()))));
  result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Box_250_8TeV",dir.c_str()))));
  //result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Born_10_8TeV",dir.c_str()))));
  //result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Born_25_8TeV",dir.c_str()))));
  //result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dipho_Born_250_8TeV",dir.c_str()))));
  //result.insert(pair<string,TTree*>("bkg",(TTree*)file->Get(Form("%s/dyjetsll_50_8TeV",dir.c_str()))));

  return result;
}

void FMTTree::initVariables(){
  // initialise vector vars
  vector<string> systematics=getsystematics();
  for (unsigned int s=0; s<systematics.size(); s++){
    massSyst_.push_back(pair<float,float>(0.,0.));
    bdtouputSyst_.push_back(pair<float,float>(0.,0.));
    weightSyst_.push_back(pair<float,float>(0.,0.));
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

int FMTTree::icCat(int cat){
  if (cat>=0 || cat<getnIncCategories()) return 0;
  else if (cat==getnIncCategories()) return 2;
  else if (cat==getnIncCategories+1) return 1;
  else return cat+1-getnIncCategories();
}

void FMTTree::FillHist(string type, int sideband, double mh){
  int cat = icCat(category_);
  float val = tmvaGetVal((mass_-mh)/mh,bdtoutput_);
  if (sideband==0) th1fs_[Form("th1f_%s_BDT_grad_%5.1f_cat%d",type.c_str(),mh,cat)]->Fill(val,weight_);
  else if (sideband<0) th1fs_[Form("th1f_%s_%dlow_BDT_grad_%5.1f_cat%d",type.c_str(),sideband*-1,mh,cat)]->Fill(val,weight_);
  else if (sideband>0) th1fs_[Form("th1f_%s_%dhigh_BDT_grad_%5.1f_cat%d",type.c_str(),sideband,mh,cat)]->Fill(val,weight_);
}

void FMTTree::FillSigHist(string proc, double mh){
  int cat = icCat(category_);
	if (mass_>=(1.-getsidebandWidth())*mh && mass_<=(1.+getsidebandWidth())){
		float val = tmvaGetVal((mass_-mh)/mh,bdtoutput_);
		th1fs_[Form("th1f_sig_BDT_grad_%s_%5.1f_cat%d",proc.c_str(),mh,cat)]->Fill(val,weight_);
	}
}

void FMTTree::FillSystHist(string proc, double mh){
	vector<string> systs = getsystematics();
	double cutLow = (1.-getsidebandWidth())*mH;
	double cutHigh = (1.+getsidebandWidth())*mH;
	
	for (unsigned int s=0; s<systematics.size(); s++){
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
			if (mass>=cutLow && mass<=cutHigh) {
				th1fs_[Form("th1f_sig_BDT_grad_%s_%5.1f_cat%d_%s%s01_sigma",proc.c_str(),mh,cat,systematics[s].c_str(),shift[t].c_str())]->Fill(val,weight);
			}
		}
	}
}

string FMTTree::getProc(string name){
 return name.substr(name.find("h")-2,3);
}

int FMTTree::getMH(string name){
  return boost::lexical_cast<int>(name.substr(name.find("m")+1,3)); 
}

void FMTTree::FillMassDatasets(){
	massVar_->setVal(mass_);
	dataSet_.add(RooArgList(*massVar));
}

void FMTTree::run(string filename){
  
  inFile_ = TFile::Open(filename.c_str());
  map<string,TTree*> sigTrees = getSignalTrees(file,"full_mva_trees");
  map<string,TTree*> bkgTrees = getBackgroundTrees(file,"full_mva_trees");
  map<string,TTree*> dataTrees = getDataTrees(file,"full_mva_trees");
  
  vector<pair<int,map<string,TTree*> > > allTrees;
  allTrees.push_back(pair<int,map<string,TTree*> >(-1,sigTrees));
  allTrees.push_back(pair<int,map<string,TTree*> >(0,dataTrees));
  allTrees.push_back(pair<int,map<string,TTree*> >(1,bkgTrees));
  
  vector<double> masses = getMHMasses();
  vector<string> systematics = getsystematics();

	// loop over trees
  for (vector<pair<int,map<string,TTree*> > >::iterator typeIt=allTrees.begin(); typeIt=allTrees.end(); typeIt++){
    int type = typeIt->first;
    map<string,TTree*> map = typeIt->second;
    for (map<string,TTree*>::iterator mapIt=map.begin(); mapIt!=map.end(); mapIt++){
      setBranchVariables(mapIt->second);
      for (int entry=0; entry<(mapIt->second)->GetEntries(); e++){
        (mapIt->second)->GetEntry(entry);
        // Data and Bkg
				if (type==0){
					FillMassDatasets();
				}
        if (type>=0){
          for (vector<double>::iterator mIt=masses.begin(); mIt=masses.end(); mIt++){
            // in signal region
            if (mass_>=(1.-sidebandWidth_)*(*mIt) && mass_<=(1.+sidebandWidth_)*(*mIt)){
              FillHist(mapIt->first,0,*mIt);
            }
            // in sidebands
            int nL = (getNsidebandsUandD(*mIt)).first;
            int nH = (getNsidebandsUandD(*mIt)).second;
            vector<double> lEdge = getLowerSidebandEdges(*mIt);
            vector<double> hEdge = getUpperSidebandEdges(*mIt);

            for (int l=0; l<nL; l++){
              if (mass_>=lEdge[l+1] && mass_<=lEdge[l]){
                FillHist(mapIt->first,-1*(l+1+getnumberOfSidebandGaps()),*mIt);
              }
            }
            for (int h=0; h<nH; h++){
              if (mass_>=hEdge[h] && mass_<=hEdge[h+1]){
                FillHist(mapIt->first,(l+1+getnumberOfSidebandGaps()),*mIt);
              }
            }
          }
        }
        // Signal
        else {
          string proc = getProc(mapIt->first);
          double mH = boost::lexical_cast<double>(getMH(mapIt->first));
					FillSigHist(proc,mH);
					FillSystHist(proc,mH);
        }
      }
    }
  }


}
