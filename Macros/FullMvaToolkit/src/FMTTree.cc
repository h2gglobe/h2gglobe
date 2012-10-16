#include "../interface/FMTTree.h"

using namespace std;

FMTTree::FMTTree(double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, int nIncCategories, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose):
  FMTBase(double intLumi, bool is2011, int mHMinimum, int mHMaximum, double mHStep, double massMin, double massMax, int nDataBins, double signalRegionWidth, double sidebandWidth, int numberOfSidebands, int numberOfSidebandsForAlgos, int numberOfSidebandGaps, double massSidebandMin, double massSidebandMax, int nIncCategories, bool includeVBF, int nVBFCategories, bool includeLEP, int nLEPCategories, vector<string> systematics, bool rederiveOptimizedBinEdges, vector<map<int, vector<double> > > AllBinEdges, bool verbose)
  {
    initVariables();
  }

FMTTree::~FMTTree(){
}

map<string,TTree*> FMTTree::getSignalTrees(TFile *file, string dir=""){
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

map<string,TTree*> FMTTree::getDataTrees(TFile *file, string dir=""){
  map<string,TTree*> result;
  result.insert(pair<string,TTree*>("data",(TTree*)file->Get(Form("%s/Data",dir.c_str()))));
  return result;
}

map<string,TTree*> FMTTree::getBackgroundTrees(TFile *file, string dir=""){

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

void FMTTree::icCat(int cat){
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

string FMTTree::getProc(string name){
 return name.substr(name.find("h")-2,3);
}

int FMTTree::getMH(string name){
  return boost::lexical_cast<int>(name.substr(name.find("m")+1,3)); 
}

void FMTTree::run(string filename){
  
  TFile *file = TFile::Open(filename.c_str());
  map<string,TTree*> sigTrees = getSignalTrees(file,"full_mva_trees");
  map<string,TTree*> bkgTrees = getBackgroundTrees(file,"full_mva_trees");
  map<string,TTree*> dataTrees = getDataTrees(file,"full_mva_trees");
  
  vector<pair<int,map<string,TTree*> > > allTrees;
  allTrees.push_back(pair<int,map<string,TTree*> >(-1,sigTrees));
  allTrees.push_back(pair<int,map<string,TTree*> >(0,dataTrees));
  allTrees.push_back(pair<int,map<string,TTree*> >(1,bkgTrees));
  
  vector<double> masses = getMHMasses();
  vector<double> mcMasses = getMCMasses();
  // loop over trees
  for (vector<pair<int,map<string,TTree*> > >::iterator typeIt=allTrees.begin(); typeIt=allTrees.end(); typeIt++){
    int type = typeIt->first;
    map<string,TTree*> map = typeIt->second;
    for (map<string,TTree*>::iterator mapIt=map.begin(); mapIt!=map.end(); mapIt++){
      setBranchVariables(mapIt->second);
      for (int entry=0; entry<(mapIt->second)->GetEntries(); e++){
        (mapIt->second)->GetEntry(entry);
        // Data and Bkg
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
          int mH = getMH(mapIt->first);
          // in signal region
          if (mass_>=(1.-sidebandWidth_)*(*mIt) && mass_<=(1.+sidebandWidth_)*(*mIt)){
            FillHist(mapIt->first,0,*mIt);
          }
        }
      }
      
    }

  }


}
