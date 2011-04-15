#define LDEBUG 0 

#include "LoopAll.h"

#include "TGraph.h"
#include <iostream>
#include <math.h>

#include "Util.h"
#include "TRandom.h"

#include <TLorentzVector.h>
#include <TVector3.h>

#include "stdlib.h"

using namespace std;

#include "Tools.h"


#include "TH2.h"
#include "TStyle.h" 
#include "TCanvas.h"

#include "GeneralFunctions_cc.h"

#include "PhotonAnalysis/PhotonAnalysisReducedOutputTree.h"
#include "PhotonAnalysis/PhotonAnalysisFunctions_cc.h"

LoopAll::LoopAll(TTree *tree) {
  
  rooContainer = new RooContainer();

}

LoopAll::~LoopAll() {
  if (!fChain)
    return;

  delete fChain->GetCurrentFile();
}
          


Int_t LoopAll::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


Long64_t LoopAll::LoadTree(Long64_t entry) {

  // Set the environment to read one entry
  if (!fChain) return -5;
  Int_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->IsA() != TChain::Class()) 
    return centry;
  TChain *chain = (TChain*)fChain;
  Notify();
  
  return centry;
}

Bool_t LoopAll::Notify() {
  return kTRUE;
}

void LoopAll::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}


void LoopAll::Loop(Double_t a) {
  
  if (fChain == 0) 
    return;
  
  Int_t nentries = Int_t(fChain->GetEntriesFast());
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) 
      break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
  }
}

void LoopAll::InitHistos(){

  for(int ind=0; ind<sampleContainer.size(); ind++) {
    HistoContainer temp(ind);
    histoContainer.push_back(temp);
  }

}

void LoopAll::InitReal(Int_t typerunpass) {

  // Set branch addresses
  typerun=typerunpass;
  hfile = (TFile*)gROOT->FindObject(utilInstance->histFileName); 
  if (hfile) 
    hfile->Close();

  //for(int ind=0; ind<sampleContainer.size(); ind++) {
  // HistoContainer temp(ind);
  //  histoContainer.push_back(temp);
 // }

  if(LDEBUG) cout << "doing InitRealPhotonAnalysis" << endl;
  InitRealPhotonAnalysis(typerun);
  if(LDEBUG) cout << "finished InitRealPhotonAnalysis" << endl;

  if (utilInstance->makeOutputTree) 
    utilInstance->outputFile->cd();
  cout<< "LoopAll::InitReal END" <<endl;
  
}

void LoopAll::TermReal(Int_t typerunpass) {

  typerun=typerunpass;

  TermRealPhotonAnalysis(typerun);
  
  if (utilInstance->makeOutputTree){ 
    utilInstance->outputFile->cd();
    utilInstance->outputParReductions++;
    utilInstance->outputTreePar->Fill();
    utilInstance->outputTreePar->Write();
  }
}


void LoopAll::Init(Int_t typerunpass, TTree *tree) {
  
  // Set branch addresses
  typerun=typerunpass;
  
  fChain = tree;
  if (tree == 0) 
    return;

  fChain->SetMakeClass(1);

  mySetBranchAddressRedPhotonAnalysis();

  Notify();
}


void LoopAll::setUtilInstance(Util *ut) {
  utilInstance = ut;
}

void LoopAll::Loop(Int_t a) {
  
  makeOutputTree = utilInstance->makeOutputTree;
  if (makeOutputTree) {
    outputTree = utilInstance->outputTree;
    outputFile = utilInstance->outputFile;
    if(outputFile) 
      outputFile->cd();
  }

  Int_t nentries = 0;
  if(fChain) 
    nentries = Int_t(fChain->GetEntriesFast());

  Int_t nbytes = 0, nb = 0;

  outputEvents=0;

  int hasoutputfile=0;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    
    if(jentry%10000==0) 
      cout << "Entry: "<<jentry<<endl;
    
    if(LDEBUG) 
      cout<<"call LoadTree"<<endl;
    
    Int_t ientry = LoadTree(jentry);
  
    if (ientry < 0) 
      break;
   
      
    if(utilInstance->typerun == 1) {


        //nb = fChain->GetEntry(jentry,1);
	
      
    } else {
      //HERE NEED TO DO A GLOBAL GETENTRY
      nb=0;
    }

    nbytes += nb;

    if(LDEBUG) 
      cout<<"Call FillandReduce "<<endl;
      
    hasoutputfile = FillAndReduce(jentry);
    if(LDEBUG) 
      cout<<"Called FillandReduce "<<endl;
  }

  if(hasoutputfile) {
    if(outputFile) {
      outputFile->cd();
      if (LDEBUG)
	cout<<"LoopAll_cc writing outputTree"<<endl;
      outputTree->Write(0,TObject::kWriteDelete);
    }

    
    if(outputFile) {
      outputFile->cd();
      if(utilInstance->TreesPar[a]) {

	std::vector<std::string> *parameters = new std::vector<std::string>;
	std::string *job_maker = new std::string;
	Int_t tot_events, sel_events, type, version, reductions;
	Int_t red_events[20];
	
	utilInstance->TreesPar[a]->SetBranchAddress("tot_events", &tot_events);
	utilInstance->TreesPar[a]->SetBranchAddress("sel_events", &sel_events);
	utilInstance->TreesPar[a]->SetBranchAddress("type", &type);
	utilInstance->TreesPar[a]->SetBranchAddress("version", &version);
	utilInstance->TreesPar[a]->SetBranchAddress("parameters", &parameters);
	utilInstance->TreesPar[a]->SetBranchAddress("job_maker", &job_maker);
	if (utilInstance->TreesPar[a]->FindBranch("reductions")) {
	  utilInstance->TreesPar[a]->SetBranchAddress("reductions", &reductions);
	  utilInstance->TreesPar[a]->SetBranchAddress("red_events", &red_events);
	}

	utilInstance->TreesPar[a]->GetEntry(0);

	if (a == 0) {
	  utilInstance->outputParTot_Events = tot_events;
	  utilInstance->outputParSel_Events = sel_events;
	} else {
	  utilInstance->outputParTot_Events += tot_events;
	  utilInstance->outputParSel_Events += sel_events;
	}

	utilInstance->outputParType = type;
	utilInstance->outputParVersion = version;
	utilInstance->outputParParameters = &(*parameters);
	utilInstance->outputParJobMaker = job_maker;

	if (!utilInstance->TreesPar[a]->FindBranch("reductions")) {
	  utilInstance->outputParReductions = 0;
	  for (int i=0; i<20; i++) {
	    utilInstance->outputParRed_Events[i] = -1;
	  }
	  utilInstance->outputParRed_Events[0] += (int)countersred[1];
	} else {
	  utilInstance->outputParReductions = reductions;
	  utilInstance->outputParRed_Events[reductions] += (int)countersred[1];
	}
      } else {
	std::cerr << "Cannot write Parameter tree." << std::endl;
      }
    }
  }
  
  int oldnentries=nentries;
  if (nentries == utilInstance->sel_events) {
    nentries = utilInstance->tot_events;
  }

  if(countersred[1] || oldnentries==0) {
    printf("red: %d_%d \n",(int)countersred[0], (int) countersred[1]);
  } else { 
    printf("norm: %d \n",(int)counters[0]);
  }
}

void LoopAll::myWriteFits() {

  
  hfile = new TFile(utilInstance->histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  hfile->cd();
  rooContainer->Save();
}
void LoopAll::myWritePlot() {

  hfile = new TFile(utilInstance->histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  hfile->cd();
  for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
    histoContainer[ind].Save();
  }
  if (utilInstance->makeOutputTree) 
    utilInstance->outputFile->cd();
}

void LoopAll::myWriteCounters() {
  if(LDEBUG) std::cout<<"LoopAll::myWriteCounters - START"<<std::endl;
  
  const int samples = sampleContainer.size();
   
  if(LDEBUG) std::cout<<"samples = "<<samples<<std::endl;
  
  stringstream fileLinesInfo[samples];
  stringstream fileLinesCatFile[samples];
  stringstream fileLinesCatCounts[samples];
  stringstream fileLinesCatEvents[samples];
  stringstream fileLinesCatSigma[samples];
  stringstream fileLinesCatEff[samples][3];
  
  if(LDEBUG) std::cout<<"initialize streams"<<std::endl;

  double counterNumerator_tot[samples];
  double counterDenominator_tot[samples][3];
  double counterEvents_tot[samples];
  double counterSigma_tot[samples];
    
  if (counterContainer.size() < 1)
    return;

  TString msName = utilInstance->histFileName.ReplaceAll(".root", 5, ".csv", 4);
  FILE *file;
  TString s2 = msName+".csv";
  file = fopen(s2, "w");

  if(LDEBUG) std::cout<<"msName is "<< msName <<std::endl;
  if(LDEBUG) std::cout<<"s2 is "<< s2 <<std::endl;
  
  for(unsigned int i=0; i<counterContainer[0].mapSize(); i++) {
    if(LDEBUG) std::cout<<"counterContainer[0].mapSize() is "<< counterContainer[0].mapSize() <<std::endl;
    fprintf(file, "Number, Sample, Counter Name, Counted, TotalEvents, Sigma, Selection 1, Eff, Selection 2, Eff, Selection 3, Eff");
    fprintf(file, ",indexfiles, namefile, scale, lumi, intlum, weight");
    fprintf(file, "\n");
    
    int ncat = counterContainer[0].ncat(i);

    fprintf(file, "Number, Sample, Counter Name, Categories");
    for (unsigned int iCat=0; iCat<ncat; iCat++) 
      fprintf(file, "Cat %d Counts", iCat);
    fprintf(file, ", TOT Cat Counts");
    for (unsigned int iCat=0; iCat<ncat; iCat++)
      fprintf(file, ", Cat %d Tot Events", iCat);
    fprintf(file, ", TOT Cat Tot Events");
    for (unsigned int iCat=0; iCat<ncat; iCat++)
      fprintf(file, ", Cat %d Tot Sigma", iCat);
    fprintf(file, ", TOT Cat Tot Sigma");
    for (unsigned int iEff=0; iEff<3; iEff++) {
      fprintf(file, ", Denominator Name");
      for (unsigned int iCat=0; iCat<ncat; iCat++)
	      fprintf(file, ", Cat %d Eff.", iCat);
      fprintf(file, ", TOT Cat Eff.");
    }
    fprintf(file, ",, indexfiles, namefile, scale, lumi, intlum, weight");
    fprintf(file, "\n");
    
    for (int c=0; c<ncat; c++) {
      if(LDEBUG) std::cout<<"cat is "<< c <<std::endl;
      for (unsigned int j=0; j<counterContainer.size(); j++) {
        if (c == 0) {
	  counterNumerator_tot[j] = 0;
	  counterEvents_tot[j] = 0;
	  counterSigma_tot[j] = 0;
	      }
	      
	      for (unsigned int iEff=0; iEff<3; iEff++)
	        counterDenominator_tot[j][iEff] = 0;
	      
	      int indexfiles = sampleContainer[j].itype;
	      //indexinfo =; // check matteo
	      float weight = sampleContainer[j].weight;
	      
	      if (c == 0) {
	        fileLinesCatFile[j] << j << sampleContainer[j].filesshortnam << counterContainer[j].name(i) << ncat;
	        //fileLinesInfo[j] << indexfiles << files << scale << lumireal << intlumi << weight; // check matteo
	        fileLinesInfo[j] << indexfiles << sampleContainer[j].scale << sampleContainer[j].lumireal << utilInstance->intlumi << sampleContainer[j].weight;
	      }

	      float counts = counterContainer[j][i][c];
	      
	      // CHECK MATTEO
	      counterNumerator_tot[j] += counts;
	      counterEvents_tot[j] += counts*weight;
	      counterSigma_tot[j] += counts*weight/utilInstance->intlumi;
      }
      if(LDEBUG) std::cout<<"end cat loop"<<std::endl;
    }
    if(LDEBUG) std::cout<<"end sample loop"<<std::endl;
  }
  fclose(file);
  if(LDEBUG) std::cout<<"LoopAll::myWriteCounters - END"<<std::endl;
}

int LoopAll::FillAndReduce(int jentry) {

  int hasoutputfile = 0;
  if(utilInstance->typerun == 1) {
    hasoutputfile = 1;
    if(LDEBUG) 
      cout<<"call myReduce"<<endl;
    myGetEntryPhotonRedAnalysis(utilInstance, jentry);
    myReducePhotonAnalysis(utilInstance, jentry);
    if(LDEBUG) 
      cout<<"called myReduce"<<endl;
  } else if (utilInstance->typerun == 0) {
    hasoutputfile = 0;
    if(LDEBUG) 
      cout<<"call myFillHist -- really?"<<endl;
    myFillHistPhotonAnalysis(utilInstance, jentry);
    if(LDEBUG) 
      cout<<"called myFillHist"<<endl;
  } else if (utilInstance->typerun == 2) {
    hasoutputfile = 0;
    if(LDEBUG) 
      cout<<"call myFillHistRed"<<endl;
    myGetEntryPhotonRedAnalysis(utilInstance, jentry);
    myFillHistPhotonAnalysisRed(utilInstance, jentry);
    if(LDEBUG) 
      cout<<"called myFillHistRed"<<endl;
  } else if (utilInstance->typerun == 3) {
    hasoutputfile = 0;
    if(LDEBUG) 
      cout<<"call myFillHistRed"<<endl;
    myGetEntryPhotonRedAnalysis(utilInstance, jentry);
    myStatPhotonAnalysis(utilInstance, jentry);
    if(LDEBUG) 
      cout<<"called myFillHistStat"<<endl;
  }

  
  return hasoutputfile;
}
void LoopAll::BookHisto
(
  int h2d
 ,int typplot
 ,int typeplotall
 ,int histoncat
 ,int nbinsx
 ,int nbinsy
 ,float lowlim
 ,float highlim
 ,float lowlim2
 ,float highlim2
 ,char *name
){

    for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
      if (nbinsy == 0)
	histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim);
      if (nbinsy != 0)
	histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim, nbinsy, lowlim2, highlim2);
    }
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
void LoopAll::AddCut(
  char *cutnamesc
 ,int ncatstmp 
 ,int ifromright 
 ,int ifinalcut 
 ,float *cutValuel
 ,float *cutValueh
) {

  if(LDEBUG) cout<<"InitCuts START"<<endl;

    Cut* this_cut = new Cut();
    this_cut->name = cutnamesc;
    this_cut->fromright = ifromright;
    this_cut->finalcut = ifinalcut;
    this_cut->ncat = ncatstmp;
    this_cut->cut.clear();
    this_cut->cutintervall.clear();
    this_cut->cutintervalh.clear();
    if(LDEBUG) cout<<"this_cut filled  "<<endl;
    for (int j=0; j<ncatstmp; j++) {
      if(ifromright == 2) {
	this_cut->cutintervall.push_back((float) cutValuel[j]);
	this_cut->cutintervalh.push_back((float) cutValueh[j]);
      } else {
	this_cut->cut.push_back((float) cutValuel[j]);
      }
    }
    if(LDEBUG) cout<<"push back cut  "<<endl;
    cutContainer.push_back(*this_cut);
    if(LDEBUG) cout<<"pushed back cut  "<<endl;
  if(LDEBUG) cout<<"InitCuts END"<<endl;
}

void LoopAll::InitCounters(){

  if(LDEBUG) cout<<"InitCounts BEGIN"<<endl;

  for(unsigned int i=0; i<sampleContainer.size(); i++)
     counterContainer.push_back(CounterContainer(i));

  if(LDEBUG) cout<<"InitCounts END"<<endl;
}
void LoopAll::AddCounter
(
  int countersncat
 ,char *countername
 ,char *denomname0
 ,char *denomname1
 ,char *denomname2
){

    
    std::string* counternames_str = new std::string(countername);
    if(LDEBUG) cout<<" counternames_str"<<counternames_str<< endl; 
    //counternames_str.assign(counternames);
    if(LDEBUG) cout<<" *counternames_str"<<*counternames_str<< endl; 

    std::string* denomname0_str = new std::string(denomname0);
    if(LDEBUG) cout<<" denomname0_str"<<denomname0_str<< endl; 

    //std::string denomname1_str;
    std::string* denomname1_str =new std::string(denomname1);
    //denomname1_str.assign(denomname1);
    if(LDEBUG) cout<<" denomname1_str"<<denomname1_str<< endl; 

    std::string* denomname2_str = new std::string(denomname2);
    if(LDEBUG) cout<<" denomname2_str"<<denomname2_str<< endl; 

    for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
      if(LDEBUG) cout<<"adding to "<<ind<<" sampleContainer"<< endl; 
      counterContainer[ind].Add(*counternames_str
				,countersncat
				,*denomname0_str
				,*denomname1_str
				,*denomname2_str);
      if(LDEBUG) cout<<"added to "<<ind<<" sampleContainer"<< endl; 
    }
  //}
}

 
int LoopAll::ApplyCut(int icut, float var, int icat) {
  
  //returns 0 if not initialized
  if(cutContainer[icut].useit==0) return 1;

  if(cutContainer[icut].fromright==2) {
    if(var<cutContainer[icut].cutintervall[icat] || var>cutContainer[icut].cutintervalh[icat]) return 0;
  }
  else if (cutContainer[icut].fromright==1) {
    if(var>cutContainer[icut].cut[icat]) return 0;
  }
  else if (cutContainer[icut].fromright==0) {
    if(var<cutContainer[icut].cut[icat]) return 0;
  }
  return 1;
}

int LoopAll::ApplyCut(std::string cutname, float var, int icat) {
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name == cutname) {
      return ApplyCut(i, var, icat);
    }
  }
  //std::cout<<"ApplyCut: attention cutname "<<cutname<<" not found"<<endl;
  return 0;
}
 
void LoopAll::FillHist(std::string name, float y) {
  FillHist(name, 0, y);
}

void LoopAll::FillHist2D(std::string name, float x, float y) {
  FillHist2D(name, 0, x, y);
}

void LoopAll::FillHist(std::string name, int category, float y) {
  histoContainer[utilInstance->current].Fill(name, category, y);
}

void LoopAll::FillHist2D(std::string name, int category, float x, float y) {
  histoContainer[utilInstance->current].Fill2D(name, category, x, y);
}

void LoopAll::FillCounter(std::string name) {
  FillCounter(name, 0);
}

void LoopAll::FillCounter(std::string name, int category) {
  counterContainer[utilInstance->current].Fill(name, category);
}
