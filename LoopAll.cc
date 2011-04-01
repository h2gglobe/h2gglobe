#define DEBUG  0

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

LoopAll::LoopAll(TTree *tree) {}

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


void LoopAll::Loop(int a) {
  
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

void LoopAll::InitReal(Int_t typerunpass) {

  // Set branch addresses
  typerun=typerunpass;
  hfile = (TFile*)gROOT->FindObject(utilInstance->histFileName); 
  if (hfile) 
    hfile->Close();

  for(int ind=0; ind<sampleContainer.size(); ind++) {
    HistoContainer temp(ind);
    histoContainer.push_back(temp);
  }

  if(DEBUG) cout << "doing InitRealPhotonAnalysis" << endl;
  InitRealPhotonAnalysis(typerun);
  if(DEBUG) cout << "finished InitRealPhotonAnalysis" << endl;

  if (utilInstance->makeOutputTree) 
    utilInstance->outputFile->cd();
}

void LoopAll::TermReal(Int_t typerunpass) {

  typerun=typerunpass;

  TermRealPhotonAnalysis(typerun);
  
  if (utilInstance->makeOutputTree) 
    utilInstance->outputFile->cd();
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

void LoopAll::Loop() {
  
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
    
    if(DEBUG) 
      cout<<"call LoadTree"<<endl;
    
    Int_t ientry = LoadTree(jentry);
  
    if (ientry < 0) 
      break;
   
    if(DEBUG) 
      cout<<"Loopall_cc RRRRRRRRRRRR calling GetEntry"<<endl;
      
    if(utilInstance->typerun == 1) {
      nb = fChain->GetEntry(jentry);
      
      if(DEBUG) 
	cout<<"Loopall_cc RRRRRRRRRRR called GetEntry"<<endl;
    } else {
      //HERE NEED TO DO A GLOBAL GETENTRY
      nb=0;
    }

    nbytes += nb;

    if(DEBUG) 
      cout<<"Call FillandReduce "<<endl;
      
    hasoutputfile = FillAndReduce(jentry);
    if(DEBUG) 
      cout<<"Called FillandReduce "<<endl;
  }

  if(hasoutputfile) {
    if(outputFile) {
      outputFile->cd();
      if (DEBUG)
	cout<<"LoopAll_cc writing outputTree"<<endl;
      outputTree->Write(0,TObject::kWriteDelete);
    }

    
    if(outputFile) {
      outputFile->cd();
      if(utilInstance->TreesPar[0]) {
	std::vector<std::string> *parameters = new std::vector<std::string>;
	Int_t tot_events, sel_events, type, version, reductions;
	Int_t red_events[20];

	utilInstance->TreesPar[0]->SetBranchAddress("tot_events", &tot_events);
	utilInstance->TreesPar[0]->SetBranchAddress("sel_events", &sel_events);
	utilInstance->TreesPar[0]->SetBranchAddress("type", &type);
	utilInstance->TreesPar[0]->SetBranchAddress("version", &version);
	utilInstance->TreesPar[0]->SetBranchAddress("parameters", &parameters);
	if (utilInstance->TreesPar[0]->FindBranch("reductions")) {
	  utilInstance->TreesPar[0]->SetBranchAddress("reductions", &reductions);
	  utilInstance->TreesPar[0]->SetBranchAddress("red_events", &red_events);
	}
	
	TTree* newtree = new TTree("global_variables", "Global Parameters");
	newtree->Branch("tot_events", &tot_events, "tot_events/I");
	newtree->Branch("sel_events", &sel_events, "sel_events/I");
	newtree->Branch("type", &type, "type/I");
	newtree->Branch("version", &version, "version/I");
	newtree->Branch("parameters", "std::vector<string>", &parameters);
	newtree->Branch("reductions", &reductions, "reductions/I");
	newtree->Branch("red_events", &red_events, "red_events[reductions]/I");
	
	utilInstance->TreesPar[0]->GetEntry(0);

	if (!utilInstance->TreesPar[0]->FindBranch("reductions")) {
	  reductions = 1;
	  for (int i=0; i<20; i++) {
	    red_events[i] = -1;
	  }
	  red_events[0] = (int)countersred[1];
	} else {
	  red_events[reductions] = (int)countersred[1];
	  reductions++;
	}

	newtree->Fill();
	newtree->Write();
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
  if(DEBUG) std::cout<<"LoopAll::myWriteCounters - START"<<std::endl;
  
  const int samples = sampleContainer.size();
   
  if(DEBUG) std::cout<<"samples = "<<samples<<std::endl;
  
  stringstream fileLinesInfo[samples];
  stringstream fileLinesCatFile[samples];
  stringstream fileLinesCatCounts[samples];
  stringstream fileLinesCatEvents[samples];
  stringstream fileLinesCatSigma[samples];
  stringstream fileLinesCatEff[samples][3];
  
  if(DEBUG) std::cout<<"initialize streams"<<std::endl;

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

  if(DEBUG) std::cout<<"msName is "<< msName <<std::endl;
  if(DEBUG) std::cout<<"s2 is "<< s2 <<std::endl;
  
  for(unsigned int i=0; i<counterContainer[0].mapSize(); i++) {
    if(DEBUG) std::cout<<"counterContainer[0].mapSize() is "<< counterContainer[0].mapSize() <<std::endl;
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
      if(DEBUG) std::cout<<"cat is "<< c <<std::endl;
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
      if(DEBUG) std::cout<<"end cat loop"<<std::endl;
    }
    if(DEBUG) std::cout<<"end sample loop"<<std::endl;
  }
  fclose(file);
  if(DEBUG) std::cout<<"LoopAll::myWriteCounters - END"<<std::endl;
}

int LoopAll::FillAndReduce(int jentry) {

  int hasoutputfile = 0;
  if(utilInstance->typerun == 1) {
    hasoutputfile = 1;
    if(DEBUG) 
      cout<<"call myReduce"<<endl;
    myReducePhotonAnalysis(utilInstance, jentry);
    if(DEBUG) 
      cout<<"called myReduce"<<endl;
  } else if (utilInstance->typerun == 0) {
    hasoutputfile = 0;
    if(DEBUG) 
      cout<<"call myFillHist -- really?"<<endl;
    myFillHistPhotonAnalysis(utilInstance, jentry);
    if(DEBUG) 
      cout<<"called myFillHist"<<endl;
  } else if (utilInstance->typerun == 2) {
    hasoutputfile = 0;
    if(DEBUG) 
      cout<<"call myFillHistRed"<<endl;
    myFillHistPhotonAnalysisRed(utilInstance, jentry);
    if(DEBUG) 
      cout<<"called myFillHistRed"<<endl;
  }
  
  return hasoutputfile;
}

void LoopAll::BookHistos() {
  int typplotall, idummystart;
  int iread, h2d, typplot, histoncat, nbinsx, nbinsy, Nvar;
  float lowlim, highlim, lowlim2, highlim2;
  char name[100];
  FILE *file;
  file = fopen("plotvariables.dat", "r");
  int dummy = fscanf(file,"%d plot=%d idummy=%d", &Nvar, &typplotall, &idummystart);
  
  for (int i=0; i<Nvar; i++) {
    dummy = fscanf(file,"%d htyp=%d plot=%d ncat=%d %d %d %f %f %f %f %s", &iread, &h2d, &typplot, &histoncat, &nbinsx, &nbinsy, &lowlim, &highlim, &lowlim2, &highlim2, name);

    for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
      if (nbinsy == 0)
	histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim);
      if (nbinsy != 0)
	histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim, nbinsy, lowlim2, highlim2);
    }
  }
}

void LoopAll::InitCuts() {
  if(DEBUG) cout<<"InitCuts START"<<endl;

  int Ncuts;
  int iread=0, ncatstmp, ifromright, ifinalcut;
  float cutValuel[100], cutValueh[100], cutValue[100];
  FILE *file;
  char cutnamesc[100];
  file = fopen("cuts.dat", "r");
  int dummy = fscanf(file,"%d\n",&Ncuts);
  if(DEBUG) cout<<"reading "<<Ncuts<<" cuts"<<endl;
  for (int i=0; i<Ncuts; i++) {
    if(DEBUG) cout<<"reading cut "<<i<<endl;
    dummy = fscanf(file,"%d  %s  ncat=%d  dir=%d  fin=%d", &iread, cutnamesc, &ncatstmp, &ifromright, &ifinalcut);
    if(DEBUG) cout<<"reading cut:  " << cutnamesc <<endl;
    if(DEBUG) cout<<"ncatstmp: " << ncatstmp <<endl;
    if(DEBUG) cout<<"fromright:  " << ifromright <<endl;
    if(DEBUG) cout<<"ifinalcut:  " << ifinalcut <<endl;
    for (int j=0; j<ncatstmp; j++) {
      if(ifromright == 2) {
	int dummy = fscanf(file,"%f %f", &cutValuel[j], &cutValueh[j]);
      } else {
	int dummy = fscanf(file,"%f", &cutValue[j]);
  if(DEBUG) cout<<"cutValue[j]:  " << cutValue[j] <<endl;
      }
    }

    Cut* this_cut = new Cut();
    this_cut->name = cutnamesc;
    this_cut->fromright = ifromright;
    this_cut->finalcut = ifinalcut;
    this_cut->ncat = ncatstmp;
    this_cut->cut.clear();
    this_cut->cutintervall.clear();
    this_cut->cutintervalh.clear();
    if(DEBUG) cout<<"this_cut filled  "<<endl;
    for (int j=0; j<ncatstmp; j++) {
      if(ifromright == 2) {
	this_cut->cutintervall.push_back((float) cutValuel[j]);
	this_cut->cutintervalh.push_back((float) cutValueh[j]);
      } else {
	this_cut->cut.push_back((float) cutValue[j]);
      }
    }
    if(DEBUG) cout<<"push back cut  "<<endl;
    cutContainer.push_back(*this_cut);
    if(DEBUG) cout<<"pushed back cut  "<<endl;
  }
  if(DEBUG) cout<<"InitCuts END"<<endl;
}

void LoopAll::InitCounters() {  
  if(DEBUG) cout<<"InitCounts BEGIN"<<endl;

  for(unsigned int i=0; i<sampleContainer.size(); i++)
     counterContainer.push_back(CounterContainer(i));

  int Ncounters, ncounterscat;
  int iread, countersncat;
  float lowlim, highlim, lowlim2, highlim2;
  FILE *file;
  file = fopen("counters.dat", "r"); 
  int dummy = fscanf(file, "%d", &Ncounters);
  ncounterscat=0;
  char counternames[100];
  char denomname0[100];
  char denomname1[100];
  char denomname2[100];
   
  for(int i=0; i<Ncounters; i++) {
    //char* counternames = new char[1024];
    //char* denomname0 = new char[1024];
    //char* denomname1 = new char[1024];
    //char* denomname2 = new char[1024];
    if(DEBUG) cout<<"reading "<<i<<" counter"<< endl; 
    int dummy = fscanf(file,"%d ncat=%d %s %s %s %s",&iread, &countersncat, counternames, denomname0, denomname1, denomname2);
    if(DEBUG) cout<<"read "<<i<<" counter and dummy is "<<dummy<< endl; 
    
    std::string* counternames_str = new std::string(counternames);
    if(DEBUG) cout<<i<<" counternames_str"<<counternames_str<< endl; 
    //counternames_str.assign(counternames);
    if(DEBUG) cout<<i<<" *counternames_str"<<*counternames_str<< endl; 

    std::string* denomname0_str = new std::string(denomname0);
    if(DEBUG) cout<<i<<" denomname0_str"<<denomname0_str<< endl; 

    //std::string denomname1_str;
    std::string* denomname1_str =new std::string(denomname1);
    //denomname1_str.assign(denomname1);
    if(DEBUG) cout<<i<<" denomname1_str"<<denomname1_str<< endl; 

    std::string* denomname2_str = new std::string(denomname2);
    if(DEBUG) cout<<i<<" denomname2_str"<<denomname2_str<< endl; 

    for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
      if(DEBUG) cout<<"adding to "<<i<<" sampleContainer"<< endl; 
      counterContainer[ind].Add(*counternames_str, countersncat, *denomname0_str, *denomname1_str, *denomname2_str);
      //void Add(std::string, int, std::string, std::string, std::string);
      if(DEBUG) cout<<"added to "<<i<<" sampleContainer"<< endl; 
    }
  }
  if(DEBUG) cout<<"InitCounts END"<<endl;
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
