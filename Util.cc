#include "Util.h"
#include "LoopAll.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#define DEBUG 0

Util::Util() {
  loops = new LoopAll();
  loops->setUtilInstance(this);
#include "branchdef/newclonesarray.h"
}

void Util::SetTypeRun(int t, const char* name) {
  typerun = t;
  
  if (t == 1) {
    makeOutputTree = 1;
    outputFileName = TString(name);
    histFileName  = TString("");
  } else {
    makeOutputTree = 0;
    outputFileName = TString("");
    histFileName  = TString(name);
  }
  
  if (makeOutputTree) {
    cout << "CREATE " << outputFileName<<endl;
    outputFile=new TFile(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)
	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == 0 || typerun == 2) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == 1) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    
    hasoutputfile=1;
    outputTree = new TTree("event","Reduced tree"); 
    outputTreePar = new TTree("global_variables","Parameters");
    
    loops->PhotonAnalysisReducedOutputTree();
  } 
}

void Util::SetOutputNames(const char* name, const char* name2) {

  if (name2 != "")
    makeOutputTree = 1;
  histFileName  = TString(name);
  outputFileName = TString(name2);
}

void Util::ReadInput(int t) {
  // FIXME make this variable global ?
  float intlumi;
  typerun = t;

  char* outFilNam = new char[1024]; 
  char* histFilName = new char[1024]; 
  FILE *file;
  int dummy; 
  
  if (typerun == 1) {
    makeOutputTree = 1;
    file = fopen("filestoreduce.dat", "r");
    dummy = fscanf(file, "%d %s\n", &nfiles, outFilNam);
    outputFileName = TString(outFilNam);
    histFileName  = TString("");
  } else {
    file = fopen("inputfiles.dat", "r");
    dummy = fscanf(file, "%d intL=%f rtree=%d %s %s\n", &nfiles, &intlumi, &makeOutputTree, outFilNam, histFilName);
    makeOutputTree = 0;
    outputFileName = TString("");
    histFileName  = TString(histFilName);
  }
  
  if (makeOutputTree) {
    cout << "CREATE " << outputFileName<<endl;
    outputFile=new TFile(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == 0 || typerun == 2) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == 1) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    
    hasoutputfile=1;
    outputTree = new TTree("event","Reduced tree"); 
    outputTreePar = new TTree("global_variables","Parameters");
    
    loops->PhotonAnalysisReducedOutputTree();
  } 
  
  
  for (int i=0; i<nfiles; i++) {
    //lumireal[i]=0.;
    files[i] = new char[1024];
    char* filesshortnam = new char[1024];
    int histoplotit, nred, histoindfromfiles;
    long long int ntot;
    float lumi, xsec, kfactor, scale;
    if(typerun==1) {
      dummy = fscanf(file,"typ=%d Fil=%s\n", &itype[i], files[i]);
    } else {
      dummy = fscanf(file,"typ=%d ind=%d draw=%d Nam=%s Fil=%s tot=%lld red=%d lum=%f xsec=%f kfac=%f scal=%f\n", &itype[i], &histoindfromfiles, &histoplotit, filesshortnam, files[i], &ntot, &nred, &lumi, &xsec, &kfactor, &scale);
    }

    cout<<"after reading: files[i] "<<i<<" "<<files[i]<<endl;

    //look in map for type as a key already
    int sample_is_defined = -1;
    for (unsigned int s=0; s<loops->sampleContainer.size(); s++) {
      if (itype[i] == loops->sampleContainer[s].itype) {
	sample_is_defined = s;
	break;
      }
    }

    if (sample_is_defined != -1) {
      loops->sampleContainer[sample_is_defined].ntot += ntot;
      loops->sampleContainer[sample_is_defined].nred += nred;
      loops->sampleContainer[sample_is_defined].computeWeight(intlumi);
    } else {
      loops->sampleContainer.push_back(SampleContainer());
      loops->sampleContainer.front().ntot = ntot;
      loops->sampleContainer.front().nred = nred;
      loops->sampleContainer.front().histoplotit = histoplotit;
      loops->sampleContainer.front().filesshortnam = filesshortnam;
      loops->sampleContainer.front().lumi = lumi;
      loops->sampleContainer.front().xsec = xsec;
      loops->sampleContainer.front().kfactor = kfactor;
      loops->sampleContainer.front().scale = scale;
      loops->sampleContainer.front().computeWeight(intlumi);
    }
  }
  
  fclose(file);
}

void Util::AddFile(char* name,int type) {
  if(DEBUG) cout << "Adding file:  " << name << " of type " << type << endl;
  files[nfiles] = name;
  itype[nfiles] = type;
  
  /*
  //look in map for type as a key already
  map<int,int>::iterator it;
  if(DEBUG) cout << "map ok" << endl;
  it = type2HistVal.find(type);
  if(DEBUG) cout << "map checked" << endl;
  if (it == type2HistVal.end()) { 
    type2HistVal[type] = ntypes;
    if(DEBUG) cout << "container to add" << endl;
    ntypes++;
    if(DEBUG) cout << "container added" << endl;
  }
  */
  nfiles++;
}

void Util::LoopAndFillHistos(TString treename) {

  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal "<<endl;

  loops->InitReal(typerun);
  
  while (i<nfiles) {
    
    // MATTEO CHECK
    for (unsigned int s=0; s<loops->sampleContainer.size(); s++) {
      if (loops->sampleContainer[s].itype == itype[i])
	current = s;
      break;
    }

    cout<<"LoopAndFillHistos: opening " << i << " " << files[i]<<endl;

    Files[i] = TFile::Open(files[i]);
    tot_events=1;
    sel_events=1;
    if(typerun == 1) { //this is a reduce job

      if(Files[i])
	TreesPar[i]=(TTree*) Files[i]->Get("global_variables");

      if(TreesPar[i]) {
	TBranch        *b_tot_events;
	TBranch        *b_sel_events;
	TreesPar[i]->SetBranchAddress("tot_events",&tot_events, &b_tot_events);
	TreesPar[i]->SetBranchAddress("sel_events",&sel_events, &b_sel_events);
	b_tot_events->GetEntry(0);
	b_sel_events->GetEntry(0);
      } else {
	cout<<"REDUCE JOB, no global_variables tree ... the C-step job must have crashed ... SKIP THE FILE"<<endl;
	tot_events=0;
	sel_events=0;
      }
    }

    if(tot_events!=0) {

      if(Files[i])
	Trees[i]=(TTree*) Files[i]->Get(treename);

      loops->Init(typerun, Trees[i]);
    }

    loops->Loop();

    if(tot_events != 0) {
      Trees[i]->Delete("");
    }
    
    if(Files[i])
      Files[i]->Close();
    
    i++;
  }
 
  loops->TermReal(typerun);
}

void Util::WriteHist() {
  loops->myWritePlot();
}

void Util::WriteCounters() {
  loops->myWriteCounters();
}
