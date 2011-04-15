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

void Util::CallAddCut(
  char  *cutnamesc
 ,int   ncatstmp 
 ,int   ifromright 
 ,int   ifinalcut 
 ,float *cutValuel
 ,float *cutValueh
){ loops->AddCut(
	 cutnamesc
 	,ncatstmp 
 	,ifromright
 	,ifinalcut 			
 	,cutValuel
 	,cutValueh	
   ); 
}

void Util::CallInitHistos(){ loops->InitHistos(); }
void Util::CallBookHisto(
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
){  loops->BookHisto( 
	h2d
 	,typplot
 	,typeplotall
 	,histoncat
 	,nbinsx
 	,nbinsy
 	,lowlim
 	,highlim
 	,lowlim2
 	,highlim2
 	,name
    );

}
void Util::CallInitCounters(){ loops->InitCounters(); }
void Util::CallAddCounter(
  int countersncat
 ,char *countername
 ,char *denomname0
 ,char *denomname1
 ,char *denomname2
){  loops->AddCounter(
	countersncat	
	,countername	
	,denomname0
	,denomname1
	,denomname2
    );
}

void Util::SetTypeRun(int t, const char* name) {
  typerun = t;
  cout << "Type Run =" <<t <<endl; 
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

    outputParParameters = new std::vector<std::string>;
    outputParJobMaker = new std::string;
    outputTreePar->Branch("tot_events", &outputParTot_Events, "tot_events/I");
    outputTreePar->Branch("sel_events", &outputParSel_Events, "sel_events/I");
    outputTreePar->Branch("type", &outputParType, "type/I");
    outputTreePar->Branch("version", &outputParVersion, "version/I");
    outputTreePar->Branch("parameters", "std::vector<string>", &outputParParameters);
    outputTreePar->Branch("job_maker", &outputParJobMaker, "job_maker/C");
    //outputTreePar->Branch("job_maker", "std::string", &outputParJobMaker);
    outputTreePar->Branch("reductions", &outputParReductions, "reductions/I");
    outputTreePar->Branch("red_events", &outputParRed_Events, "red_events[reductions]/I");

  } 
}

void Util::ReadInput(int t) {
  // FIXME make this variable global ?
  typerun = t;

  
  if (typerun == 1) {
    makeOutputTree = 1;
  } else {
    makeOutputTree = 0;
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
  
} 
void Util::DefineSamples(
   const char *filesshortnam
  ,int type
  ,int histtoindfromfiles
  ,int histoplotit
  ,int nred
  ,long long ntot
  ,float intlumi
  ,float lumi
  ,float xsec
  ,float kfactor
  ,float scale
  ){


    //look in map for type as a key already
    int sample_is_defined = -1;
    for (unsigned int s=0; s<loops->sampleContainer.size(); s++) {
      if (type == loops->sampleContainer[s].itype) {
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

void Util::AddFile(std::string name,int type) {
  if(DEBUG) cout << "Adding file:  " << name << " of type " << type << endl;
  files.push_back(name);
  itype.push_back(type);
  
  nfiles++;
}

void Util::LoopAndFillHistos(TString treename) {


  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal "<<endl;

  loops->InitReal(typerun);
   
  Files.resize(files.size());  
  Trees.resize(files.size());  
  TreesPar.resize(files.size());  

  std::vector<std::string>::iterator it;
  std::vector<TTree*>::iterator it_tree;
  std::vector<TTree*>::iterator it_treepar;
  std::vector<TFile*>::iterator it_file;

  it 	  	= files.begin();
  it_file 	= Files.begin();
  it_tree	= Trees.begin();
  it_treepar    = TreesPar.begin();  

  for (;it!=files.end()
       ;it_file++,it_tree++,it_treepar++,it++){ 
 
    this->current = i;
    cout<<"LoopAndFillHistos: opening " << i << " " << files[i]<<endl;

    *it_file = TFile::Open((*it).c_str());
    //Files[i] = TFile::Open(files[i]);
    tot_events=1;
    sel_events=1;
    if(typerun == 1) { //this is a reduce job

      if(*it_file)
	*it_treepar=(TTree*) (*it_file)->Get("global_variables");

      if(*it_treepar) {
	TBranch        *b_tot_events;
	TBranch        *b_sel_events;
	(*it_treepar)->SetBranchAddress("tot_events",&tot_events, &b_tot_events);
	(*it_treepar)->SetBranchAddress("sel_events",&sel_events, &b_sel_events);
	b_tot_events->GetEntry(0);
	b_sel_events->GetEntry(0);
      } else {
	cout<<"REDUCE JOB, no global_variables tree ... the C-step job must have crashed ... SKIP THE FILE"<<endl;
	tot_events=0;
	sel_events=0;
      }
    }

    if(tot_events!=0) {

      if(*it_file)
	*it_tree=(TTree*) (*it_file)->Get(treename);

      loops->Init(typerun, *it_tree);
    }

    loops->Loop(i);

    if(tot_events != 0) {
      (*it_tree)->Delete("");
    }
   
    // EDIT - Cannot close the first file since it is in use after 
    // file 0 
    if(*it_file && i>0)
      (*it_file)->Close();
    
    i++;
  }
  //now close the first File
  if(Files[0]) Files[0]->Close();

  loops->TermReal(typerun);
}

void Util::WriteHist() {
  loops->myWritePlot();
}

void Util::WriteFits() {
  loops->myWriteFits();
}

void Util::WriteCounters() {
  loops->myWriteCounters();
}
