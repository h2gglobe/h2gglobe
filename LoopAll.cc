#define LDEBUG 0

#include "LoopAll.h"
#include "ErrorCodes.h"

#include <iostream>
#include <iterator>
#include <math.h>
#include <ctime>
#include "stdlib.h"

using namespace std;

#include "BaseAnalysis.h"

// ------------------------------------------------------------------------------------
BaseAnalysis* LoopAll::AddAnalysis(BaseAnalysis* baseAnalysis) {
  
  analyses.push_back(baseAnalysis);
  
  return baseAnalysis;
}

// ------------------------------------------------------------------------------------
void LoopAll::SetTypeRun(int t, const char* name) {
  typerun = t;
  if (typerun == kReduce || typerun == kFillReduce ) {
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
    outputFile=TFile::Open(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)
	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == kFill ) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == kReduce  || typerun == kFillReduce) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    
    hasoutputfile=1;
    outputFile->cd();
    outputTree = new TTree("event","Reduced tree"); 
    
    // book output braches
    Branches(outputBranchNames);
    for (size_t i=0; i<analyses.size(); i++) {
      analyses[i]->ReducedOutputTree(*this,outputTree);
    }

    outputTreePar = new TTree("global_variables","Parameters");
    outputParParameters = new std::vector<std::string>;
    outputParJobMaker = new std::string;
    outputTreePar->Branch("tot_events", &outputParTot_Events, "tot_events/I");
    outputTreePar->Branch("sel_events", &outputParSel_Events, "sel_events/I");
    outputTreePar->Branch("type", &outputParType, "type/I");
    outputTreePar->Branch("version", &outputParVersion, "version/I");
    outputTreePar->Branch("parameters", "std::vector<string>", &outputParParameters);
    outputTreePar->Branch("jobmaker", &outputParJobMaker, "jobmaker/C");
    //outputTreePar->Branch("job_maker", "std::string", &outputParJobMaker);
    outputTreePar->Branch("reductions", &outputParReductions, "reductions/I");
    outputTreePar->Branch("red_events", &outputParRed_Events, "red_events[reductions]/I");

  } 
  outputTreeLumi = new TTree("lumi","Lumi info tree"); 
  outputTreeLumi->Branch("run", &run, "run/I");
  outputTreeLumi->Branch("lumis", &lumis, "lumis/I");
}

// ------------------------------------------------------------------------------------
void LoopAll::ReadInput(int t) {
  // FIXME make this variable global ?
  typerun = t;

  
  if (typerun == kReduce || typerun == kFillReduce ) {
    makeOutputTree = 1;
  } else {
    makeOutputTree = 0;
  }
  
  if (makeOutputTree) {
    cout << "CREATE " << outputFileName<<endl;
    outputFile=TFile::Open(outputFileName,"recreate");
    outputFile->cd();
    if(outputFile) 
      if (DEBUG)	cout<<"output file defined"<<endl;
  }
  
  int hasoutputfile=0;
  
  if(typerun == kFill ) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == kReduce || typerun == kFillReduce ) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    
    hasoutputfile=1;
    outputTree = new TTree("event","Reduced tree"); 
    outputTreePar = new TTree("global_variables","Parameters");
    
    for (size_t i=0; i<analyses.size(); i++) {
      analyses[i]->ReducedOutputTree(*this,outputTree);
    }
    
    
  } 
  outputTreeLumi = new TTree("lumi","Lumi info tree"); 
  outputTreeLumi->Branch("run", &run, "run/I");
  outputTreeLumi->Branch("lumis", &lumis, "lumis/I");

  
} 

// ------------------------------------------------------------------------------------
SampleContainer & LoopAll::DefineSamples(const char *filesshortnam,
					 int type,
					 int histtoindfromfiles,
					 int histoplotit,
					 int nred,
					 long long ntot,
					 float intlumi,
					 float lumi,
					 float xsec,
					 float kfactor,
					 float scale,
					 bool ignoreEvWeight,
					 int forceVersion,
					 bool addnevents,
					 TString pileup
					 ) {
  
  // set intlumi in LoopAll
  intlumi_=intlumi;
  //look in map for type as a key already
  int sample_is_defined = -1;
  for (unsigned int s=0; s<sampleContainer.size(); s++) {
    if (type == sampleContainer[s].itype) {
      sample_is_defined = s;
      break;
    }
  }
  
  if (sample_is_defined != -1) {
    if( addnevents ) {
      sampleContainer[sample_is_defined].ntot += ntot;
      sampleContainer[sample_is_defined].nred += nred;
      sampleContainer[sample_is_defined].computeWeight(intlumi);
    }
    return sampleContainer[sample_is_defined];
  }

  sampleContainer.push_back(SampleContainer((ignoreEvWeight?0:&weight)));
  sampleContainer.back().itype = type;
  sampleContainer.back().ntot = ntot;
  sampleContainer.back().nred = nred;
  sampleContainer.back().histoplotit = histoplotit;
  sampleContainer.back().filesshortnam = filesshortnam;
  sampleContainer.back().lumi = lumi;
  sampleContainer.back().xsec = xsec;
  sampleContainer.back().kfactor = kfactor;
  sampleContainer.back().scale = scale;
  sampleContainer.back().forceVersion = forceVersion;
  sampleContainer.back().computeWeight(intlumi);
  if( pileup != "" ) { 
	  sampleContainer.back().pileup = pileup;
  }
  
  
  return sampleContainer.back();
}

// ------------------------------------------------------------------------------------
void LoopAll::AddFile(std::string name,int type) {
  if(DEBUG) 
    cout << "Adding file:  " << name << " of type " << type << endl;
  files.push_back(name);
  itype.push_back(type);
  
  nfiles++;
}

void LoopAll::MergeContainers(){
  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal " << endl;
  
  InitReal(typerun);
  
  int numberOfFiles = files.size();  
  Files.resize(numberOfFiles);  

  std::vector<std::string>::iterator it;
  std::vector<TFile*>::iterator it_file;
  
  it = files.begin();
  it_file = Files.begin();
  // Get names of objects inside RooContainer
  std::vector<std::string> histogramNames = rooContainer->GetTH1FNames();
  std::vector<std::string> datasetNames = rooContainer->GetDataSetNames();
		
  // Loop Over the files and get the relevant pieces to Merge:
  for (;it!=files.end()
	 ;it_file++,it++){
    
	  *it_file = TFile::Open((*it).c_str());
    (*it_file)->cd();
    std::cout << "Combining Current File " << i << " / " << numberOfFiles << " - " << (*it) << std::endl;

    for (std::vector<std::string>::iterator it_hist=histogramNames.begin()
	   ;it_hist!=histogramNames.end()
	   ;it_hist++) {
			
      TH1F *histExtra = (TH1F*) (*it_file)->Get(Form("th1f_%s",it_hist->c_str()));
      rooContainer->AppendTH1F(*it_hist,histExtra);	
      //delete histExtra;
    }
		
    RooWorkspace *work = (RooWorkspace*) (*it_file)->Get("cms_hgg_workspace");
    for (std::vector<std::string>::iterator it_data=datasetNames.begin()
	   ;it_data!=datasetNames.end()
	   ;it_data++) {

      RooDataSet *dataExtra = (RooDataSet*) work->data(Form("%s",it_data->c_str()));
      if( dataExtra == 0 ) {
	      std::cout << "skipping "<< it_data->c_str() << " " << dataExtra << std::endl;
	      continue;
      }
      rooContainer->AppendDataSet(*it_data,dataExtra);	
      //delete dataExtra;
    }

    delete work;
    //std::cout << "Finished Combining File - " << (*it) << std::endl;

    (*it_file)->Close();
    i++;
    //delete tmpFile;			
  } 
  TermReal(typerun);
  Term();
   
}
// ------------------------------------------------------------------------------------
void LoopAll::LoopAndFillHistos(TString treename) {

  int i=0;
  
  if (DEBUG)
    cout<<"LoopAndFillHistos: calling InitReal " << endl;
  
  InitReal(typerun);
 
  int numberOfFiles = files.size(); 
  Files.resize(numberOfFiles);  
  Trees.resize(numberOfFiles);  
  LumiTrees.resize(numberOfFiles);  
  TreesPar.resize(numberOfFiles);  

  std::vector<std::string>::iterator it;
  std::vector<TTree*>::iterator it_tree;
  std::vector<TTree*>::iterator it_treelumi;
  std::vector<TTree*>::iterator it_treepar;
  std::vector<TFile*>::iterator it_file;
  
  it = files.begin();
  it_file = Files.begin();
  it_tree	= Trees.begin();
  it_treelumi	= LumiTrees.begin();
  it_treepar = TreesPar.begin();  
  

  cout << "SAMPLE CONTAINER SIZE " << sampleContainer.size() <<endl;
  for (;it!=files.end()
	       ;it_file++,it_tree++,it_treelumi++,it_treepar++,it++,i++){ 
    
    this->current = i;
	
    int type = itype[i];
    int which_sample = -1;
    for (int sample=0;sample<sampleContainer.size();sample++){
      if (type == sampleContainer[sample].itype){
	which_sample = sample;
      }
    }
    this->current_sample_index = which_sample;

    cout<<"LoopAndFillHistos: opening file " << i+1 << " / " << numberOfFiles << " : " << files[i]<<endl;
    
    for(int itry=0; itry<3; ++itry) {
	    *it_file = TFile::Open((*it).c_str(),"TIMEOUT=60");
	    if( *it_file == 0 ){
		    std::cerr << "Error opening file " << (*it).c_str() << " attempt " << itry << std::endl;
	    } else {
		    break;
	    }
    }
    if( *it_file == 0 ){
	    std::cerr << "Error opening file " << (*it).c_str() << std::endl;
	    exit(H2GG_ERR_FILEOP);
    }
    
    //Files[i] = TFile::Open(files[i]);
    tot_events=1;
    sel_events=1;
    
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
      cout<<"NO global_variables tree ... the C-step job must have crashed ... SKIP THE FILE"<<endl;
      tot_events=0;
      sel_events=0;
      continue;
    }
    
    *it_treelumi = (TTree*) (*it_file)->Get("lumi");
    if( *it_treelumi != 0 && outputFile ) {
      StoreProcessedLumis( *it_treelumi );
    }
    
    if(typerun == kReduce || typerun == kFillReduce) { //this is a reduce job
      
      // Cannot mix PU histograms from different samples
      //assert(sampleContainer.size()==1);

      if (type!=0 && outputFile){
	  if (i==0) {// this is the first file and so need to clone the pileup histo
	  	pileup = (TH1D*) ((*it_file)->Get("pileup"))->Clone();
	  }
	  else {	
	  	pileup->Add((TH1D*) ((*it_file)->Get("pileup")));
	  }
      }      
      
      if (type == 0)
	if (typerun==kReduce || typerun == kFillReduce )
	  pileup  =  new TH1D("pileup", "pileup", 100, 0, 100); 
    }
    
    if(tot_events!=0) {
      
      if(*it_file)
        *it_tree=(TTree*) (*it_file)->Get(treename);
      
      Init(typerun, *it_tree);
    }

    if( *it_tree != 0 && (*it_tree)->GetEntries() != 0 ) {
	Loop(i);
    }
    
    if(tot_events != 0 && (*it_tree) != 0  ) {
      (*it_tree)->Delete("");
    }
    
    // EDIT - Cannot close the first file since it is in use after 
    // file 0 
    if (i>0)
      if((*it_file)->IsOpen())
        (*it_file)->Close();
  }
  
  TermReal(typerun);
  Term();
  
  //now close the first File
  if( !Files.empty() &&  Files[0]->IsOpen())
    Files[0]->Close();
}

// ------------------------------------------------------------------------------------
void LoopAll::StoreProcessedLumis(TTree * tree){
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumis",&lumis);
  for( int ii=0; ii<tree->GetEntries(); ++ii) {
    tree->GetEntry(ii);
    if( !CheckLumiSelection(run,lumis)  ) { continue; }
    if( outputTreeLumi ) outputTreeLumi->Fill();
  }
  if(typerun != kFill && outputFile) {
    outputFile->cd();
    outputTreeLumi->Write(0,TObject::kWriteDelete);
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::Term(){
  if (outputFile)
    if (outputFile->IsOpen())
      outputFile->Close();
}

// ------------------------------------------------------------------------------------
LoopAll::LoopAll(TTree *tree) :
	counters(4,0.), countersred(4,0.), checkBench(0), sqrtS(8)
{  
#include "branchdef/newclonesarray.h"

#ifndef __CINT__
#include "branchdef/branchdict.h"
  DefineUserBranches();
#endif

  rooContainer       = new RooContainer();
  signalNormalizer   = 0;
  /// signalNormalizer   = new Normalization_8TeV();

  rooContainer->BlindData();	// 2012 requires that we Blind our data
  // Best Set Global parameters accesible via python to defauls
  
  funcReader_dipho_MIT = 0;
  /// signalNormalizer->FillSignalTypes();

  runZeeValidation = false;
  makeDummyTrees = false;
  usePFCiC = true;
  applyEcalIsoPresel = false;
  pfisoOffset=2.5;
  cicVersion="7TeV";
  pho_r9_cic = &pho_r9[0];
  pho_idmva_cached = false;
}

// ------------------------------------------------------------------------------------
LoopAll::~LoopAll() {
  if (!fChain)
    return;
  //delete fChain->GetCurrentFile();
}

// ------------------------------------------------------------------------------------
Int_t LoopAll::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}


// ------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------------------
Bool_t LoopAll::Notify() {
  return kTRUE;
}

// ------------------------------------------------------------------------------------
void LoopAll::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

// ------------------------------------------------------------------------------------
void LoopAll::InitHistos(){

  for(int ind=0; ind<sampleContainer.size(); ind++) {
    SampleContainer thisSample = (SampleContainer) sampleContainer.at(ind);
    HistoContainer temp(ind,thisSample.filesshortnam);
    temp.setScale(1.);
    /// temp.setScale(thisSample.weight);
    histoContainer.push_back(temp);
  }

  HistoContainer temp(sampleContainer.size(),"tot");
  temp.setScale(1.);
  histoContainer.push_back(temp);
}

// ------------------------------------------------------------------------------------
void LoopAll::InitTrees(std::string dirName){

  std::vector<TreeContainer> vTemp;
  for(int ind=0; ind<sampleContainer.size(); ind++) {
    SampleContainer thisSample = (SampleContainer) sampleContainer.at(ind);
    TreeContainer temp(ind,thisSample.filesshortnam, dirName);
    vTemp.push_back(temp);
  }

  treeContainer[dirName] = vTemp;
}

// ------------------------------------------------------------------------------------
void LoopAll::InitReal(Int_t typerunpass) {

  // Set branch addresses
  typerun = typerunpass;
  hfile = (TFile*)gROOT->FindObject(histFileName); 
  if (hfile) 
    hfile->Close();

  //for(int ind=0; ind<sampleContainer.size(); ind++) {
  // HistoContainer temp(ind);
  //  histoContainer.push_back(temp);
  // }

  if(LDEBUG) cout << "doing InitRealPhotonAnalysis" << endl;
  for (size_t i=0; i<analyses.size(); i++) {
    analyses[i]->Init(*this);
  }
  
  if(LDEBUG) cout << "finished InitRealPhotonAnalysis" << endl;

  // Initialize all MVA -> This is done inside the Analysese now
  // SetAllMVA();

  if (makeOutputTree) 
    outputFile->cd();
  cout<< "LoopAll::InitReal END" <<endl;
  
}

// ------------------------------------------------------------------------------------
void LoopAll::TermReal(Int_t typerunpass) {

  typerun=typerunpass;

  for (size_t i=0; i<analyses.size(); i++) {
    analyses[i]->Term(*this);
  }

  if (makeOutputTree){ 
    outputFile->cd();
    assert( outputTree->GetEntries() == countersred[3] );
    outputTree->Write(0,TObject::kWriteDelete);
    outputParReductions++;
    outputTreePar->Fill();
    outputTreePar->Write(0,TObject::kWriteDelete);
    outputTreeLumi->Write(0,TObject::kWriteDelete);
    pileup->Write(0,TObject::kWriteDelete);
    for(int ii=0; ii<globalHistos.size(); ++ii) {
	    globalHistos[ii]->Write(0,TObject::kWriteDelete);
    }
  }
}


// ------------------------------------------------------------------------------------
void LoopAll::Init(Int_t typerunpass, TTree *tree) {
  
  // Set branch addresses
  typerun=typerunpass;
  
  fChain = tree;
  if (tree == 0) 
    return;

  fChain->SetMakeClass(1);
  
  // set-up inputs
  inputBranches.clear();
  GetBranches(inputBranchNames, inputBranches);
  if( typerun != kReduce && typerun != kFillReduce ) {
    for (size_t i=0; i<analyses.size(); i++) {
      analyses[i]->GetBranches(fChain, inputBranches);
    }
  }
  SetBranchAddresses(inputBranchNames);
  

  Notify();
}


// ------------------------------------------------------------------------------------
void bookGlobalCounters( TTree * intree, TTree * outTree, 
			 std::vector<std::string> & globalCountersNames, 
			 std::vector<int> & globalCounters, 
			 std::vector<int> & fileGlobalCounters )
{
  // ((TBranch *) global_variables->GetListOfBranches()->At(0)->GetAddress()
  TObjArray * inbranches = intree->GetListOfBranches();
  for( int ib=0; ib<inbranches->GetEntries(); ++ib ) {
    TBranch * br = (TBranch *)inbranches->At(ib);
    if( br->GetAddress() == 0 ) { 
      TString title = br->GetTitle();
      if( title.EndsWith("/I") 
	  && find( globalCountersNames.begin(), globalCountersNames.end(), br->GetName() ) == globalCountersNames.end() ) {
	globalCountersNames.push_back( std::string(br->GetName()) );
      }
    }
  }
  globalCounters.resize( globalCountersNames.size(), 0 );
  fileGlobalCounters.clear();
  fileGlobalCounters.resize( globalCountersNames.size(), 0 );
  for(size_t ic=0; ic<globalCountersNames.size(); ++ic ) {
    if( ! outTree->FindBranch( globalCountersNames[ic].c_str() ) ) {
      std::cerr << "Booking global counter " << globalCountersNames[ic] << std::endl;
      outTree->Branch( globalCountersNames[ic].c_str(), &globalCounters[ic], (globalCountersNames[ic]+"/I").c_str() );
    }
    std::cerr << "Reading global counter " << globalCountersNames[ic] << std::endl;
    intree->SetBranchAddress( globalCountersNames[ic].c_str(), &fileGlobalCounters[ic] );
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::Loop(Int_t a) {
  
  //makeOutputTree = makeOutputTree;
  if (makeOutputTree) {
    // outputTree = outputTree;
    // outputFile = outputFile;
    if(outputFile) 
      outputFile->cd();
  }

  Int_t nentries = 0;
  if(fChain) 
    nentries = Int_t(fChain->GetEntriesFast());


  Int_t nbytes = 0, nb = 0;

  outputEvents=0;

  int hasoutputfile=0;
  // Call the Reset Analysis at start of new file
  for (size_t i=0; i<analyses.size(); i++) {
    analyses[i]->ResetAnalysis(); 
  }

  time_t tfilestart,tfileend;
  tfilestart = time(0);
  
  /// read global parameters
  std::vector<std::string> *parameters = new std::vector<std::string>;
  std::string *job_maker = new std::string;
  Int_t red_events[20];
  
  TreesPar[a]->SetBranchAddress("tot_events", &tot_events);
  TreesPar[a]->SetBranchAddress("sel_events", &sel_events);
  TreesPar[a]->SetBranchAddress("type", &type);
  TreesPar[a]->SetBranchAddress("version", &version);
  TreesPar[a]->SetBranchAddress("parameters", &parameters);
  TreesPar[a]->SetBranchAddress("jobmaker", &job_maker);
  if (TreesPar[a]->FindBranch("reductions")) {
    TreesPar[a]->SetBranchAddress("reductions", &reductions);
    TreesPar[a]->SetBranchAddress("red_events", &red_events);
  }
  if(typerun == kReduce || typerun == kFillReduce  ){
    bookGlobalCounters( TreesPar[a], outputTreePar, globalCountersNames, globalCounters, fileGlobalCounters );
    assert( globalCountersNames.size() == globalCounters.size() && globalCounters.size() == fileGlobalCounters.size() );
  }
  TreesPar[a]->GetEntry(0);

  if( sampleContainer[current_sample_index].forceVersion > 0 ) { 
	  version = sampleContainer[current_sample_index].forceVersion;
  }
  
  // Loop over events
  if(checkBench > 0) {
	  stopWatch.Start();
  }
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    
    if(jentry%10000==0) {
      cout << "Entry: "<<jentry << " / "<<nentries <<  " "  ;
      copy(countersred.begin(), countersred.end(), std::ostream_iterator<float>(cout, "_") );
      cout << endl;
    }
    if(makeDummyTrees) continue;
    
    //if(checkBench > 0 && jentry%checkBench == 0 ) {
    //  stopWatch.Stop();
    //  float cputime  = stopWatch.CpuTime();
    //  float realtime = stopWatch.RealTime();
    //  stopWatch.Start(false);
    //  if( realtime > benchStart*100. && cputime / realtime < benchThr ) {
    //	std::cout << 
    //	  "\n\n\n\nAbourting exection.\n"
    //	  "Sorry: too inefficient to continue (cputime " << cputime << " realtime " << realtime << ")" 
    //		  << std::endl;
    //	exit(H2GG_ERR_DUTYC);
    //  }
    //}

    if(LDEBUG) 
      cout<<"call LoadTree"<<endl;
    
    Int_t ientry = LoadTree(jentry);
  
    if (ientry < 0) 
      break;
    
    if(typerun == kFill ) {
      nb=0;
    }

    nbytes += nb;

    if(LDEBUG) 
      cout<<"Call FillandReduce "<<endl;
      
    hasoutputfile = this->FillAndReduce(jentry);
    if(LDEBUG) 
      cout<<"Called FillandReduce "<<endl;
  }
  tfileend = time(0);
  std::cout << "Average time per event: " << (float) difftime(tfileend,tfilestart)/nentries << std::endl;

  //  if(hasoutputfile) {
  if(typerun == kReduce || typerun == kFillReduce  ){
    if(outputFile) {
      outputFile->cd();
      if (LDEBUG)
	cout<<"LoopAll_cc writing outputTree"<<endl;
      outputTree->Write(0,TObject::kWriteDelete);
    }

    //   if(outputFile) {
    outputFile->cd();
    //    if(TreesPar[a]) {
        
    if (a == 0) {
      outputParTot_Events = tot_events;
      outputParSel_Events = sel_events;

    } else {
      outputParTot_Events += tot_events;
      outputParSel_Events += sel_events;
    }
    for(size_t ii=0; ii<globalCountersNames.size(); ++ii) {
      globalCounters[ii] += fileGlobalCounters[ii];
      cout << "globalCounters Tot - " << globalCounters[ii]<<endl;
    }
	
    outputParType = type;
    outputParVersion = version;
    outputParParameters = &(*parameters);
    outputParJobMaker = job_maker;
        
    if (!TreesPar[a]->FindBranch("reductions")) {
      outputParReductions = 0;
      for (int i=0; i<20; i++) {
	outputParRed_Events[i] = 0;
      }
      outputParRed_Events[0] += (int)countersred[3];
    } else {
      outputParReductions = reductions;
      outputParRed_Events[reductions] += (int)countersred[3];
    }
  }
  
  int oldnentries=nentries;
  if (nentries == sel_events) {
    nentries = tot_events;
  }
  
  if(countersred[1] || oldnentries==0) {
    //printf("red: %d_%d \n",(int)countersred[0], (int) countersred[1]);
    cout << "red: ";
    copy(countersred.begin(), countersred.end(), std::ostream_iterator<float>(cout, "_") );
    cout << endl;
  } else { 
    // printf("norm: %d \n",(int)counters[0]);
    cout << "norm: "; 
    copy(countersred.begin(), countersred.end(), std::ostream_iterator<float>(cout, "_") );
    cout << endl;
  }
  if(checkBench > 0) {
	  stopWatch.Stop();
  }

}

// ------------------------------------------------------------------------------------
void LoopAll::WriteFits() {
  
  hfile = TFile::Open(histFileName, "RECREATE", "Globe ROOT file with histograms");
  
  hfile->cd();
  //hfile->cd();
  for (std::vector<TMacro*>::iterator it = configFiles.begin(); it!=configFiles.end(); it++){
    (*it)->Write();
  }
  
  rooContainer->Save();
  hfile->Close();
}

void LoopAll::StoreConfigFile(std::string configfilename) {
	 
  TMacro *mac = new TMacro(configfilename.c_str(),configfilename.c_str());
  configFiles.push_back(mac);
	 
}
// ------------------------------------------------------------------------------------
void LoopAll::WriteHist() {

  hfile = TFile::Open(histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  for(unsigned int ind=0; ind<histoContainer.size(); ind++) {
    histoContainer[ind].Save();
  }

  for (std::map<std::string, std::vector<TreeContainer> >::iterator ii = treeContainer.begin(); ii!=treeContainer.end(); ++ii) {
    for(unsigned int ind=0; ind<(*ii).second.size(); ind++) {
      ((*ii).second)[ind].Save(hfile);
      
    }
  }

  outputTreeLumi->Write();

  WritePI();

  hfile->Close();
      
  if (makeOutputTree) 
    outputFile->cd();

  //WriteCounters();
}


void LoopAll::WritePI() {
  Int_t Nvar;
  Int_t h2d[10000], typplot[10000], histoncat[10000], histoindfromfiles[10000];
  Int_t histoncatindtonames[10000];
  Int_t nbinsx[10000], nbinsy[10000];
  Float_t lowlim[10000], highlim[10000], lowlim2[10000], highlim2[10000];
  Int_t itype[10000], histoind[10000], infoind[10000], histoplotit[10000];
  Int_t ntot[10000], nred[10000];
  Float_t lumi[10000], xsec[10000], kfactor[10000], scale[10000];
  Int_t typplotall = 0;
  Int_t plothistoplotitPI[10000];

  hfile->cd();
  plotvartree = new TTree("plotvariables","globe plotvariables provenance information");
  
  plotvartree->Branch("Nvar", &Nvar, "Nvar/I");
  plotvartree->Branch("typplotall", &typplotall, "typplotall/I");
  plotvartree->Branch("doplot", &plothistoplotitPI, "plothistoplotitPI[Nvar]/I");
  plotvartree->Branch("h2d", &h2d, "h2d[Nvar]/I");
  plotvartree->Branch("typplot", &typplot, "typplot[Nvar]/I");
  plotvartree->Branch("histoncat", &histoncat, "histoncat[Nvar]/I");
  plotvartree->Branch("histoncatindtonames", &histoncatindtonames, "histoncatindtonames[Nvar]/I");
  plotvartree->Branch("nbinsx", &nbinsx, "nbinsx[Nvar]/I");
  plotvartree->Branch("nbinsy", &nbinsy, "nbinsy[Nvar]/I");
  plotvartree->Branch("lowlim", &lowlim, "lowlim[Nvar]/F");
  plotvartree->Branch("highlim", &highlim, "highlim[Nvar]/F");
  plotvartree->Branch("lowlim2", &lowlim2, "lowlim2[Nvar]/F");
  plotvartree->Branch("highlim2", &highlim2, "highlim2[Nvar]/F");

  Nvar = histoContainer[0].size();
  for (int i=0; i<Nvar; i++) {
    h2d[i]       = histoContainer[0].getDimension(i) - 1;
    plothistoplotitPI[i]    = 1;
    typplot[i]   = 1;
    histoncat[i] = histoContainer[0].ncat(i);
    histoncatindtonames[i] = -1;
    nbinsx[i]    = histoContainer[0].nbins(i, true);
    nbinsy[i]    = histoContainer[0].nbins(i, false);
    lowlim[i]    = histoContainer[0].min(i, true);
    highlim[i]   = histoContainer[0].max(i, true);
    lowlim2[i]   = histoContainer[0].min(i, false);
    highlim2[i]  = histoContainer[0].max(i, false);
  }
  
  TClonesArray* tca_xaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("xaxislabels", "TClonesArray", &tca_xaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_xaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_xaxislabels->At(iplot))->SetString(histoContainer[0].axisName(iplot, true).c_str());
  }

  TClonesArray* tca_yaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("yaxislabels", "TClonesArray", &tca_yaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_yaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_yaxislabels->At(iplot))->SetString(histoContainer[0].axisName(iplot, false).c_str());
  }

  TClonesArray* tca_plotvarnames = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("plotvarnames", "TClonesArray", &tca_plotvarnames, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_plotvarnames)[iplot]) TObjString(); 
    ((TObjString *)tca_plotvarnames->At(iplot))->SetString(histoContainer[0].getName(iplot).c_str());
  }

  /*
    plotvartree->Branch("Nvarcats", &Nvarcats, "Nvarcats/I");
    plotvartree->Branch("catid", &catid, "catid[Nvarcats]/I");
    plotvartree->Branch("ncats", &ncats, "ncats[Nvarcats]/I");
    tca_plotvarcatnames = new TClonesArray("TObjString",Ncatvar);
    plotvartree->Branch("plotvarcatnames", "TClonesArray", &tca_plotvarcatnames, 32000, 0);
    int catvartemp=0;
    for(int i=0; i<Nvarcats; i++) {
    for(int j=0; j<ncats[i]; j++) {
    new ((*tca_plotvarcatnames)[catvartemp]) TObjString(); 
    ((TObjString *)tca_plotvarcatnames->At(catvartemp))->SetString(catnames[i][j]);
    catvartemp++;
    } 
    } 
    std::cout << "Ncatvar: " << Ncatvar << std::endl;
    std::cout << "catvartemp: " << catvartemp << std::endl;
  */

  plotvartree->Fill();
  plotvartree->Write(0,TObject::kWriteDelete);
  inputfiletree = new TTree("inputfiles","globe inputfiles provenance information");
  int nfiles = files.size();
  //int nindfiles = type2HistVal.size();
  int nindfiles = sampleContainer.size();
  inputfiletree->Branch("nfiles", &nfiles, "nfiles/I");
  inputfiletree->Branch("nindfiles", &nindfiles, "nindfiles/I");
  inputfiletree->Branch("intlumi", &intlumi_, "intlumi/F");
  inputfiletree->Branch("makeOutputTree", makeOutputTree, "makeOutputTree/I");
  TClonesArray* tca_histfilename = new TClonesArray("TObjString",1);
  inputfiletree->Branch("histfilename", "TClonesArray", &tca_histfilename, 32000, 0);
  new ((*tca_histfilename)[0]) TObjString(); 
  ((TObjString *)tca_histfilename->At(0))->SetString(histFileName);
  inputfiletree->Branch("itype", &itype, "itype[nfiles]/I");
  inputfiletree->Branch("histoind", &histoindfromfiles, "histoindfromfiles[nfiles]/I");
  inputfiletree->Branch("infoind", &infoind, "infoind[nindfiles]/I");
  inputfiletree->Branch("histoplotit", &histoplotit, "histoplotit[nfiles]/I");
  inputfiletree->Branch("ntot", &ntot, "ntot[nfiles]/I");
  inputfiletree->Branch("nred", &nred, "nred[nfiles]/I");
  inputfiletree->Branch("lumi", &lumi, "lumi[nfiles]/F");
  inputfiletree->Branch("xsec", &xsec, "xsec[nfiles]/F");
  inputfiletree->Branch("kfactor", &kfactor, "kfactor[nfiles]/F");
  inputfiletree->Branch("scale", &scale, "scale[nfiles]/F");

  TClonesArray* tca_inshortnames = new TClonesArray("TObjString", sampleContainer.size());
  inputfiletree->Branch("inshortnames", "TClonesArray", &tca_inshortnames, 32000, 0);
  TClonesArray* tca_infilenames = new TClonesArray("TObjString", sampleContainer.size());
  inputfiletree->Branch("infilenames", "TClonesArray", &tca_infilenames, 32000, 0);

  for (unsigned int i=0; i<sampleContainer.size(); i++) {
    itype[i] = sampleContainer[i].itype;
    histoind[i] = sampleContainer[i].ind;
    infoind[i] = sampleContainer[i].ind;
    histoplotit[i] = sampleContainer[i].histoplotit;
    ntot[i] = sampleContainer[i].ntot; 
    nred[i] = sampleContainer[i].nred;
    lumi[i] = sampleContainer[i].lumi;
    xsec[i] = sampleContainer[i].xsec;
    kfactor[i] = sampleContainer[i].kfactor;
    scale[i] = sampleContainer[i].scale;
    new ((*tca_inshortnames)[i]) TObjString(); 
    ((TObjString *)tca_inshortnames->At(i))->SetString(sampleContainer[i].filesshortnam.c_str());
    
    new ((*tca_infilenames)[i]) TObjString(); 
    ((TObjString *)tca_infilenames->At(i))->SetString("n.a.");
  }

  inputfiletree->Fill();
  inputfiletree->Write(0,TObject::kWriteDelete);
}

// ------------------------------------------------------------------------------------
void LoopAll::WriteCounters() {

  // number of files
  int nindfiles=currentindexfiles;
  
  TString msName = histFileName;
  msName.ReplaceAll(".root", ".csv");

  FILE *file;
  TString s1 = msName;
  file = fopen(s1, "w");

  // FIXME size
  for (Int_t var=0; var<counterContainer[0].size(); var++) {
    fprintf(file, ",Sample,Counter Name,Categories,");
    Int_t cats = counterContainer[0].ncat(var);
    for (Int_t cat=0; cat<cats; cat++)
      fprintf(file, "Cat %d Counts,", cat);    
    fprintf(file, "TOT Cat Counts, , ,");
    for (Int_t cat=0; cat<cats; cat++)
      fprintf(file, "Cat %d Tot Events,", cat);
    fprintf(file, "TOT Cat Tot Events, , ,");
    for (Int_t cat=0; cat<cats; cat++)
      fprintf(file, "Cat %d Sigma,", cat);
    fprintf(file, "TOT Cat Sigma, , ,");
    for (Int_t den=0; den<3; den++) {
      fprintf(file, "Denominator Name,");
      for (Int_t cat=0; cat<cats; cat++)
	fprintf(file, "Cat %d Counts,", cat);    
      fprintf(file, "TOT Cat Counts, , ,");
    }
    fprintf(file, "indexfiles,scale,lumi,intLumi,weight\n");
    
    for (Int_t c=0; c<counterContainer.size(); c++) {
      fprintf(file, ",%s,%s,%d,", 
	      sampleContainer[c].filesshortnam.c_str(), 
	      counterContainer[0].name(var).c_str(),
	      cats);

      float tot = counterContainer[c].tot(var);
      for (Int_t cat=0; cat<cats; cat++) {
	fprintf(file, "%f,", counterContainer[c][var][cat]);
      }
      fprintf(file, "%f, , ,", tot);

      float weight = sampleContainer[c].weight();
      for (Int_t cat=0; cat<cats; cat++) {
	fprintf(file, "%f,", counterContainer[c][var][cat]*weight);
      }
      fprintf(file, "%f, , ,", tot*weight);

      // FIXME Sigma ???
      for (Int_t cat=0; cat<cats; cat++) {
	fprintf(file, "%f,", counterContainer[c][var][cat]*weight);
      }
      fprintf(file, "%f, , ,", tot*weight);

      for (Int_t den=0; den<3; den++) {
	fprintf(file, "%s ,", counterContainer[c].denomName(var, den).c_str());
	for (Int_t cat=0; cat<cats; cat++)
	  fprintf(file, "%f,", counterContainer[c].efficiency(var, cat, den));    
	fprintf(file, "%f, , ,", counterContainer[c].efficiency(var, den));
      }
      fprintf(file, "%d,%f,%f,%f,%f\n",
	      sampleContainer[c].itype,
	      sampleContainer[c].scale,
	      sampleContainer[c].lumi,
	      intlumi_,
	      weight);
    }

    fprintf(file, "\n\n\n\n");
  }

  fclose(file);
}





// ------------------------------------------------------------------------------------
int LoopAll::FillAndReduce(int jentry) {

  int hasoutputfile = 0;

  //count all events
  countersred[0]++;

  //
  // read all inputs 
  //
  if(!makeDummyTrees){
    GetEntry(inputBranches, jentry);
  }

  //b_run->GetEntry(jentry);
  //b_lumis->GetEntry(jentry);
  if(!CheckLumiSelection(run,lumis)){
    return hasoutputfile;
  }
  countersred[1]++;

  // 
  // call skimming methods before reading data
  // 
  for (size_t i=0; i<analyses.size(); i++) {
    if( ! analyses[i]->SkimEvents(*this, jentry) ) {
      return hasoutputfile;
    }
  }
  countersred[2]++;

  //
  // reduction step
  //
  if( typerun == kReduce || typerun == kFillReduce ) {
    hasoutputfile = 1;
    
    // compute additional quantites
    for (size_t i=0; i<analyses.size(); i++) {
      analyses[i]->FillReductionVariables(*this, jentry);
    }
    
    // (pre-)select events
    for (size_t i=0; i<analyses.size(); i++) {
      if( ! analyses[i]->SelectEventsReduction(*this, jentry) ) {
	return hasoutputfile;
      }
    }

    // fill output tree
    outputEvents++;
    if(LDEBUG) 
      cout<<"before fill"<<endl;
    outputTree->Fill();
    if(LDEBUG) 
      cout<<"after fill"<<endl;
    // flush the output tree
    if(outputEvents==100) {
      outputEvents=0;
      outputTree->Write(0,TObject::kWriteDelete);
    }

  }

  // count selected events
  countersred[3]++;

  //
  // analysis step
  //
  if( typerun == kFill || typerun == kFillReduce ) {
    if(makeDummyTrees) return hasoutputfile;
    // event selection
    for (size_t i=0; i<analyses.size(); i++) {
      if( ! analyses[i]->SelectEvents(*this, jentry) ) {
	return hasoutputfile;
      }
    }
    // final analysis
    for (size_t i=0; i<analyses.size(); i++) {
      if (analyses[i]->Analysis(*this, jentry)) { 
      	FillTreeContainer();
      }
    }
  }
  
  return hasoutputfile;
}

// ------------------------------------------------------------------------------------
void LoopAll::GetBranches(std::map<std::string,int> & names, std::set<TBranch *> & branches)
{
  for(std::map<std::string,int>::iterator it=names.begin(); it!=names.end(); ++it ) {
    const std::string & name = (*it).first;
    int typ = (*it).second;
    branch_info_t & info = branchDict[ name ];
    if( info.branch == 0 ) {
      std::cerr << "no branch '"<< name << "'" << std::endl;
      assert( 0 );
    }
    if ( itype[current]==0 && typ==1 || itype[current]!=0 && typ == 2 ){ continue; }
    /// if ( itype[current]!=0 || typ!=1 ){
    *(info.branch) = fChain->GetBranch( name.c_str() );
    if (*(info.branch) == NULL)
	    cerr << "WARNING: in LoopAll::GetBranches(..): got null pointer for branch '" << name << "', typ=" << typ << endl;
    
    branches.insert( *(info.branch) );
    // }
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::SetBranchAddresses(std::map<std::string,int> & names) {
  for(std::map<std::string,int>::iterator it=names.begin(); it!=names.end(); ++it ) {
    const std::string & name = (*it).first;
    int typ = (*it).second;
    branch_info_t & info = branchDict[ name ];
    if( info.read == 0 ){
      std::cerr << "no read function for branch '"<< name << "'" << std::endl;
      assert( 0 );
    }
    if ( itype[current]==0 && typ==1 || itype[current]!=0 && typ == 2 ){ continue; }
    /// if ( itype[current]!=0 || typ!=1)
    (this->*(info.read)) (fChain);
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::Branches(std::list<std::string> & names) {
  for(std::list<std::string>::iterator it=names.begin(); it!=names.end(); ++it ) {
    std::cerr << __FILE__ << ":" << __LINE__ << " " << *it << " " << outputTree << std::endl;  
    const std::string & name = *it;
    branch_info_t & info = branchDict[ name ];
    if( info.write == 0 ){
      std::cerr << "ERROR: no write function for branch '"<< name << "', "
		<< "check the branch names in the file specified with the 'outputBranches' datacard (e.g. reduction_output.dat)" << std::endl;

      assert( 0 );
    }
    (this->*(info.write)) (outputTree);
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::GetEntry(std::set<TBranch *> & branches, int jentry)
{
    for(std::set<TBranch *>::iterator it=branches.begin(); it!=branches.end(); ++it ) {
	if(LDEBUG){
	    std::cout<<"getting entry:  "<<(*it)->GetName()<<std::endl;
	}
    if( (*it)->GetReadEntry() != jentry ) {  (*it)->GetEntry(jentry); }
  }
}
// ------------------------------------------------------------------------------------
void LoopAll::BookTreeBranch(std::string name, int type, std::string dirName){
  for(unsigned int ind=0; ind<treeContainer[dirName].size(); ind++) {
    treeContainer[dirName][ind].AddTreeBranch(name,type);
  }
}
// ------------------------------------------------------------------------------------
void LoopAll::BookHisto(int h2d,
			int typplot,
			int typeplotall,
			int histoncat,
			int nbinsx,
			int nbinsy,
			float lowlim,
			float highlim,
			float lowlim2,
			float highlim2,
			const char *name,
			const char *xaxis, 
			const char* yaxis) {

  for(unsigned int ind=0; ind<histoContainer.size(); ind++) {
    if (nbinsy == 0)
      histoContainer[ind].Add(const_cast<char*>(name), const_cast<char*>(xaxis), const_cast<char*>(yaxis), histoncat, nbinsx, lowlim, highlim);
    if (nbinsy != 0)
      histoContainer[ind].Add(const_cast<char*>(name), const_cast<char*>(xaxis), const_cast<char*>(yaxis), histoncat, nbinsx, lowlim, highlim, nbinsy, lowlim2, highlim2);
  }
}

// ------------------------------------------------------------------------------------
void LoopAll::AddCut(char *cutnamesc, int ncatstmp, int ifromright, int ifinalcut, float *cutValuel, float *cutValueh) {

  if(LDEBUG) cout<<"InitCuts START"<<endl;
  
  Cut* this_cut = new Cut();
  this_cut->name = cutnamesc;
  this_cut->fromright = ifromright;
  this_cut->finalcut = ifinalcut;
  this_cut->ncat = ncatstmp;
  this_cut->cut.clear();
  this_cut->cutintervall.clear();
  this_cut->cutintervalh.clear();
  this_cut->mycutvar=0;
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

// ------------------------------------------------------------------------------------
void LoopAll::InitCounters(){
  if(LDEBUG) cout<<"InitCounts BEGIN"<<endl;

  for(unsigned int i=0; i<sampleContainer.size(); i++)
    counterContainer.push_back(CounterContainer(i));

  if(LDEBUG) cout<<"InitCounts END"<<endl;
}

/////// ------------------------------------------------------------------------------------
/////void LoopAll::AddCounter(int countersncat,
/////			 char *countername,
/////			 char *denomname0,
/////			 char *denomname1,
/////			 char *denomname2) {
/////    
/////  std::string* counternames_str = new std::string(countername);
/////  if(LDEBUG) cout<<" counternames_str"<<counternames_str<< endl; 
/////  //counternames_str.assign(counternames);
/////  if(LDEBUG) cout<<" *counternames_str"<<*counternames_str<< endl; 
/////  
/////  std::string* denomname0_str = new std::string(denomname0);
/////  if(LDEBUG) cout<<" denomname0_str"<<denomname0_str<< endl; 
/////  
/////  //std::string denomname1_str;
/////  std::string* denomname1_str =new std::string(denomname1);
/////  //denomname1_str.assign(denomname1);
/////  if(LDEBUG) cout<<" denomname1_str"<<denomname1_str<< endl; 
/////  
/////  std::string* denomname2_str = new std::string(denomname2);
/////  if(LDEBUG) cout<<" denomname2_str"<<denomname2_str<< endl; 
/////  
/////  for(unsigned int ind=0; ind<sampleContainer.size(); ind++) {
/////    if(LDEBUG) cout<<"adding to "<<ind<<" sampleContainer"<< endl; 
/////    counterContainer[ind].Add(*counternames_str
/////			      ,countersncat
/////			      ,*denomname0_str
/////			      ,*denomname1_str
/////			      ,*denomname2_str);
/////    if(LDEBUG) cout<<"added to "<<ind<<" sampleContainer"<< endl; 
/////  }
/////}


void LoopAll::AddCounter(int countersncat,
			 const char *countername,
			 const char *denomname0,
			 const char *denomname1,
			 const char *denomname2) {

  std::string* counternames_str = new std::string(countername);
  if(LDEBUG) cout<<" counternames_str"<<counternames_str<< endl; 

  if(LDEBUG) cout<<" *counternames_str"<<*counternames_str<< endl; 
  
  std::string* denomname0_str = new std::string(denomname0);
  if(LDEBUG) cout<<" denomname0_str"<<denomname0_str<< endl; 
  

  std::string* denomname1_str =new std::string(denomname1);

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
}



 
// ------------------------------------------------------------------------------------
int LoopAll::ApplyCut(int icut, float var, int icat) {
  
  //returns 0 if not initialized
  //MARCO if(cutContainer[icut].useit==0) return 1;
  if(cutContainer[icut].ncat<2)icat=0;
  
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

// ------------------------------------------------------------------------------------
int LoopAll::ApplyCut(std::string cutname, float var, int icat) {
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name == cutname) {
      return ApplyCut(i, var, icat);
    }
  }

  std::cout<<"ApplyCut: attention cutname "<<cutname<<" not found"<<std::endl;
  return 0;
}
// ------------------------------------------------------------------------------------
void LoopAll::FillTreeContainer(std::string dir){	// To Be Called after each jentry
  if( ! dir.empty() ) {
    treeContainer[dir][current_sample_index].FillTree();
    return;
  }
  for (std::map<std::string, std::vector<TreeContainer> >::iterator ii = treeContainer.begin(); ii!=treeContainer.end(); ++ii) {
    (((*ii).second)[current_sample_index]).FillTree();
  }
}
// ------------------------------------------------------------------------------------
void LoopAll::FillTree(std::string name, float x, std::string dirName){
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillFloat(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}
// ------------------------------------------------------------------------------------
void LoopAll::FillTree(std::string name, double x, std::string dirName){
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillDouble(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}
// ------------------------------------------------------------------------------------
void LoopAll::FillTree(std::string name, int x, std::string dirName){
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillInt(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}
// ------------------------------------------------------------------------------------
void LoopAll::FillTree(std::string name, unsigned int x, std::string dirName){
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillUInt(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}

void LoopAll::FillTree(std::string name, std::string x, std::string dirName) {
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillString(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}
void LoopAll::FillTree(std::string name, bool x, std::string dirName) {
  if (treeContainer.find(dirName) != treeContainer.end())
    treeContainer[dirName][current_sample_index].FillBool(name, x);
  else
    std::cout << "Tree type " << dirName << " not defined" << std::endl;
}
 
// ------------------------------------------------------------------------------------
void LoopAll::FillHist(std::string name, float y) {
  FillHist(name, 0, y);
}

// ------------------------------------------------------------------------------------
void LoopAll::FillHist2D(std::string name, float x, float y) {
  FillHist2D(name, 0, x, y);
}

// ------------------------------------------------------------------------------------
void LoopAll::FillHist(std::string name, int category, float y, float wt ) {
  histoContainer[current_sample_index].Fill(name, category, y, wt);
  histoContainer.back().Fill(name, category, y, wt);
}
// ------------------------------------------------------------------------------------
void LoopAll::FillHist2D(std::string name, int category, float x, float y, float wt ) {
  histoContainer[current_sample_index].Fill2D(name, category, x, y, wt);
  histoContainer.back().Fill2D(name, category, x, y, wt);
}

// ------------------------------------------------------------------------------------
//// // ------------------------------------------------------------------------------------
//// void LoopAll::FillCounter(std::string name, float weight) {
//// 	FillCounter(name, 0, weight);
//// }

// ------------------------------------------------------------------------------------
void LoopAll::FillCounter(std::string name, float weight, int category ) 
{
  if( counterContainer.size() > current_sample_index ) {
    counterContainer[current_sample_index].Fill(name, category, weight);
  }
}

// ----------------------------------------------------------------------------------------------------------------------
bool LoopAll::CheckLumiSelection( int run, int lumi )
{
  if( (typerun == kReduce && current_sample_index > sampleContainer.size() ) ||
      ! sampleContainer[current_sample_index].hasLumiSelection ) {
    return true;
  }

  std::vector<std::pair<int,int> > &run_lumis = sampleContainer[current_sample_index].goodLumis[run];
  for(std::vector<std::pair<int,int> >::iterator it=
	run_lumis.begin(); it!=run_lumis.end(); ++it ) {
    if( lumi >= it->first && lumi <= it->second ) {
      return true;
    }
  }
  return false;
}

// ----------------------------------------------------------------------------------------------------------------------
bool LoopAll::CheckEventList( int run, int lumi, int event  )
{
  if(typerun == kReduce || ! sampleContainer[current_sample_index].hasEventList ){
    return true;
  }

  std::vector<std::pair<int,int> > &run_events = sampleContainer[current_sample_index].eventList[run];
  for(std::vector<std::pair<int,int> >::iterator it=
	run_events.begin(); it!=run_events.end(); ++it ) {
    if( lumi == it->first && event == it->second ) {
      return true;
    }
  }
  return false;
}


// ----------------------------------------------------------------------------------------------------------------------
Normalization_8TeV * LoopAll::normalizer()
{
    if( signalNormalizer == 0 ) { 
	signalNormalizer = new Normalization_8TeV();
	signalNormalizer->Init(sqrtS);
	signalNormalizer->FillSignalTypes();
    }
    return signalNormalizer;
}

// ----------------------------------------------------------------------------------------------------------------------
float LoopAll::GetCutValue(TString cutname, int icat, int highcut) {
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name == cutname) {
      if(cutContainer[i].ncat<2) icat=0;
      if(cutContainer[i].fromright==2){
        if(highcut){
          return cutContainer[i].cutintervalh[icat];
        }else{
          return cutContainer[i].cutintervall[icat];
        }
      }else{
        return cutContainer[i].cut[icat];
      }
    }
  }
  std::cout<<"GetCutValue: attention cutname "<<cutname<<" not found"<<std::endl;
  return 0.;
}


int LoopAll::ApplyCut(int icut, int icat) {

  if(cutContainer[icut].fromright==2) {
    if(*(cutContainer[icut].mycutvar)<cutContainer[icut].cutintervall[icat] || *(cutContainer[icut].mycutvar)>cutContainer[icut].cutintervalh[icat]) return 0; 
  }
  else if (cutContainer[icut].fromright==1) {
    if(*(cutContainer[icut].mycutvar)>cutContainer[icut].cut[icat]) return 0; 
  }
  else if (cutContainer[icut].fromright==0) {
    if(*(cutContainer[icut].mycutvar)<cutContainer[icut].cut[icat]) return 0; 
  }
  return 1;
}

int LoopAll::ApplyCut(TString cutname, int icat) {

  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name == cutname) {
      return ApplyCut(i, icat);
    }
  } 
  //std::cout<<"ApplyCut: attention cutname "<<cutname<<" not found"<<endl;
  return 0;
}


int LoopAll::ApplyCut(int icut, int * passcategory) {
  for (int i=0; i<cutContainer[icut].ncat; i++) {
    passcategory[i]=ApplyCut(icut, i);
  }
  return cutContainer[icut].ncat;
}

int LoopAll::ApplyCut(TString cutname, int * passcategory) { //returns the number of categories
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name==cutname) {
      return ApplyCut(i, passcategory);
    }
  }
  cout<<"ApplyCut: attention cutname "<<cutname<<" not found"<<endl;
  return 0;
}

int LoopAll::ApplyCuts(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed) {

  int ncats=0;
  int passcuts=1;
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].finalcut==cutset) {
      //if(cutContainer[i].useit) 
      {

	if(cutContainer[i].ncat>1) {
	  if(ncats==0) ncats=cutContainer[i].ncat;
	  if(cutContainer[i].ncat!=ncats) {
	    cout<<"ApplyCuts: attention inconsistent number of categories for cutset "<<cutset<<endl;
	    return 0;
	  }
	}
	ncutsapplied++;
	int icatuse=icat;
	if(cutContainer[i].ncat<=1) {
	  icatuse=0;
	}
	if(ApplyCut(i,icatuse)) {
	  ncutspassed++;
	}
	else {
	  ncutsfailed++;
	  passcuts=0;
	}
      }
    }
  }
  return passcuts;
}

int LoopAll::ApplyCutsFill(int icat, int cutset, int & ncutsapplied, int & ncutspassed,  int & ncutsfailed, float histweight, float countweight) {

  int ntmpcuts=0;
  int tmppasscut[100];
  int indexcut[100];
  for(int i=0; i<100; i++) {
    tmppasscut[i]=0;
  }

  int ncats=0;
  int passcuts=1;
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].finalcut==cutset) {
      {
	
	//cout<<"ApplyCutsFill "<<cutContainer[i].finalcut<<" "<<cutContainer[i].useit<<" "<<cutContainer[i].name<<endl;
	
	if(cutContainer[i].ncat>1) {
	  if(ncats==0) ncats=cutContainer[i].ncat;
	  if(cutContainer[i].ncat!=ncats) {
	    cout<<"ApplyCuts: attention inconsistent number of categories for cutset "<<cutset<<endl;
	    return 0;
	  }
	}
	
	ncutsapplied++;
	
	int icatuse=icat;
	if(cutContainer[i].ncat<=1) {
	  icatuse=0;
	}
	//std::cout << i << " " << icatuse << std::endl;
	//std::cout << cutContainer[i].name << std::endl;
	if(ApplyCut(i,icatuse)) {
	  ncutspassed++;
	  tmppasscut[ntmpcuts]=1;
	}
	else {
	  ncutsfailed++;
	  passcuts=0;
	}
	indexcut[ntmpcuts]=i;
	ntmpcuts++;
      }
    }
  }
  //cout<<"ApplyCutsFill "<<ncutsapplied<<" "<<ncutsfailed<<" "<<ncutspassed<<" "<<endl;

  //// //n-1 cut histograms
  if(ncutsfailed==1||ncutsfailed==0) {
    for(int i=0; i<ntmpcuts; i++) {
      if(ncutsfailed==0||tmppasscut[i]==0) {
  	FillHist(cutContainer[indexcut[i]].name+"_nminus1", icat, *(cutContainer[indexcut[i]].mycutvar), histweight);
  	FillCounter(cutContainer[indexcut[i]].name+"_nminus1", countweight, icat);
  	
      }
    } 
  }
  
  for(int i=0; i<ntmpcuts; i++) {
    FillHist(cutContainer[indexcut[i]].name+"_sequential", icat, *(cutContainer[indexcut[i]].mycutvar), histweight);
    FillCounter(cutContainer[indexcut[i]].name+"_sequential", countweight, icat);
    if(tmppasscut[i]==0) break;
  }

  //remember to divide by the first for efficiencies

  return passcuts;
}

void LoopAll::FillCutPlots(int icat, int cutset, std::string postfix, float histweight, float countweight)
{
    for (unsigned int i=0; i<cutContainer.size(); i++) {
	if(cutContainer[i].finalcut!=cutset) { continue; }
	FillHist(cutContainer[i].name+postfix, icat, *(cutContainer[i].mycutvar), histweight);
	FillCounter(cutContainer[i].name+postfix, countweight, icat);
    }

}

int LoopAll::ApplyCuts(int icat, int cutset) {
  int ncutsapplied=0;
  int ncutspassed=0;
  int ncutsfailed=0;
  return ApplyCuts(icat, cutset, ncutsapplied, ncutspassed, ncutsfailed);
}

int LoopAll::ApplyCutsFill(int icat, int cutset, float histweight, float countweight) {
  int ncutsapplied=0;
  int ncutspassed=0;
  int ncutsfailed=0;
  return ApplyCutsFill(icat, cutset, ncutsapplied, ncutspassed, ncutsfailed, histweight, countweight);
}

//DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
int LoopAll::ApplyCuts(int cutset, int * passcategory, int * ncutsapplied, int * ncutspassed,  int * ncutsfailed) { //returns the number of categories

  int ncats=0;

  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].finalcut==cutset) {
      //if(cutContainer[i].useit) 
      {
	if(cutContainer[i].ncat>1) {
	  if(ncats==0) ncats=cutContainer[i].ncat;
	  if(cutContainer[i].ncat!=ncats) {
	    cout<<"ApplyCut: attention inconsistent number of categories for cutset "<<cutset<<endl;
	    return 0;
	  }
        }
      }
    }
  }

  for (int j=0; j<ncats; j++) {
    ApplyCuts(j, cutset, ncutsapplied[j], ncutspassed[j],  ncutsfailed[j]);
  }

  return ncats;
}

//DON'T USE THE FOLLOWING ONE, IT MAY BE CONFUSING
int LoopAll::ApplyCuts(int cutset, int * passcategory) { //returns the number of categories
  int ncutsapplied[100];
  int ncutspassed[100];  
  int ncutsfailed[100];
  return ApplyCuts(cutset, passcategory,ncutsapplied, ncutspassed,ncutsfailed);
}


int LoopAll::SetCutVariables(int i, float * variables) {
  cutContainer[i].mycutvar=variables;
  //cout<<"TEST CUT "<<i<<" "<<cutContainer[i].name<<" "<<*(cutContainer[i].mycutvar)<<endl;
}

int LoopAll::SetCutVariables(TString cutname, float * variables) {
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name==cutname) {
      return SetCutVariables(i,variables);
    }
  }
  cout<<"SetCutVariables: attention cutname "<<cutname<<" not found"<<endl;
  return 0;
}



void LoopAll::myPrintCounters() {

  /*
  // number of files
  int nindfiles=mp->nindfiles;
  
  TString msName = histFileName;
  msName.ReplaceAll(".root", ".csv");

  cout<<"Print the counters for intL: "<<mp->intlumi<<endl;

  // first oocalc file ("all in a row")
  FILE *filenew;
  TString s1 = msName+"_new.csv";
  filenew=fopen(s1,"w");
  //
  
  // second oocalc file (multi rows)
  FILE *file;
  TString s2 = msName+".csv";
  file=fopen(s2, "w");
  //
  
  // third file: ascii
  FILE *fileascii;
  TString s3 = msName+".tzt";
  fileascii=fopen(s3, "w");
  //
  
  // first oocalc file ("all in a row")
  fprintf(file,"indexfiles, namefile, scale, lumi, intlum, weight, sample,");
  //
  
  // first oocalc file ("all in a row")
  for (int i=0; i<Ncounters; i++) {
  if(counterprint[i]) {
  int nCategories = countersncat[i];
  if(nCategories == 0)
  nCategories = 1;
  //
  // multiplicity of lines for the different categories
  for(int iCategory = 0; iCategory<nCategories; iCategory++)
  fprintf(file,",Name, Counted, TotalEvents, Sigma, Selection 1, Eff, Selection 2, Eff, Selection 3, Eff,");
  //
  }
  }
  fprintf(file,"\n");
  //
  
  // first oocalc file ("all in a row")
  stringstream fileLines[nindfiles];
  bool fillInit = true;
  //
  
  for (int i=0; i<Ncounters; i++) {
  if(counterprint[i]) {
      
  // second oocalc file (multi rows)
  fprintf(filenew,"Number, Sample, Counter Name, Counted, TotalEvents, Sigma, Selection 1, Eff, Selection 2, Eff, Selection 3, Eff");
  fprintf(filenew,",indexfiles, namefile, scale, lumi, intlum, weight, xsec, ntot, nred, kfac");
  fprintf(filenew,"\n");
  //
      
  // third file: ascii
  fprintf(fileascii,"#############################################\n");
  fprintf(fileascii,"Number \t Sample");
  fprintf(fileascii,"\n");
  //
      
  int ncat = countersncat[i];
  int c=-1;
  while (c<ncat) {
  if(c==-1) c=0;
  float counts;
  float denominatorCounts[3];
  TString counterNames[3];
  double counterEfficiencies[3];
	
  for (int ind=0; ind<nindfiles; ind++) {
	  
  int indexfiles=mp->histoind[ind];
  int indexinfo=mp->infoind[ind];
  double weight=mp->weightind[ind];
  //double scale=mp->intlumi / mp->lumireal[indexfiles]  * mp->scale[indexinfo];
	  
  // first oocalc file ("all in a row")
  // fill file infos only once
  if(fillInit){
  fileLines[ind]
  << indexfiles                   << ","
  << mp->files[indexinfo]         << ","
  << mp->scale[indexinfo]         << ","
  << mp->lumireal[indexfiles]     << ","
  << mp->intlumi                  << ","
  << weight                       << ","
  << mp->filesshortnam[indexinfo] << ",";
  }
  //
  	  
  if (ncat==0) {
  counts=counters[indexfiles][i];
  }
  else {
  counts=counterscat[indexfiles][counterscatind[i]][c];
  }
	  
  // Three denominators for partial efficiencies
  for(unsigned int iDen=0; iDen<3; iDen++) {
  if(counterdenom[i][iDen] < 0) {
  denominatorCounts[iDen] = -1; // no counter selected
  counterNames[iDen] = "---";
  } else {
  if(countersncat[counterdenom[i][iDen]] == 0) {
  denominatorCounts[iDen]=counters[indexfiles][counterdenom[i][iDen]];
  } else if(countersncat[counterdenom[i][iDen]] != ncat) {
  // take the first category for the denominator if numerator and denominator categories do not correspond
  denominatorCounts[iDen]=counterscat[indexfiles][counterscatind[counterdenom[i][iDen]]][0];
  } else {
  denominatorCounts[iDen]=counterscat[indexfiles][counterscatind[counterdenom[i][iDen]]][c];
  }
  counterNames[iDen] = counternames[counterdenom[i][iDen]];
  }
  counterEfficiencies[iDen]= (denominatorCounts[iDen]>0 ? counts/denominatorCounts[iDen] : -1);
  }
	  
  // first oocalc file ("all in a row")
  fileLines[ind]
  << ","  << counternames[i]                        << ","  // counter Name
  << counts                                         << ","  // counter Counts
  << counts*weight/mp->scale[indexinfo]             << ","  // Total Events = counts x weight / scale
  << counts*weight/mp->scale[indexinfo]/mp->intlumi << ","  // Sigma (=xsec) = counts x weight / intlumi
  << counterNames[0].Data()                         << ","  // First denominator
  << counterEfficiencies[0]                         << ","  // Real Eff vs First denominator
  << counterNames[1].Data()                         << ","  // Second denominator
  << counterEfficiencies[1]                         << ","  // Real Eff vs Second denominator
  << counterNames[2].Data()                         << ","  // Third denominator
  << counterEfficiencies[2]                         << ","; // Real Eff vs Third denominator
  //
  // second oocalc file (multi rows)
  fprintf(filenew,"%d,%s,%s,%f,%f,%f,%s,%f,%s,%f,%s,%f",
  i, // counter Number
  mp->filesshortnam[indexinfo], // file short name
  counternames[i], // counter Name
  counts, // counter Counts
  counts*weight/mp->scale[indexinfo], // Total Events = counts x weight / scale
  counts*weight/mp->scale[indexinfo]/mp->intlumi, // Sigma (=xsec) = counts x weight / intlumi
  counterNames[0].Data(), // First denominator
  counterEfficiencies[0], // Real Eff vs First denominator
  counterNames[1].Data(), // Second denominator
  counterEfficiencies[1], // Real Eff vs Second denominator
  counterNames[2].Data(), // Third denominator
  counterEfficiencies[2]  // Real Eff vs Third denominator
  );
  fprintf(filenew,",%d,%s,%f,%f,%f,%f\n",indexfiles,mp->files[indexinfo], mp->scale[indexinfo],
  mp->lumireal[indexfiles],mp->intlumi,weight);
  //
  // third file: ascii
  fprintf(fileascii,"---------------------------------------------\n");
  fprintf(fileascii,"%d \t %s \n Counter Name: %s \n",
  i, // counter Number
  mp->filesshortnam[indexinfo], // file short name
  counternames[i] // counter Name
  );
  fprintf(fileascii,"---------------------------------------------\n");
  fprintf(fileascii,"Counted = %6.3f \n",
  counts // counter Counts
  );
  fprintf(fileascii,"  Total = %6.3f \n",
  counts*weight/mp->scale[indexinfo] // Total Events = counts x weight / scale
  );
  fprintf(fileascii,"  Sigma = %6.3f \n",
  counts*weight/mp->scale[indexinfo]/mp->intlumi // Sigma (=xsec) = counts x weight / intlumi
  );
  fprintf(fileascii,"  Scale = %6.3f \n",
  mp->scale[indexinfo] // Scale
  );
  fprintf(fileascii,"   Lumi = %6.3f \n",
  mp->lumireal[indexfiles] // Luminosity
  );
  fprintf(fileascii,"intLumi = %6.3f \n",
  mp->intlumi // Integrated Luminosity 
  );
  fprintf(fileascii," Weight = %6.3f \n",
  weight // Weight
  );
  fprintf(fileascii,"@@@\n");
  for(unsigned int iDen=0; iDen<3; iDen++) {
  if(counterNames[iDen] != "---")
  fprintf(fileascii,"Efficiency: %s / %s = %1.5f \n",
  counternames[i], // this counter
  counterNames[iDen].Data(), // i-th denominator
  counterEfficiencies[iDen]  // Real Eff vs i-th denominator
  );
  }
  }
  // second oocalc file (multi rows)
  fprintf(filenew,"\n");
  //
  // third file: ascii
  fprintf(fileascii,"\n");
  //
  // stringstream.str.c_str	
  c++;
  }
  // second oocalc file (multi rows)
  fprintf(filenew,"\n");
  //
  // third file: ascii
  fprintf(fileascii,"\n");
  //
      
  fillInit = false; 
  }// counter print if
    
  }
  
  // first oocalc file ("all in a row")
  for (int iFile=0; iFile<nindfiles; iFile++)
  fprintf(file,"%s,\n",fileLines[iFile].str().c_str());
  fclose(file);
  //
  // second oocalc file (multi rows)
  fclose(filenew);
  //
  // third file: ascii
  fclose(fileascii);
  //
  */
}

void LoopAll::myPrintCountersNew() {
  /*
  // number of files
  int nindfiles=mp->nindfiles;
  
  cout<<"Print the counters for intL: "<<mp->intlumi<<endl;


  TString a(mp->histFilNam);

  // first oocalc file ("all in a row")
  FILE *file;
  file=fopen(a+"counters_out_new.csv","w");
  //
  
  // second oocalc file (multi rows)
  FILE *filenew;
  filenew=fopen(a+"counters_out.csv","w");
  //
  
  // third file: ascii
  FILE *fileascii;
  fileascii=fopen(a+"counters_out.txt","w");
  //
  
  for (int i=0; i<Ncounters; i++) {
    
  // first oocalc file ("all in a row")
  stringstream fileLinesCatFile[nindfiles];
  stringstream fileLinesCatCounts[nindfiles];
  stringstream fileLinesCatEvents[nindfiles];
  stringstream fileLinesCatSigma[nindfiles];
  stringstream fileLinesCatEff[nindfiles][3];
  //
    
  // second oocalc file (multi rows)
  stringstream fileLinesCategories[nindfiles];
  stringstream fileLinesEff[nindfiles][3];
  stringstream fileLinesInfo[nindfiles];
  //
    
  if(i==0) {
  for (int ind=0; ind<nindfiles; ind++) {
  int indexfiles=mp->histoind[ind];
  int indexinfo=mp->infoind[ind];
  cout<<"DEBUGDEBUG mp->scale[indexinfo]  "<<ind<<" "<<indexfiles<<" "<<indexinfo<<" "<<mp->scale[indexinfo]<<endl;
  }
  }
    
  if(counterprint[i]) {
      
  int ncat = countersncat[i];
      
  // first oocalc file ("all in a row")
  fprintf(file,"Number, Sample, Counter Name, Categories");
  //
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(file,", Cat %d Counts",iCat);
  fprintf(file,", TOT Cat Counts");
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(file,", Cat %d Tot.Events",iCat);
  fprintf(file,", TOT Cat Tot.Events");
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(file,", Cat %d Sigma",iCat);
  fprintf(file,", TOT Cat Sigma");
  for(unsigned int iEff=0; iEff<3; iEff++) {
  fprintf(file,", Denominator Name");
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(file,", Cat %d Eff.",iCat);
  fprintf(file,", TOT Cat Eff.");
  }
  fprintf(file,",, indexfiles, namefile, scale, lumi, intlum, weight");
  fprintf(file,"\n");
  //
      
  // second oocalc file (multi rows)
  fprintf(filenew,"Number, Sample, Counter Name");
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(filenew,", Category, Counted, TotalEvents, Sigma");
  fprintf(filenew,",, Total Counted, ,indexfiles, namefile, scale, lumi, intlum, weight, xsec, ntot, nred, kfac");
  fprintf(filenew,"\n");
  //
  double counterNumerator_tot[nindfiles];
  double counterDenominator_tot[nindfiles][3];
  double counterEvents_tot[nindfiles];
  double counterSigma_tot[nindfiles];
  //
      
  // third file: ascii
  fprintf(fileascii,"#############################################\n");
  fprintf(fileascii,"Number \t Sample");
  fprintf(fileascii,"\n");
  //
      
  int c=-1;
  while (c<ncat) {
  if(c==-1) c=0;
  float counts;
  float denominatorCounts[3];
  TString counterNames[3];
  double counterEfficiencies[3];
	
  for (int ind=0; ind<nindfiles; ind++) {
	  
  // reset only once
  if(c==0) {
  counterNumerator_tot[ind] = 0.;
  counterEvents_tot[ind] = 0.;
  counterSigma_tot[ind] = 0.;
  }
  //
	  
  for(unsigned int iEff=0; iEff<3; iEff++)
  counterDenominator_tot[ind][iEff] = 0.;
	  
  int indexfiles=mp->histoind[ind];
  int indexinfo=mp->infoind[ind];
  double weight=mp->weightind[ind];
  //double scale=mp->intlumi / mp->lumireal[indexfiles]  * mp->scale[indexinfo];
	  

  // fill file infos only once
  if(c==0){
  // first oocalc file ("all in a row")
  fileLinesCatFile[ind]
  << i                            << "," // counter Number
  << mp->filesshortnam[indexinfo] << "," // file short name
  << counternames[i]              << "," // counter Name
  << ncat;                               // number of categories
  // second oocalc file (multi rows)
  fileLinesCategories[ind]
  << i                            << ","  // counter Number
  << mp->filesshortnam[indexinfo] << ","  // file short name
  << counternames[i]              << ","; // counter Name
  }
  // fill file infos only once
  if(c==0){
  fileLinesInfo[ind]
  << indexfiles               << ","
  << mp->files[indexinfo]     << ","
  << mp->scale[indexinfo]     << ","
  << mp->lumireal[indexfiles] << ","
  << mp->lumireal[indexinfo] << ","
  << mp->intlumi              << ","
  << weight                   << ","
  <<  mp->xsec[indexinfo]                  << ","
  <<  mp->ntot[indexinfo]                  << ","
  <<  mp->nred[indexinfo]                  << ","   
  <<  mp->kfactor[indexinfo]                
  ;
  //
  }
  //
  	  
  if (ncat==0) {
  counts=counters[indexfiles][i];
  }
  else {
  counts=counterscat[indexfiles][counterscatind[i]][c];
  }
	  
  // Three denominators for partial efficiencies
  counterNumerator_tot[ind] += counts;
  counterEvents_tot[ind]    += counts*weight;
  counterSigma_tot[ind]     += counts*weight/mp->intlumi;
  for(unsigned int iDen=0; iDen<3; iDen++) {
  if(counterdenom[i][iDen] < 0) {
  denominatorCounts[iDen] = -1; // no counter selected
  counterNames[iDen] = "---";
  } else {
  if(countersncat[counterdenom[i][iDen]] == 0) {
  denominatorCounts[iDen]=counters[indexfiles][counterdenom[i][iDen]];
  } else if(countersncat[counterdenom[i][iDen]] != ncat) {
  // take the first category for the denominator if numerator and denominator categories do not correspond
  denominatorCounts[iDen]=counterscat[indexfiles][counterscatind[counterdenom[i][iDen]]][0];
  } else {
  denominatorCounts[iDen]=counterscat[indexfiles][counterscatind[counterdenom[i][iDen]]][c];
  }
  counterNames[iDen] = counternames[counterdenom[i][iDen]];
  }
  counterEfficiencies[iDen] = (denominatorCounts[iDen]>0 ? counts/denominatorCounts[iDen] : -1);
  counterDenominator_tot[ind][iDen] += (denominatorCounts[iDen]>0 ? denominatorCounts[iDen] : 0);
  }
	  
  if(c==0){
  // first oocalc file ("all in a row")
  for(unsigned int iEff=0; iEff<3; iEff++)
  //	      if(counterNames[iEff]!="---")
  fileLinesCatEff[ind][iEff]
  << counterNames[iEff] << ","; // counter Name i-th denominator
  //
  // second oocalc file (multi rows)
  for(unsigned int iEff=0; iEff<3; iEff++)
  fileLinesEff[ind][iEff]
  << ",Efficiency " << iEff+1 << ","  // i-th Efficiency
  << counterNames[iEff]       << ","; // i-th denominator name
  //
  }
	  
  // first oocalc file ("all in a row")
  fileLinesCatCounts[ind]
  << counts << ","; // counter Counts
  fileLinesCatEvents[ind]
  << counts*weight << ","; // Total Events = counts x weight
  fileLinesCatSigma[ind]
  << counts*weight/mp->intlumi << ","; // Sigma (=xsec) = counts x weight / intlumi
  for(unsigned int iEff=0; iEff<3; iEff++)
  //	    if(counterNames[iEff]!="---")
  fileLinesCatEff[ind][iEff]
  << counterEfficiencies[iEff] << ","; // Eff i-th denominator
  //
	  
  // second oocalc file (multi rows)
  fileLinesCategories[ind]
  << c                         << ","  // category
  << counts                    << ","  // counter Counts
  << counts*weight             << ","  // Total Events = counts x weight
  << counts*weight/mp->intlumi << ","; // Sigma (=xsec) = counts x weight / intlumi
  for(unsigned int iEff=0; iEff<3; iEff++)
  fileLinesEff[ind][iEff]
  << c                         << ","  // category
  << denominatorCounts[iEff]   << ","  // denominator counter Counts
  << counterEfficiencies[iEff] << ","; // Eff vs i-th denominator
  //
	  
  // third file: ascii
  fprintf(fileascii,"---------------------------------------------\n");
  fprintf(fileascii,"%d \t %s \n Counter Name: %s \n",
  i, // counter Number
  mp->filesshortnam[indexinfo], // file short name
  counternames[i] // counter Name
  );
  fprintf(fileascii,"---------------------------------------------\n");
  fprintf(fileascii,"Counted = %6.3f \n",
  counts // counter Counts
  );
  fprintf(fileascii,"  Total = %6.3f \n",
  //counts*weight/mp->scale[indexinfo] // Total Events = counts x weight / scale
  counts*weight // Total Events = counts x weight / scale
  );
  fprintf(fileascii,"  Sigma = %6.3f \n",
  //counts*weight/mp->scale[indexinfo]/mp->intlumi // Sigma (=xsec) = counts x weight / intlumi
  counts*weight/mp->intlumi // Sigma (=xsec) = counts x weight / intlumi
  );
  float scaledummy=1.;
  fprintf(fileascii,"  Scale = %6.3f \n",
  scaledummy // Scale
  //mp->scale[indexinfo] // Scale
  );
  fprintf(fileascii,"   Lumi = %6.3f \n",
  mp->lumireal[indexfiles] // Luminosity
  );
  fprintf(fileascii,"intLumi = %6.3f \n",
  mp->intlumi // Integrated Luminosity 
  );
  fprintf(fileascii," Weight = %6.3f \n",
  weight // Weight
  );
  fprintf(fileascii,"@@@\n");
  for(unsigned int iDen=0; iDen<3; iDen++) {
  if(counterNames[iDen] != "---")
  fprintf(fileascii,"Efficiency: %s / %s = %1.5f \n",
  counternames[i], // this counter
  counterNames[iDen].Data(), // i-th denominator
  counterEfficiencies[iDen]  // Real Eff vs i-th denominator
  );
  }
	  
	  
  }
	
  // third file: ascii
  fprintf(fileascii,"\n");
  //
	
  c++;
	
  }
      
  // Loop on Files
  for (int iFile=0; iFile<nindfiles; iFile++) {
  // first oocalc file ("all in a row")
  fprintf(file,"%s, %s %f, %s %f, %s %f",
  fileLinesCatFile[iFile].str().c_str(),   // File Name and counter name + number of categories
  fileLinesCatCounts[iFile].str().c_str(), // Categories Counts
  counterNumerator_tot[iFile],             // TOT Categories Counts
  fileLinesCatEvents[iFile].str().c_str(), // Categories Events
  counterEvents_tot[iFile],                // TOT Categories Events = counts x weight
  fileLinesCatSigma[iFile].str().c_str(),  // Categories Sigma
  counterSigma_tot[iFile]                  // TOT Categories Sigma (=xsec) = counts x weight / intlumi
  );
  for(unsigned int iEff=0; iEff<3; iEff++)
  fprintf(file,", %s %f",
  fileLinesCatEff[iFile][iEff].str().c_str(), // Categories Efficiencies
  (counterDenominator_tot[iFile][iEff] !=0 ? counterNumerator_tot[iFile]/counterDenominator_tot[iFile][iEff] : -1) // overall efficiency
  );
  fprintf(file,",, %s",
  fileLinesInfo[iFile].str().c_str()
  );
  fprintf(file,"\n");
  //
	
  // second oocalc file (multi rows)
  fprintf(filenew,"%s, %f,, %s\n",
  fileLinesCategories[iFile].str().c_str(),
  counterNumerator_tot[iFile],
  fileLinesInfo[iFile].str().c_str());
  //
  fprintf(filenew,", , Denominator Name");
  for(unsigned int iCat = 0; iCat<ncat; iCat++)
  fprintf(filenew,", Category, Counted, Efficiency");
  // Total
  fprintf(filenew,",, Total Counted, Overall Efficiency, ");
  fprintf(filenew,"\n");
  //
  for(unsigned int iEff=0; iEff<3; iEff++) {
  fprintf(filenew,"%s",fileLinesEff[iFile][iEff].str().c_str());
  fprintf(filenew,",%f, %f \n",
  counterDenominator_tot[iFile][iEff], // overall counts
  (counterDenominator_tot[iFile][iEff] !=0 ? counterNumerator_tot[iFile]/counterDenominator_tot[iFile][iEff] : -1) // overall efficiency
  );
  }
  fprintf(filenew,"\n");
  //
	
  } // Loop on Files
  //

  // first oocalc file ("all in a row")
  fprintf(file,"\n\n");
  //
      
  // second oocalc file (multi rows)
  fprintf(filenew,"\n");
  //
      
  // third file: ascii
  fprintf(fileascii,"\n");
  //
      
  } // counter print if
    
  }
  
  //// first oocalc file ("all in a row")
  //for (int iFile=0; iFile<nindfiles; iFile++)
  //  fprintf(file,"%s,\n",fileLines[iFile].str().c_str());
  //fclose(file);
  ////

  // second oocalc file (multi rows)
  fclose(filenew);
  //
  
  // third file: ascii
  fclose(fileascii);
  //
  */
}

void LoopAll::AddCut2(char *cutnamesc, int ncatstmp, int ifromright, int ifinalcut, float *cutValuel, float *cutValueh, int iread, int plot, int bins, float xmin, float xmax, char* xaxis, char* yaxis) {

  AddCut(cutnamesc, ncatstmp, ifromright, ifinalcut, cutValuel, cutValueh);

  char a[100];
  sprintf(a, "%s_nminus1", cutnamesc);
  BookHisto(0, 0, 0, ncatstmp, bins, 0, xmin, xmax, 0., 0., a, xaxis, yaxis);
  AddCounter(ncatstmp, a, "", "", "");

  sprintf(a, "%s_sequential", cutnamesc);
  BookHisto(0, 0, 0, ncatstmp, bins, 0, xmin, xmax, 0, 0, a, xaxis, yaxis);
  AddCounter(ncatstmp, a, "", "", "");
}

void LoopAll::SetSubJob(bool issubjob){

  is_subjob=issubjob;
}

