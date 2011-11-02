#define LDEBUG 0 

#include "LoopAll.h"
#include "Tools.h"

#include <iostream>
#include <iterator>
#include <math.h>
#include "stdlib.h"

using namespace std;

#include "GeneralFunctions_cc.h"
#include "BaseAnalysis.h"

// ------------------------------------------------------------------------------------
BaseAnalysis* LoopAll::AddAnalysis(BaseAnalysis* baseAnalysis) {
  
  analyses.push_back(baseAnalysis);
  
  return baseAnalysis;
}

// ------------------------------------------------------------------------------------
void LoopAll::SetTypeRun(int t, const char* name) {
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
    outputFile=TFile::Open(outputFileName,"recreate");
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
    
    // book output braches
    Branches(outputBranchNames);
    for (size_t i=0; i<analyses.size(); i++) {
	    analyses[i]->ReducedOutputTree(*this,outputTree);
    }

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

  
  if (typerun == 1) {
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
  
  if(typerun == 0 || typerun == 2) {
    if(DEBUG) cout<<"typerun is " << typerun <<endl;
    hasoutputfile=0;
  } else if(typerun == 1) {
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
			    bool addnevents
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

  sampleContainer.push_back(SampleContainer());
  sampleContainer.back().itype = type;
  sampleContainer.back().ntot = ntot;
  sampleContainer.back().nred = nred;
  sampleContainer.back().histoplotit = histoplotit;
  sampleContainer.back().filesshortnam = filesshortnam;
  sampleContainer.back().lumi = lumi;
  sampleContainer.back().xsec = xsec;
  sampleContainer.back().kfactor = kfactor;
  sampleContainer.back().scale = scale;
  sampleContainer.back().computeWeight(intlumi);

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
	std::cout << "Combining Current File - " << (*it) << std::endl;

	for (std::vector<std::string>::iterator it_hist=histogramNames.begin()
	    ;it_hist!=histogramNames.end()
	    ;it_hist++) {
			
		TH1F *histExtra = (TH1F*) (*it_file)->Get(Form("th1f_%s",it_hist->c_str()));
		rooContainer->AppendTH1F(*it_hist,histExtra);	
	}
		
	RooWorkspace *work = (RooWorkspace*) (*it_file)->Get("cms_hgg_workspace");
	for (std::vector<std::string>::iterator it_data=datasetNames.begin()
	    ;it_data!=datasetNames.end()
	    ;it_data++) {

		RooDataSet *dataExtra = (RooDataSet*) work->data(Form("%s",it_data->c_str()));
		rooContainer->AppendDataSet(*it_data,dataExtra);	
	}

	std::cout << "Finished Combining File - " << (*it) << std::endl;
	delete work;
	(*it_file)->Close();
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
  
  for (;it!=files.end()
         ;it_file++,it_tree++,it_treelumi++,it_treepar++,it++){ 
    
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
      
      Init(typerun, *it_tree);
    }

    Loop(i);

    if(tot_events != 0 && (*it_tree) != 0  ) {
      (*it_tree)->Delete("");
    }

    *it_treelumi = (TTree*) (*it_file)->Get("lumi");
    if( *it_treelumi != 0 && outputFile ) {
      StoreProcessedLumis( *it_treelumi );
    }

    // EDIT - Cannot close the first file since it is in use after 
    // file 0 
    if (i>0)
      if((*it_file)->IsOpen())
        (*it_file)->Close();
    
    i++;
  }
  
  //now close the first File
  if( !Files.empty() &&  Files[0]->IsOpen())
    Files[0]->Close();

  TermReal(typerun);
  Term();
}

// ------------------------------------------------------------------------------------
void LoopAll::StoreProcessedLumis(TTree * tree){
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumis",&lumis);
  for( int ii=0; ii<tree->GetEntries(); ++ii) {
    tree->GetEntry(ii);
    if( !CheckLumiSelection(run,lumis)  ) { continue; }
    if( outputFile ) outputTreeLumi->Fill();
  }
  if(outputFile) {
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
	counters(4,0.), countersred(4,0.)
{  
#include "branchdef/newclonesarray.h"

#ifndef __CINT__
#include "branchdef/branchdict.h"
	DefineUserBranches();
#endif

  rooContainer = new RooContainer();
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
    temp.setScale(thisSample.weight);
    histoContainer.push_back(temp);
  }

  HistoContainer temp(sampleContainer.size(),"tot");
  temp.setScale(1.);
  histoContainer.push_back(temp);
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
    if(outputFile) outputTreeLumi->Write(0,TObject::kWriteDelete);
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
  for (Int_t jentry=0; jentry<nentries;jentry++) {
    
    if(jentry%10000==0) {
      cout << "Entry: "<<jentry << " / "<<nentries <<  " "  ;
      copy(countersred.begin(), countersred.end(), std::ostream_iterator<float>(cout, "_") );
      cout << endl;
    }
    
    if(LDEBUG) 
      cout<<"call LoadTree"<<endl;
    
    Int_t ientry = LoadTree(jentry);
  
    if (ientry < 0) 
      break;
    
    if(typerun != 1) {
      nb=0;
    }

    nbytes += nb;

    if(LDEBUG) 
      cout<<"Call FillandReduce "<<endl;
      
    hasoutputfile = this->FillAndReduce(jentry);
    if(LDEBUG) 
      cout<<"Called FillandReduce "<<endl;
  }

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
        
        std::vector<std::string> *parameters = new std::vector<std::string>;
        std::string *job_maker = new std::string;
        Int_t tot_events, sel_events, type, version, reductions;
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
	bookGlobalCounters( TreesPar[a], outputTreePar, globalCountersNames, globalCounters, fileGlobalCounters );
	assert( globalCountersNames.size() == globalCounters.size() && globalCounters.size() == fileGlobalCounters.size() );

        TreesPar[a]->GetEntry(0);
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
 //     } else {
  //      std::cerr << "Cannot write Parameter tree." << std::endl;
 //     }
 //   }
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
}

// ------------------------------------------------------------------------------------
void LoopAll::WriteFits() {
  
	hfile = TFile::Open(histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  hfile->cd();
  rooContainer->Save();
  hfile->Close();
}

// ------------------------------------------------------------------------------------
void LoopAll::WriteHist() {

	hfile = TFile::Open(histFileName, "RECREATE", "Globe ROOT file with histograms");

  hfile->cd();
  for(unsigned int ind=0; ind<histoContainer.size(); ind++) {
    histoContainer[ind].Save();
  }
  outputTreeLumi->Write();

  hfile->Close();
      
  if (makeOutputTree) 
    outputFile->cd();

}


// ------------------------------------------------------------------------------------
void LoopAll::WriteCounters() {
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

  TString msName = histFileName.ReplaceAll(".root", 5, ".csv", 4);
  FILE *file;
  TString s2 = msName;/// +".csv";
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
	        fileLinesInfo[j] << indexfiles << sampleContainer[j].scale << sampleContainer[j].lumireal << intlumi << sampleContainer[j].weight;
	      }

	      float counts = counterContainer[j][i][c];
	      
	      // CHECK MATTEO
	      counterNumerator_tot[j] += counts;
	      counterEvents_tot[j] += counts*weight;
	      counterSigma_tot[j] += counts*weight/intlumi;
      }
      if(LDEBUG) std::cout<<"end cat loop"<<std::endl;
    }
    if(LDEBUG) std::cout<<"end sample loop"<<std::endl;
  }
  fclose(file);
  if(LDEBUG) std::cout<<"LoopAll::myWriteCounters - END"<<std::endl;
}

// ------------------------------------------------------------------------------------
int LoopAll::FillAndReduce(int jentry) {

  int hasoutputfile = 0;

  //count all events
  countersred[0]++;
  //
  // read all inputs 
  //
  GetEntry(inputBranches, jentry);

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
    // event selection
    for (size_t i=0; i<analyses.size(); i++) {
      if( ! analyses[i]->SelectEvents(*this, jentry) ) {
    	  return hasoutputfile;
      }
    }
    // final analysis
    for (size_t i=0; i<analyses.size(); i++) {
      analyses[i]->Analysis(*this, jentry); 
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
        if ( itype[current]!=0 || typ!=1){
	  *(info.branch) = fChain->GetBranch( name.c_str() );
          if (*(info.branch) == NULL)
            cerr << "WARNING: in LoopAll::GetBranches(..): got null pointer for branch '" << name << "', typ=" << typ << endl;

	  branches.insert( *(info.branch) );
	}
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
        if ( itype[current]!=0 || typ!=1)
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
		std::cerr << "no write function for branch '"<< name << "'" << std::endl;
		assert( 0 );
	}
	(this->*(info.write)) (outputTree);
    }
}

// ------------------------------------------------------------------------------------
void LoopAll::GetEntry(std::set<TBranch *> & branches, int jentry)
{
  for(std::set<TBranch *>::iterator it=branches.begin(); it!=branches.end(); ++it ) {
    if( (*it)->GetReadEntry() != jentry ) {  (*it)->GetEntry(jentry); }
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
			char *name) {

  for(unsigned int ind=0; ind<histoContainer.size(); ind++) {
    if (nbinsy == 0)
      histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim);
    if (nbinsy != 0)
      histoContainer[ind].Add(name, histoncat, nbinsx, lowlim, highlim, nbinsy, lowlim2, highlim2);
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

// ------------------------------------------------------------------------------------
void LoopAll::AddCounter(int countersncat,
			 char *countername,
			 char *denomname0,
			 char *denomname1,
			 char *denomname2) {
    
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
}

 
// ------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------------------
int LoopAll::ApplyCut(std::string cutname, float var, int icat) {
  for (unsigned int i=0; i<cutContainer.size(); i++) {
    if(cutContainer[i].name == cutname) {
      return ApplyCut(i, var, icat);
    }
  }
  
  //std::cout<<"ApplyCut: attention cutname "<<cutname<<" not found"<<endl;
  return 0;
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
  histoContainer.back().Fill(name, category, y, wt);
}

//// // ------------------------------------------------------------------------------------
//// void LoopAll::FillCounter(std::string name, float weight) {
//// 	FillCounter(name, 0, weight);
//// }

// ------------------------------------------------------------------------------------
void LoopAll::FillCounter(std::string name, float weight, int category ) {
	counterContainer[current_sample_index].Fill(name, category, weight);
}

// ----------------------------------------------------------------------------------------------------------------------
bool LoopAll::CheckLumiSelection( int run, int lumi )
{
	if(typerun == kReduce || ! sampleContainer[current_sample_index].hasLumiSelection ){
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
