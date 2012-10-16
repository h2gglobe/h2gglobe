#include <typeinfo>
//#include <algorithm>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TMacro.h"
#include "TObjString.h"

#include "../python/createCorrectedBackgroundModel.C"

#include "../interface/FMTPlots.h"
#include "../interface/FMTSigInterp.h"
#include "../interface/FMTSetup.h"

using namespace std;
namespace po = boost::program_options;

FMTSetup::FMTSetup(string filename):
	all_(false),
	fit_(false),
  catByHand_(false),
	rebin_(false),
	skipRebin_(false),
	justRebin_(false),
	binEdges_(false),
  dumpDatFile_(false),
	interp_(false),
	datacards_(false),
	diagnose_(false),
	web_(false),
	blinding_(true),
	runCombine_(false),
	checkHistos_(false),
	safeMode_(true),
  noPlot_(true),
	runSB_(false),
	cleaned(false)
{
  //if (filename!="0") system(Form("cp %s %s_beforeFMT.root",filename.c_str(),filename.c_str()));
  	
}

FMTSetup::~FMTSetup(){
	if (!cleaned) delete rebinner;
	cout << "Exiting..." << endl;
  //system(Form("cp %s %s_afterFMT.root",filename_.c_str(),filename_.c_str()));
  cout << "Original file " << filename_ << " backed up to " << Form("%s_beforeFMT.root",filename_.c_str()) << endl;
  cout << "Original file " << filename_ << " updated." << endl;
  cout << "Complete file copied to " << Form("%s_afterFMT.root",filename_.c_str()) << endl; 

}

void FMTSetup::OptionParser(int argc, char *argv[]){

  cout << "\033[1mFullMvaToolkit -- Developed by Matthew Kenzie and Nick Wardle \033[0m \n";
  cout <<        "                  Imperial College London\n" << endl;
  cout << "\033[1mRecommend running is with following arguments:\033\[0m \n\t\t./runIt.exe -i <filename> -b -I -d -D -w <web_dir> -c \n" << endl;

  po::options_description desc("Allowed options");

  desc.add_options()
    ("help,h",                                                          "Show this message")
    ("filename,i", po::value<string>(&filename_),                       "Input file name")
    ("outfilename,o", po::value<string>(&outfilename_),                 "Output file name")
    ("fit,f",      po::value<string>(&fitString_),                      "Fit invariant mass distribution: \n" 
                                                                        "  - \tCan take mulitiple arguments which can be split by comma (110,120), or by range (110-120), or both (110,112,115-120)\n"
                                                                        "  - \tDefault is to run all \n" 
                                                                        "  - \tWill accept any double in fit range \n"
                                                                        "  - \tNOTE: this option should not be run at a MC mass without calling rebinning. ")
    ("rebin,r",    po::value<string>(&rebinString_),                    "Run optimized rebinning: \n" 
                                                                        "  - \tCan take multiple arguments which can be split by comma (110,120), or by range (110-120), or both (110,120-150)\n" 
                                                                        "  - \tDefault is to run all \n" 
                                                                        "  - \tWill accept only integers in 5GeV steps \n" 
                                                                        "  - \tNOTE: this will re-run all fits and rebinnings around this mass. This is the recommended way of executing any refit or re-rebinning. You should opt to run on the nearest MC mass. E.g. to refit and rebin 112.5 use --rebin 115")
    ("catByHand,H",                                                     "Categorize events by hand")
    ("skipRebin,N",  																										"Skip the rebinning stage")
		("justRebin,J",																											"Just extract bin edges don't do anything else")
    ("getBinEdges,B",																										"Use bin edges from mvaanalysis")
    ("dumpDatFile,F",po::value<string>(&dumpDatFil_),                   "Save a new .dat file. For example if you want to save the bin edges so they can be read in later.")
    ("bkgModel,b",  																										"Correct the background model")
    ("interp,I",    																										"Run signal interpolation")
    ("datacards,d", 																										"Produce datacards")
    ("diagnose,D",  																										"Run full diagnostics (makes plots) - increases running time")
    ("www,w",     	po::value<string>(&webDir_),												"Publish to web - increases running time")
    ("runCombine,C", 																										"Run Higgs combine tool (please make sure you have checked out HiggsAnalyis/CombinedLimit and compile it somewhere in your CMSSW area)")
    ("mHMin,l",			po::value<int>(&tempmHMin_),													"Set lower bound (GeV) for Higgs mH")
    ("mHMax,u",			po::value<int>(&tempmHMax_),													"Set upper bound (GeV) for Higgs mH")
    ("mHStep,s",		po::value<double>(&tempmHStep_),											"Set bin size (GeV) for Higgs mH")
		("unblind,E",	  																											"Unblind analysis - data will be plotted (default is off)")
		("checkHistos,c",              																	  	"Run check on histograms in file")
    ("verbose,v",                                                       "Increase output level")
		("useDat,U", po::value<string>(&datFil_)->default_value("0"),				"Get options from .dat file not TFile")
    ("doPlot,P",                                                        "Remake plots")
    ("safeModeOff,S",                                                        "Do not run in safe mode (will delete plots)")
    ("setLumi,L",po::value<double>(&userLumi_),                         "Set luminosity (by default read from datafiles.dat)")
    ("is2011,7",                                                        "Is data run 2011 (7TeV) - default is 2012 (8TeV)")
		("runSB,R",																													"Make S/B plots - in development (please don\'t use yet)")
    ;
  
	po::positional_options_description p;
	p.add("filename",-1);

	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::store(po::command_line_parser(argc,argv).options(desc).positional(p).run(),vm);
	po::notify(vm);

  if (vm.count("help")){
    cout << desc << endl;
    exit(1);
  }
  if (!vm.count("filename")){ cerr << "WARNING -- A FILENAME MUST BE PROVIDED" << endl; exit(1);}
  if (!vm.count("outfilename")) { outfilename_=filename_+"_processed.root"; }

 	ReadRunConfig();
  // FIXME HARD CODE
  //setnLEPCategories(0);
  //setincludeLEP(false);
	fitMasses_ = getVecFromString<double>(fitString_);
  rebinMasses_ = getVecFromString<int>(rebinString_);
  organiseVectors(rebinMasses_,fitMasses_);

  if (fitMasses_.size()>0)        fit_=true;
  if (rebinMasses_.size()>0)      rebin_=true;
  if (vm.count("catByHand"))      catByHand_=true;
	if (vm.count("skipRebin")) 			skipRebin_=true;
	if (vm.count("justRebin")) 			justRebin_=true;
	if (vm.count("getBinEdges")) 		binEdges_=true;
  if (vm.count("dumpDatFile"))    dumpDatFile_=true;
	if (vm.count("bkgModel")) 			bkgModel_=true;
	if (vm.count("interp")) 				interp_=true;
	if (vm.count("datacards")) 			datacards_=true;
	if (vm.count("diagnose")) 			diagnose_=true;
	if (vm.count("unblind")) 				blinding_=false;
	if (vm.count("checkHistos")) 		checkHistos_=true;
	if (vm.count("safeModeOff")) 		safeMode_=false;
  if (vm.count("www"))            web_=true;
  if (vm.count("doPlot"))         noPlot_=false;
	if (vm.count("runSB"))					runSB_=true;
	if (vm.count("runCombine"))			runCombine_=true;
  if (vm.count("verbose"))        verbose_=true;
  else                            verbose_=false;
	if (vm.count("mHMin")) 					setmHMinimum(tempmHMin_);
	if (vm.count("mHMax")) 					setmHMaximum(tempmHMax_);
	if (vm.count("mHStep")) 				setmHStep(tempmHStep_);
	if (vm.count("getBinEdges"))		setrederiveOptimizedBinEdges(false);
	else 														setrederiveOptimizedBinEdges(true);
  if (vm.count("is2011"))         setis2011(true);
  else                            setis2011(false);
  if (vm.count("setLumi"))        setintLumi(userLumi_);

	if (fitMasses_.size()==0 && rebinMasses_.size()==0 && !skipRebin_) all_=true;

  //FIXME TEST
	if (checkHistos_) checkAllHistos();

  rebinner = new FMTRebin(filename_, getintLumi(), getis2011(), getmHMinimum(), getmHMaximum(), getmHStep(), getmassMin(), getmassMax(), getnDataBins(), getsignalRegionWidth(), getsidebandWidth(), getnumberOfSidebands(), getnumberOfSidebandsForAlgos(), getnumberOfSidebandGaps(), getmassSidebandMin(), getmassSidebandMax(), getnIncCategories(),getincludeVBF(), getnVBFCategories(), getincludeLEP(), getnLEPCategories(), getsystematics(), getrederiveOptimizedBinEdges(), getAllBinEdges(),verbose_);
  rebinner->setAllBinEdges(getAllBinEdges());
	rebinner->fitter->setblind(blinding_);
	rebinner->fitter->setplot(diagnose_);
  rebinner->setcatByHand(catByHand_);
	rebinner->setjustRebin(justRebin_);

  if (diagnose_) {
    system("mkdir -p plots/png");
    system("mkdir -p plots/pdf");
    system("mkdir -p plots/macro");
  }
	printPassedOptions();
  if (dumpDatFile_) dumpDatFile(dumpDatFil_);

	cout << "MC: [";
	printVec(getMCMasses());
	cout << "] " << endl;
  cout << "Bussy" << endl;
	cout << "MH: [";
	printVec(getAllMH());
	cout << "] " << endl;

	// clean plots
	if (diagnose_ && !safeMode_) {
		cout << "WARNING! -- Removing...." << endl;
		system("rm -rf plots/pdf/*");
		system("cp plots/png/PhoPhoDraw.png plots/");
		system("cp plots/png/PhotPhotEvent.png plots/");
		system("rm -rf plots/png/*");
		system("mv plots/PhoPhoDraw.png plots/png");
		system("mv plots/PhotPhotEvent.png plots/png");
	}

}

void FMTSetup::checkAllHistos(){
  TFile *tFile = TFile::Open(filename_.c_str());
  TList *histList = tFile->GetListOfKeys();
  for (int j=0; j<histList->GetSize(); j++){
    string name = histList->At(j)->GetName(); 
    if (name.find("th1f")!=string::npos){
      TH1F *temp = (TH1F*)tFile->Get(histList->At(j)->GetName());
      checkHisto(temp);
    }
  }
  tFile->Close();
}

void FMTSetup::ReadRunConfig(){

	TMacro *mva;
	if (datFil_=="0"){ // read from tFile
		cout << "Reading run configuration from file..." << filename_ << endl;
		TFile *inFile = TFile::Open(filename_.c_str());
		mva = (TMacro*)inFile->Get("mvaanalysis");
		inFile->Close();
	}
	else { // read from dat
		cout << "Reading run configuration from file..." << datFil_ << endl;
		mva = new TMacro(datFil_.c_str(),"mvaanalysis");
	}
	TList *list = (TList*)mva->GetListOfLines();
	bool comment=false;
  for (int l=0; l<list->GetSize(); l++){
    TObjString *line = (TObjString*)list->At(l);
    string sline = Form("%s",line->GetName());
    if (sline.find("->")!=string::npos) {
      comment=!comment;
      continue;
    }
		if (sline.find("#")!=string::npos) continue;
		if (sline.find("mHMinimum=")!=string::npos)										setmHMinimum(boost::lexical_cast<int>(getOptFromConfig<double>(sline)));
		if (sline.find("mHMaximum=")!=string::npos)										setmHMaximum(boost::lexical_cast<int>(getOptFromConfig<double>(sline)));
		if (sline.find("mHStep=")!=string::npos)										 	setmHStep(getOptFromConfig<double>(sline));
		if (sline.find("massMin=")!=string::npos)										 	setmassMin(getOptFromConfig<double>(sline));
		if (sline.find("massMax=")!=string::npos)										 	setmassMax(getOptFromConfig<double>(sline));
		if (sline.find("nDataBins=")!=string::npos)										setnDataBins(getOptFromConfig<int>(sline));
		if (sline.find("signalRegionWidth=")!=string::npos)						setsignalRegionWidth(getOptFromConfig<double>(sline));
		if (sline.find("sidebandWidth=")!=string::npos)								setsidebandWidth(getOptFromConfig<double>(sline));
		if (sline.find("numberOfSidebands=")!=string::npos)						setnumberOfSidebands(getOptFromConfig<int>(sline));
		if (sline.find("numberOfSidebandsForAlgos=")!=string::npos)		setnumberOfSidebandsForAlgos(getOptFromConfig<int>(sline));
		if (sline.find("numberOfSidebandGaps=")!=string::npos)				setnumberOfSidebandGaps(getOptFromConfig<int>(sline));
		if (sline.find("massSidebandMin=")!=string::npos)							setmassSidebandMin(getOptFromConfig<double>(sline));
		if (sline.find("massSidebandMax=")!=string::npos)							setmassSidebandMax(getOptFromConfig<double>(sline));
    if (sline.find("nInclusiveCategories=")!=string::npos)        setnIncCateogies(getOptFromConfig<int>(sline));
		if (sline.find("includeVBF=")!=string::npos)									setincludeVBF(getOptFromConfig<bool>(sline));
    else                                                          setincludeVBF(false);
    if (sline.find("nVBFCategories=")!=string::npos)              setnVBFCategories(getOptFromConfig<int>(sline));
    else                                                          setnVBFCategories(0);
		if (sline.find("includeLEP=")!=string::npos)									setincludeLEP(getOptFromConfig<bool>(sline));
    else                                                          setincludeLEP(false);
    if (sline.find("nLEPCategories_")!=string::npos)              setnLEPCategories(getOptFromConfig<int>(sline));
    else                                                          setnLEPCategories(0);

    if (sline.find("doEscaleSyst=")!=string::npos)                if (getOptFromConfig<bool>(sline)) setsystematic("E_scale"); 
    if (sline.find("doEresolSyst=")!=string::npos)                if (getOptFromConfig<bool>(sline)) setsystematic("E_res"); 
    if (sline.find("doEcorrectionSyst=")!=string::npos)           if (getOptFromConfig<bool>(sline)) setsystematic("E_corr"); 
    if (sline.find("doRegressionSyst=")!=string::npos)            if (getOptFromConfig<bool>(sline)) setsystematic("regSig"); 
    if (sline.find("doPhotonIdEffSyst=")!=string::npos)           if (getOptFromConfig<bool>(sline)) setsystematic("idEff"); 
    if (sline.find("doVtxEffSyst=")!=string::npos)                if (getOptFromConfig<bool>(sline)) setsystematic("vtxEff"); 
    if (sline.find("doTriggerEffSyst=")!=string::npos)            if (getOptFromConfig<bool>(sline)) setsystematic("triggerEff"); 
    if (sline.find("doPhotonMvaIdSyst=")!=string::npos)           if (getOptFromConfig<bool>(sline)) setsystematic("phoIdMva"); 
    if (sline.find("doR9Syst=")!=string::npos)                    if (getOptFromConfig<bool>(sline)) setsystematic("r9Eff"); 
    if (sline.find("doKFactorSyst=")!=string::npos)               if (getOptFromConfig<bool>(sline)) setsystematic("kFactor");
    if (sline.find("doPdfWeightSyst=")!=string::npos)             if (getOptFromConfig<bool>(sline)) setsystematic("pdfWeight");

    if (sline.find("rederiveOptimizedBinEdges=")!=string::npos)   setrederiveOptimizedBinEdges(getOptFromConfig<bool>(sline));
    for (int m=110; m<=150; m+=5){
      if (m==145) continue;
      if (sline.find(Form("GradBinEdges_%3d=",m))!=string::npos)  setBinEdges(m,getBinEdgesFromString(getOptFromConfig<string>(sline)));
      if (sline.find(Form("VbfBinEdges_%3d=",m))!=string::npos)   setVBFBinEdges(m,getBinEdgesFromString(getOptFromConfig<string>(sline)));
      if (sline.find(Form("LepBinEdges_%3d=",m))!=string::npos)   setLEPBinEdges(m,getBinEdgesFromString(getOptFromConfig<string>(sline)));
    }
  }
  userLumi_ = getLumiFromWorkspace();
  setintLumi(userLumi_);
	printRunOptions();
}

double FMTSetup::getLumiFromWorkspace(){
  TFile *inFile = TFile::Open(filename_.c_str());
  RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  RooRealVar *tempLumi = (RooRealVar*)inWS->var("IntLumi");
  inFile->Close();
  return (tempLumi->getVal())/1000.;
}

vector<double> FMTSetup::getBinEdgesFromString(string name){
  vector<double> result;
  while (name.find(",")!=string::npos){
    result.push_back(boost::lexical_cast<double>(name.substr(0,name.find(","))));
    name = name.substr(name.find(",")+1,string::npos);
  }
  result.push_back(boost::lexical_cast<double>(name));
  return result;
}


void FMTSetup::organiseVectors(vector<int> &rebinVec, vector<double> &fitVec){
	
  if (fitVec.size()==0 && rebinVec.size()==0) return;
	else{
		for (vector<int>::iterator rebM = rebinVec.begin(); rebM != rebinVec.end(); rebM++){
			vector<double> mhMasses = getMHMasses(*rebM);
			for (vector<double>::iterator mhIt = mhMasses.begin(); mhIt != mhMasses.end(); mhIt++){
				if (find(fitVec.begin(), fitVec.end(), *mhIt) != fitVec.end()) fitVec.erase(find(fitVec.begin(),fitVec.end(),*mhIt)); 
			}
		}
		sort(fitVec.begin(),fitVec.end());
		sort(rebinVec.begin(),rebinVec.end());
	}
}

void FMTSetup::printPassedOptions(){
  cout << "Running with following options:" << endl;
  cout << "\tFile:         " << filename_ << endl;
  cout << "\tFit masses:            [";
  if (fitMasses_.size()>0) {
    for (int i=0; i<fitMasses_.size()-1; i++) cout << fitMasses_[i] << ",";
    cout << fitMasses_[fitMasses_.size()-1];
  }
  cout << "]" << endl;
  cout << "\tRebin masses:          [";
  if (rebinMasses_.size()>0) {
    for (int i=0; i<rebinMasses_.size()-1; i++) cout << rebinMasses_[i] << ",";
    cout << rebinMasses_[rebinMasses_.size()-1];
  }
  cout << "]" << endl;
  if (!fit_ && !rebin_) cout << "\tNo specifics given so running all masses" << endl;
  cout << "\tFit:                    " << fit_ << endl;
  cout << "\tRebin:                  " << rebin_ << endl;
  cout << "\tSkipRebin:              " << skipRebin_ << endl;
  cout << "\tBinEdges:               " << binEdges_ << endl;
  cout << "\tBkgModel                " << bkgModel_ << endl;
  cout << "\tInterpolate:            " << interp_ << endl;
  cout << "\tDatacards:              " << datacards_ << endl;
  cout << "\tDiagnostics:            " << diagnose_ << endl;
  cout << "\tBlinding:               " << blinding_ << endl;
  cout << "\tIs 2011 data:           " << getis2011() << endl;
  cout << "\tWeb publish:            " << web_ << endl;
  cout << "\tCheck histos:           " << checkHistos_ << endl;
  cout << "\tDump dat file:          " << dumpDatFile_ << endl;
}

void FMTSetup::CheckRunOptions(){
	if (fitMasses_.size()>0){
		cerr << "WARNING -- FMTSetup::CheckRunOptions -- ARE YOU SURE ABOUT THESE RUNNING OPTIONS?" << endl;
		if (fitMasses_.size()==1) cerr << "You have requested an indivdual fit mass: " << fitMasses_[0] << endl;
		else if (fitMasses_.size()>1) {
			cerr << "You have requested individual fit masses: [";
			for (int i=0; i<fitMasses_.size()-1; i++) cerr << fitMasses_[i] << ",";
			cerr << fitMasses_[fitMasses_.size()-1] << "]" << endl;
		}
		cerr << "The fit at the MC masses has an effect on the binning algorithm so be careful!" << endl;
		cerr << "It is recommended that you run on the nearest binning mass instead." << endl;
		cerr << "YOU HAVE BEEN WARNED! \n";
		cerr << "Press RETURN to continue (q to quit)." << endl;
    string temp;
    temp = cin.get();
    if (temp=="q" || temp=="Q") exit(1);
	}
}

void FMTSetup::runRebinning(){

	vector<int> theMasses;
	if (!skipRebin_ && rebin_) theMasses = rebinMasses_;
	else if (!skipRebin_ && all_)	 theMasses = getMCMasses();
	else return;
	for (vector<int>::iterator rebM = theMasses.begin(); rebM != theMasses.end(); rebM++){
		cout << "Running rebinning for mass " << *rebM << endl;
		cout << "UandD: ["; printVec(getUandDMCMasses(*rebM)); cout << "]" << endl;
		cout << "mH:    ["; printVec(getMHMasses(*rebM)); cout << "]" << endl;
		rebinner->executeRebinning(*rebM);
		setAllBinEdges(rebinner->getAllBinEdges());
    if (dumpDatFile_) dumpDatFile(dumpDatFil_); 
		cout << "Done rebinning" << endl;
	}
	cout << "Done all requested rebinning" << endl;
}

void FMTSetup::runFitting(){

	if (!skipRebin_ && fit_){
		for (vector<double>::iterator fitM = fitMasses_.begin(); fitM != fitMasses_.end(); fitM++){
			cout << "Running fitting for mass " << *fitM << endl;
			rebinner->fitter->redoFit(*fitM);
		}
	}
}

void FMTSetup::cleanUp(){
	// need to call FMTRebin destructor to free up file
	delete rebinner;
	cleaned=true;
}

void FMTSetup::createCorrBkgModel(){
	if (!cleaned) cleanUp();
	if (bkgModel_){
		cout << "Running createCorrectedBackgroundModel...." << endl;
		createCorrectedBackgroundModel(filename_,getnumberOfSidebands(),getsidebandWidth(),getsignalRegionWidth(),getnumberOfSidebandGaps(),getmassSidebandMin(),getmassSidebandMax(),boost::lexical_cast<double>(getmHMinimum()),boost::lexical_cast<double>(getmHMaximum()),getmHStep(),diagnose_, blinding_);
		cout << "Finished correcting background model" << endl;
	}
}

void FMTSetup::interpolateBDT(){
	if (!cleaned) cleanUp();
	if (interp_){
    cout << "Running signal interpolation...." << endl;
    FMTSigInterp *interpolater = new FMTSigInterp(filename_, getintLumi(), getis2011(), diagnose_,false,getmHMinimum(), getmHMaximum(), getmHStep(), getmassMin(), getmassMax(), getnDataBins(), getsignalRegionWidth(), getsidebandWidth(), getnumberOfSidebands(), getnumberOfSidebandsForAlgos(), getnumberOfSidebandGaps(), getmassSidebandMin(), getmassSidebandMax(), getnIncCategories(), getincludeVBF(), getnVBFCategories(), getincludeLEP(), getnLEPCategories(), getsystematics(), getrederiveOptimizedBinEdges(), getAllBinEdges(),blinding_,verbose_);
    interpolater->runInterpolation();
		delete interpolater;
	}
}

void FMTSetup::writeDataCards(){
	if (!cleaned) cleanUp();
	if (datacards_){
		cout << "Preparing to write datacards...." << endl;
		if (getis2011()){
     // cerr << "This option isn't supported yet. You will have to do this by hand. Sorry :( " << endl;
     // exit(0);
      if (blinding_) system(Form("python python/writeBinnedMvaCard_7TeV.py -i %s -p plots --makePlot --mhLow %3d.0 --mhHigh %3d.0 --mhStep %1.1f --intLumi %1.1f --blind",filename_.c_str(),getmHMinimum(),getmHMaximum(),getmHStep(),getintLumi()));
      else system(Form("python python/writeBinnedMvaCard_7TeV.py -i %s -p plots --makePlot --mhLow %3d.0 --mhHigh %3d.0 --mhStep %1.1f --intLumi %1.1f ",filename_.c_str(),getmHMinimum(),getmHMaximum(),getmHStep(),getintLumi()));
    }
    else {
      if (blinding_) system(Form("python python/writeBinnedMvaCard.py -i %s -p plots --makePlot --mhLow %3d.0 --mhHigh %3d.0 --mhStep %1.1f --intLumi %1.1f --blind",filename_.c_str(),getmHMinimum(),getmHMaximum(),getmHStep(),getintLumi()));
      else system(Form("python python/writeBinnedMvaCard.py -i %s -p plots --makePlot --mhLow %3d.0 --mhHigh %3d.0 --mhStep %1.1f --intLumi %1.1f",filename_.c_str(),getmHMinimum(),getmHMaximum(),getmHStep(),getintLumi()));
    }
	}
}

void FMTSetup::makePlots(){
  if (!cleaned) cleanUp();
  if (diagnose_ && !noPlot_){
		if (runSB_) {
			cout << "Running S/B stuff...." << endl;
			system(Form("python python/GetFakeShapeDatacards.py -i mva-datacards-grad -o fake-shape-cards -D -C -N --mhLow %3d.0 --mhHigh %3d.0 --mhStep %1.1f",getmHMinimum(),getmHMaximum(),getmHStep()));
    }
		cout << "Making plots..." << endl;
    FMTPlots *plotter = new FMTPlots(filename_, runSB_, getintLumi(), getis2011(), getmHMinimum(), getmHMaximum(), getmHStep(), getmassMin(), getmassMax(), getnDataBins(), getsignalRegionWidth(), getsidebandWidth(), getnumberOfSidebands(), getnumberOfSidebandsForAlgos(), getnumberOfSidebandGaps(), getmassSidebandMin(), getmassSidebandMax(), getnIncCategories(), getincludeVBF(), getnVBFCategories(), getincludeLEP(), getnLEPCategories(), getsystematics(), getrederiveOptimizedBinEdges(), getAllBinEdges(),blinding_,verbose_);
    vector<double> theMasses = getAllMH();
    for (vector<double>::iterator mh = theMasses.begin(); mh != theMasses.end(); mh++){
      plotter->plotAll(*mh);
    }
    plotter->makeNormPlot();
    delete plotter;
  }
}


void FMTSetup::publishToWeb(){
	if (!cleaned) cleanUp();
	if (web_){
		cout << "Publishing to web: " << webDir_ << "/plots/png/home.html" << endl;
		if (blinding_) system(Form("python python/publish_plots.py %s --blind",filename_.c_str()));
		else system(Form("python python/publish_plots.py %s --not",filename_.c_str()));
		system(Form("rm -r %s",webDir_.c_str()));
		system(Form("mkdir %s",webDir_.c_str()));
		system(Form("cp -r plots %s",webDir_.c_str()));
	}
}

void FMTSetup::runCombine(){
	if (!cleaned) cleanUp();
	if (runCombine_){
		cout << "Running combine tool..... " << endl;
		cout << Form("./python/limit.sh mva-datacards-grad/ grad $PWD") << endl;
		system(Form("./python/limit.sh mva-datacards-grad/ grad $PWD"));
	}
}
