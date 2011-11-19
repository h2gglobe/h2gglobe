//=====================================================================================
//============================   Interactive Plotting  ================================
//=====================================================================================

void myPlotInteractive(TString inputfiles="");
void myPlotInteractiveSetup(TString inputfiles="", TString anaysistag="");

void mySetTDRStyle(void);

std::pair<int,int> myGetInput(std::pair<int,int> ivar);

Int_t NFILES;
Int_t NIND;
std::pair<int,int> prevPair;
TString OPTIONS;
TString outname;
TString outnamesob;
TCanvas *ch;
TCanvas *ch2;
TCanvas *ch3;
TCanvas *ch4;
TCanvas *ch5;
THStack *myStack;
THStack *mySigStack;
THStack *myDataStack;
#define MYNINDS 50
//#define MAXINDEXFILES 500
int omitdefault[MAXINDEXFILES];//MYNINDS
int omitInd[MAXINDEXFILES];//MYNINDS
bool omitFile[MAXINDEXFILES];//MYNINDS
bool allSgnl[MAXINDEXFILES];//MYNINDS
bool mySgnl[MAXINDEXFILES];//MYNINDS
bool allBkgd[MAXINDEXFILES];//MYNINDS
bool myBkgd[MAXINDEXFILES];//MYNINDS
Int_t myColor[MAXINDEXFILES];//MYNINDS
Int_t myMarker[MAXINDEXFILES];//MYNINDS
Int_t myStackOrder[MAXINDEXFILES];//MYNINDS
Float_t myScale[MAXINDEXFILES];//MYNINDS
Float_t sigScale;
TString myLabel[MAXINDEXFILES];//MYNINDS
TString myLabeldisplay[MAXINDEXFILES];//MYNINDS
TString myHistname[MAXINDEXFILES];//MYNINDS
TString mySampleString[MAXINDEXFILES];//MYNINDS
Int_t NColsLegend;
Int_t NStackLines;
bool histStackLineIsSig[MAXINDEXFILES];
TH1F histStackLine[MAXINDEXFILES];
Int_t histStackLineColor[MAXINDEXFILES];
Float_t sethistmax, sethistmin;
Float_t textx1,textx2,texty1,texty2;
Float_t linex;
Float_t legendx1,legendx2,legendy1,legendy2;
bool RmTitles;
bool DoLog;
bool DoLegend;
bool DoBigLegend;
bool DoLegendT;
bool DoLegendR;
bool DoGridX;
bool DoGridY;
bool DoTitle;
bool DoXTitle;
bool DoYTitle;
bool DoFill;
bool DoFillSig;
bool DoReSize;
bool DoRebin;
bool DoSnB;
bool DoSvB;
bool DoSvN;
bool DoCats;
bool DoMulti;
bool DoStack;
bool DoWriteAll;
bool DoPopSig;
bool DoData;
bool DoPrintFlows;
bool DoPrintBins;
bool DoStats;
bool DoScale;
bool DoSumw2;
bool DoFitSig;
bool DoFit;
bool DoInt;
bool DoRevInt;
bool DoRedraw;
bool DoOverFlow;
bool DoUnderFlow;

Width_t lineWidth;
Int_t BigLegendInd;
Int_t nrows;
Int_t ncolumns;
Int_t NReBin;
Int_t fitXmin;
Int_t fitXmax;
Int_t xSizeMAX;
Int_t ySizeMAX;
Int_t singlexSize;
Int_t singleySize;
Int_t xSize;
Int_t ySize;
Float_t sobmax;
Int_t StartCat;
TString allext;

TString num_to_string(int myint);


//Should be in other


TString varnames[MAXVAR];
const char * varnamesc[MAXVAR];
char * varnamescread[MAXVAR];
char * xaxislabel[MAXVAR];
char * yaxislabel[MAXVAR];
Float_t highlim[MAXVAR];
Float_t lowlim[MAXVAR];
Float_t *varlinks[MAXVAR];
Int_t doplot[MAXVAR];
Int_t h2d[MAXVAR];
Int_t nbinsx[MAXVAR];
Int_t nbinsy[MAXVAR];
Float_t highlim2[MAXVAR];
Float_t lowlim2[MAXVAR];
Int_t cutindex[MAXVAR];
Int_t Nvar;
Int_t Ncatvar;

Int_t typplotall;
Int_t typplot[MAXVAR];
Int_t histoncat[MAXVAR];

Int_t histoncatindtonames[MAXVAR];

TClonesArray * tca_xaxislabels;
TClonesArray * tca_yaxislabels;
TClonesArray * tca_plotvarnames;
TClonesArray * tca_plotvarcatnames;
TClonesArray * tca_histfilename;
TClonesArray * tca_inshortnames;
TClonesArray * tca_infilenames;

Int_t Nvarcats;
Int_t catid[MAXVAR];
Int_t ncats[MAXVAR];
char * catnames[MAXVAR][MAXCAT];


TTree* plotvariables;
TTree* inputfiles;

  TString plotvarnames[1000];
  TString xaxislabels[1000];
  TString yaxislabels[1000];
  TString plotvarcatnames[40][50];

  TBranch * b_Nvar;
  TBranch * b_typplotall;
  TBranch * b_doplot;
  TBranch * b_h2d;
  TBranch * b_typplot;
  TBranch * b_histoncat;
  TBranch * b_histoncatindtonames;
  TBranch * b_nbinsx;
  TBranch * b_nbinsy;
  TBranch * b_lowlim;
  TBranch * b_highlim;
  TBranch * b_lowlim2;
  TBranch * b_highlim2;
  TBranch * b_xaxislabels;
  TBranch * b_yaxislabels;
  TBranch * b_plotvarnames;
  TBranch * b_plotvarcatnames;
  TBranch * b_Nvarcats;
  TBranch * b_catid;
  TBranch * b_ncats;

//Int_t nfiles;
  Int_t nindfiles;
//Float_t intlumi;
  TString  histfilename;
  Int_t itypePI[MAXFILES];
  Int_t histoind[MAXFILES];
  Int_t infoind[MAXFILES];
  Int_t histoplotit[MAXFILES];
  Int_t ntot[MAXFILES];
  Int_t nred[MAXFILES];
  Float_t lumi[MAXFILES];
  Float_t xsec[MAXFILES];
  Float_t kfactor[MAXFILES];
  Float_t scale[MAXFILES];
  TString  inshortnames[MAXFILES];
  TString  infilenames[MAXFILES];
  TString  infilenamesnew[MAXFILES];

  TBranch * b_nfiles;
  TBranch * b_nindfiles;
  TBranch * b_intlumi;
  TBranch * b_makeOutputTree;
  TBranch * b_histfilename;
  TBranch * b_itype;
  TBranch * b_histoind;
  TBranch * b_infoind;
  TBranch * b_histoplotit;
  TBranch * b_ntot;
  TBranch * b_nred;
  TBranch * b_lumi;
  TBranch * b_xsec;
  TBranch * b_kfactor;
  TBranch * b_scale;
  TBranch * b_inshortnames;
  TBranch * b_infilenames;



//=====================================================================================
//==========================  End Interactive Plotting  ===============================
//=====================================================================================
