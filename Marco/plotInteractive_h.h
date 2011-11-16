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
#define MAXINDEXFILES 500
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


//=====================================================================================
//==========================  End Interactive Plotting  ===============================
//=====================================================================================
