#include "TCanvas.h"
#include "TF1.h"
#include "RooPlot.h"
#include "RooVoigtian.h"

#include "../interface/FinalModelConstruction.h"

using namespace std;
using namespace RooFit;

FinalModelConstruction::FinalModelConstruction(RooRealVar *massVar, RooRealVar *MHvar, RooRealVar *intL, int mhLow, int mhHigh, string proc, int cat, int nIncCats, bool doSecMods, bool isCB):
  mass(massVar),
  MH(MHvar),
  intLumi(intL),
  mhLow_(mhLow),
  mhHigh_(mhHigh),
  proc_(proc),
  cat_(cat),
  nIncCats_(nIncCats),
  doSecondaryModels(doSecMods),
  isCutBased(isCB),
  verbosity_(0),
  systematicsSet_(false),
  rvFractionSet_(false)
{
  allMH_ = getAllMH();

  // load xs and br info from Normalization_8TeV
  norm = new Normalization_8TeV();
  TGraph *brGraph = norm->GetBrGraph();
  brSpline = graphToSpline("fbr",brGraph);
  
  string procs[8] = {"ggh","vbf","wzh","wh","zh","tth","gg_grav","qq_grav"};
  for (int i=0; i<8; i++){
    TGraph *xsGraph = norm->GetSigmaGraph(procs[i].c_str());
    RooSpline1D *xsSpline = graphToSpline(Form("fxs_%s",procs[i].c_str()),xsGraph);
    xsSplines.insert(pair<string,RooSpline1D*>(procs[i],xsSpline));
  }

  // const smearing terms for resolution systematic hard-coded here!!
  constSmearVals.push_back(0.0080838); // cat0
  constSmearVals.push_back(0.0083661); // cat1
  constSmearVals.push_back(0.0086261); // cat2
  constSmearVals.push_back(0.0134325); // cat3
  constSmearVals.push_back(0.0119792); // cat4 (all exc cats)
  // should pick up different values for the cut-based but not sure what they are!
}

FinalModelConstruction::~FinalModelConstruction(){}

vector<int> FinalModelConstruction::getAllMH(){
  vector<int> result;
  for (int m=mhLow_; m<=mhHigh_; m+=5){
    if (verbosity_>=1) cout << "FinalModelConstruction - Adding mass: " << m << endl;
    result.push_back(m);
  }
  return result;
}

void FinalModelConstruction::setSecondaryModelVars(RooRealVar *mh_sm, RooRealVar *deltam, RooAddition *mh_2, RooRealVar *width){
  MH_SM = mh_sm;
  DeltaM = deltam;
  MH_2 = mh_2;
  higgsDecayWidth = width;
  
  TGraph *brGraph = norm->GetBrGraph();
  brSpline_SM = graphToSpline("fbr_SM",brGraph,MH_SM);
  brSpline_2 = graphToSpline("fbr_2",brGraph,MH_2);
  brSpline_NW = graphToSpline("fbr_NW",brGraph,MH);
  string procs[8] = {"ggh","vbf","wzh","wh","zh","tth","gg_grav","qq_grav"};
  for (int i=0; i<8; i++){
    TGraph *xsGraph = norm->GetSigmaGraph(procs[i].c_str());
    RooSpline1D *xsSpline_SM = graphToSpline(Form("fxs_%s_SM",procs[i].c_str()),xsGraph,MH_SM);
    RooSpline1D *xsSpline_2 = graphToSpline(Form("fxs_%s_2",procs[i].c_str()),xsGraph,MH_2);
    RooSpline1D *xsSpline_NW = graphToSpline(Form("fxs_%s_NW",procs[i].c_str()),xsGraph,MH);
    xsSplines_SM.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_SM));
    xsSplines_2.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_2));
    xsSplines_NW.insert(pair<string,RooSpline1D*>(procs[i],xsSpline_NW));
  }
  secondaryModelVarsSet=true;
}

void FinalModelConstruction::getRvFractionFunc(string name){
  
  assert(allMH_.size()==rvDatasets.size());
  assert(allMH_.size()==wvDatasets.size());
  vector<double> mhValues, rvFracValues;
  for (unsigned int i=0; i<allMH_.size(); i++){
    int mh = allMH_[i];
    mhValues.push_back(mh);
    double rvN = rvDatasets[mh]->sumEntries();
    double wvN = wvDatasets[mh]->sumEntries();
    rvFracValues.push_back(rvN/(rvN+wvN));
  }
  rvFracFunc = new RooSpline1D(name.c_str(),name.c_str(),*MH,mhValues.size(),&(mhValues[0]),&(rvFracValues[0]));
  rvFractionSet_=true;
}

void FinalModelConstruction::setupSystematics(){
  
  int nuisCat = cat_;
  // this is a hack for correlation with previous 2011 ws
  // for legacy paper - this MUST BE UPDATED
  if (cat_>=nIncCats_) nuisCat = nIncCats_;
  
  vertexNuisance = new RooRealVar("CMS_hgg_nuissancedeltafracright","CMS_hgg_nuissancedeltafracright",1.,0.1,10.);
  vertexNuisance->setConstant(true);
  globalScale = new RooRealVar("CMS_hgg_globalscale","CMS_hgg_globalscale",0.,-5.,5.);
  globalScale->setConstant(true);
  categoryScale = new RooRealVar(Form("CMS_hgg_nuissancedeltamcat%d",nuisCat),Form("CMS_hgg_nuissancedeltamcat%d",nuisCat),0.,-5.,5.);
  categoryScale->setConstant(true);
  categorySmear = new RooConstVar(Form("CMS_hgg_constsmearcat%d",nuisCat),Form("CMS_hgg_constsmearcat%d",nuisCat),constSmearVals[nuisCat]);
  categoryResolution = new RooRealVar(Form("CMS_hgg_nuissancedeltasmearcat%d",nuisCat),Form("CMS_hgg_nuissancedeltasmearcat%d",nuisCat),0.0,-0.2,0.2);
  systematicsSet_=true;
}

void FinalModelConstruction::buildStdPdf(string name, int nGaussians, bool recursive){

  if (!systematicsSet_) setupSystematics();
  vector<RooAddPdf*> pdfs;
  pdfs = buildPdf(name.c_str(),nGaussians,recursive,stdSplines);
  finalPdf = pdfs[0];
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    finalPdf_SM = pdfs[1];
    finalPdf_2 = pdfs[2];
    finalPdf_NW = pdfs[3];
  }
}

void FinalModelConstruction::buildRvWvPdf(string name, int nGrv, int nGwv, bool recursive){

  if (!rvFractionSet_) getRvFractionFunc(Form("%s_rvFracFunc",name.c_str()));
  if (!systematicsSet_) setupSystematics();
  RooFormulaVar *rvFraction = new RooFormulaVar(Form("%s_rvFrac",name.c_str()),Form("%s_rvFrac",name.c_str()),"TMath::Min(@0*@1,1.0)",RooArgList(*vertexNuisance,*rvFracFunc));
  vector<RooAddPdf*> rvPdfs = buildPdf(name,nGrv,recursive,rvSplines,"_rv");
  vector<RooAddPdf*> wvPdfs = buildPdf(name,nGwv,recursive,wvSplines,"_wv");
  finalPdf = new RooAddPdf(Form("%s_%s_cat%d",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[0],*wvPdfs[0]),RooArgList(*rvFraction));
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    finalPdf_SM = new RooAddPdf(Form("%s_%s_cat%d_SM",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_SM",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[1],*wvPdfs[1]),RooArgList(*rvFraction));
    finalPdf_2 = new RooAddPdf(Form("%s_%s_cat%d_2",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_2",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[2],*wvPdfs[2]),RooArgList(*rvFraction));
    finalPdf_NW = new RooAddPdf(Form("%s_%s_cat%d_NW",name.c_str(),proc_.c_str(),cat_),Form("%s_%s_cat%d_NW",name.c_str(),proc_.c_str(),cat_),RooArgList(*rvPdfs[3],*wvPdfs[3]),RooArgList(*rvFraction));
  }
}

vector<RooAddPdf*> FinalModelConstruction::buildPdf(string name, int nGaussians, bool recursive, map<string,RooSpline1D*> splines, string add){
  
  vector<RooAddPdf*> result;

  RooArgList *gaussians = new RooArgList();
  RooArgList *coeffs = new RooArgList();
  string ext = Form("%s_cat%d%s",proc_.c_str(),cat_,add.c_str());
  
  // for SM Higgs as Background
  RooArgList *gaussians_SM = new RooArgList();
  RooArgList *coeffs_SM = new RooArgList();
  
  // for Second Higgs
  RooArgList *gaussians_2 = new RooArgList();
  RooArgList *coeffs_2 = new RooArgList();

  // for Natural Width
  RooArgList *voigtians_NW = new RooArgList();
  RooArgList *coeffs_NW = new RooArgList();

  for (int g=0; g<nGaussians; g++){
    RooAbsReal *dm = splines[Form("dm_g%d",g)];
    dm->SetName(Form("dm_g%d_%s",g,ext.c_str()));
    RooAbsReal *mean = new RooFormulaVar(Form("mean_g%d_%s",g,ext.c_str()),Form("mean_g%d_%s",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH,*dm,*globalScale,*categoryScale));
    RooAbsReal *sig_fit = splines[Form("sigma_g%d",g)];
    sig_fit->SetName(Form("sigma_g%d_%s",g,ext.c_str()));
    RooAbsReal *sigma = new RooFormulaVar(Form("sig_g%d_%s",g,ext.c_str()),Form("sig_g%d_%s",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fit,*MH,*categorySmear,*categoryResolution));
    RooGaussian *gaus = new RooGaussian(Form("gaus_g%d_%s",g,ext.c_str()),Form("gaus_g%d_%s",g,ext.c_str()),*mass,*mean,*sigma);
    gaussians->add(*gaus);
    // add secondary models as well
    if (doSecondaryModels){
      assert(secondaryModelVarsSet);
      // sm higgs as background
      RooAbsReal *dmSM = splines[Form("dm_g%d_SM",g)];
      dmSM->SetName(Form("dm_g%d_%s_SM",g,ext.c_str()));
      RooAbsReal *meanSM = new RooFormulaVar(Form("mean_g%d_%s_SM",g,ext.c_str()),Form("mean_g%d_%s_SM",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH_SM,*dmSM,*globalScale,*categoryScale));
      RooAbsReal *sig_fitSM = splines[Form("sigma_g%d_SM",g)];
      sig_fitSM->SetName(Form("sigma_g%d_%s_SM",g,ext.c_str()));
      RooAbsReal *sigmaSM = new RooFormulaVar(Form("sig_g%d_%s_SM",g,ext.c_str()),Form("sig_g%d_%s_SM",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fitSM,*MH_SM,*categorySmear,*categoryResolution));
      RooGaussian *gausSM = new RooGaussian(Form("gaus_g%d_%s_SM",g,ext.c_str()),Form("gaus_g%d_%s_SM",g,ext.c_str()),*mass,*meanSM,*sigmaSM);
      gaussians_SM->add(*gausSM);
      // second degen higgs
      RooAbsReal *dm2 = splines[Form("dm_g%d_2",g)];
      dm2->SetName(Form("dm_g%d_%s_2",g,ext.c_str()));
      RooAbsReal *mean2 = new RooFormulaVar(Form("mean_g%d_%s_2",g,ext.c_str()),Form("mean_g%d_%s_2",g,ext.c_str()),"@0+@1+@0*(@2+@3)",RooArgList(*MH_2,*dm2,*globalScale,*categoryScale));
      RooAbsReal *sig_fit2 = splines[Form("sigma_g%d_2",g)];
      sig_fit2->SetName(Form("sigma_g%d_%s_2",g,ext.c_str()));
      RooAbsReal *sigma2 = new RooFormulaVar(Form("sig_g%d_%s_2",g,ext.c_str()),Form("sig_g%d_%s_2",g,ext.c_str()),"(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2))>0. ? TMath::Sqrt(TMath::Power(@0,2)-TMath::Power(@1*@2,2)+TMath::Power(@1*(@2+@3),2)) : 0.",RooArgList(*sig_fit2,*MH_2,*categorySmear,*categoryResolution));
      RooGaussian *gaus2 = new RooGaussian(Form("gaus_g%d_%s_2",g,ext.c_str()),Form("gaus_g%d_%s_2",g,ext.c_str()),*mass,*mean2,*sigma2);
      gaussians_2->add(*gaus2);
      // natural width
      RooVoigtian *voigNW = new RooVoigtian(Form("voig_g%d_%s_NW",g,ext.c_str()),Form("voig_g%d_%s_NW",g,ext.c_str()),*mass,*mean,*higgsDecayWidth,*sigma);
      voigtians_NW->add(*voigNW);
    }
    if (g<nGaussians-1) {
      RooAbsReal *frac = splines[Form("frac_g%d",g)];
      frac->SetName(Form("frac_g%d_%s",g,ext.c_str()));
      coeffs->add(*frac);
      // add secondary models as well
      if (doSecondaryModels){
        assert(secondaryModelVarsSet);
        // sm higgs as background
        RooAbsReal *fracSM = splines[Form("frac_g%d_SM",g)];
        fracSM->SetName(Form("frac_g%d_%s_SM",g,ext.c_str()));
        coeffs_SM->add(*fracSM);
        // second degen higgs
        RooAbsReal *frac2 = splines[Form("frac_g%d_2",g)];
        frac2->SetName(Form("frac_g%d_%s_2",g,ext.c_str()));
        coeffs_2->add(*frac2);
        // natural width
        coeffs_NW->add(*frac);
      }
    }
  }
  assert(gaussians->getSize()==nGaussians && coeffs->getSize()==nGaussians-1);
  RooAddPdf *pdf = new RooAddPdf(Form("%s_%s",name.c_str(),ext.c_str()),Form("%s_%s",name.c_str(),ext.c_str()),*gaussians,*coeffs,recursive);
  result.push_back(pdf);
  
  // add secondary models as well
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    // sm higgs as background
    RooAddPdf *pdf_SM = new RooAddPdf(Form("%s_%s_SM",name.c_str(),ext.c_str()),Form("%s_%s_SM",name.c_str(),ext.c_str()),*gaussians_SM,*coeffs_SM,recursive);
    result.push_back(pdf_SM);
    // second degen higgs
    RooAddPdf *pdf_2 = new RooAddPdf(Form("%s_%s_2",name.c_str(),ext.c_str()),Form("%s_%s_2",name.c_str(),ext.c_str()),*gaussians_2,*coeffs_2,recursive);
    result.push_back(pdf_2);
    // natural width
    RooAddPdf *pdf_NW = new RooAddPdf(Form("%s_%s_NW",name.c_str(),ext.c_str()),Form("%s_%s_NW",name.c_str(),ext.c_str()),*voigtians_NW,*coeffs_NW,recursive);
    result.push_back(pdf_NW);
  }

  return result;
}

void FinalModelConstruction::setRVsplines(map<string,RooSpline1D*> splines){
  rvSplines = splines;
}

void FinalModelConstruction::setWVsplines(map<string,RooSpline1D*> splines){
  wvSplines = splines;
}

void FinalModelConstruction::setSTDsplines(map<string,RooSpline1D*> splines){
  stdSplines = splines;
}

void FinalModelConstruction::setRVdatasets(map<int,RooDataSet*> data){
  rvDatasets = data;
}

void FinalModelConstruction::setWVdatasets(map<int,RooDataSet*> data){
  wvDatasets = data;
}

void FinalModelConstruction::setSTDdatasets(map<int,RooDataSet*> data){
  stdDatasets = data;
}

void FinalModelConstruction::plotPdf(string outDir){
  system(Form("mkdir -p %s",outDir.c_str()));
  
  TCanvas *canv = new TCanvas();
  RooPlot *dataPlot = mass->frame(Range(100,160));
  for (unsigned int i=0; i<allMH_.size(); i++){
    int mh=allMH_[i];
    stdDatasets[mh]->plotOn(dataPlot,Binning(80));
    MH->setVal(mh);
    extendPdf->plotOn(dataPlot);
  }
  dataPlot->Draw();
  canv->Print(Form("%s/%s_cat%d_fits.pdf",outDir.c_str(),proc_.c_str(),cat_));
  
  RooPlot *pdfPlot = mass->frame(Range(100,160));
  for (int mh=mhLow_; mh<=mhHigh_; mh++){
    MH->setVal(mh);
    extendPdf->plotOn(pdfPlot,Normalization(1.0,RooAbsReal::RelativeExpected));
  }
  pdfPlot->Draw();
  canv->Print(Form("%s/%s_cat%d_interp.pdf",outDir.c_str(),proc_.c_str(),cat_));
  delete canv;

}

RooSpline1D* FinalModelConstruction::graphToSpline(string name, TGraph *graph){
  
  vector<double> xValues, yValues;
  for (double mh=mhLow_; mh<(mhHigh_+0.25); mh+=0.5){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*MH,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

RooSpline1D* FinalModelConstruction::graphToSpline(string name, TGraph *graph, RooAbsReal *var){
  
  vector<double> xValues, yValues;
  for (double mh=mhLow_; mh<(mhHigh_+0.25); mh+=0.5){
    xValues.push_back(mh);
    yValues.push_back(graph->Eval(mh));
  }
  RooSpline1D *res = new RooSpline1D(name.c_str(),name.c_str(),*var,xValues.size(),&(xValues[0]),&(yValues[0]));
  return res;
}

void FinalModelConstruction::getNormalization(){
 
  TGraph *temp = new TGraph();
  TF1 *pol2 = new TF1("pol","pol2",110,150);
  for (unsigned int i=0; i<allMH_.size(); i++){
    double mh = double(allMH_[i]);
    RooDataSet *data = stdDatasets[mh];
    // calcu eA as sumEntries / totalxs * totalbr * intL
    double effAcc = (data->sumEntries()/(intLumi->getVal()*norm->GetXsection(mh,proc_)*norm->GetBR(mh)));
    temp->SetPoint(i,mh,effAcc);
  }
  verbosity_ >=2 ?
    temp->Fit(pol2,"EMFEX0") :
    temp->Fit(pol2,"QEMFEX0")
  ;
  TGraph *eaGraph = new TGraph(pol2);
  RooSpline1D *eaSpline = graphToSpline(Form("fea_%s_cat%d",proc_.c_str(),cat_),eaGraph);
  RooSpline1D *xs = xsSplines[proc_];
  finalNorm = new RooFormulaVar(Form("%s_norm",finalPdf->GetName()),Form("%s_norm",finalPdf->GetName()),"@0*@1*@2",RooArgList(*xs,*brSpline,*eaSpline));
  // these are for plotting
  finalNormThisLum = new RooFormulaVar(Form("%s_normThisLumi",finalPdf->GetName()),Form("%s_normThisLumi",finalPdf->GetName()),"@0*@1*@2*@3",RooArgList(*xs,*brSpline,*eaSpline,*intLumi));
  extendPdfRel = new RooExtendPdf(Form("extend%s",finalPdf->GetName()),Form("extend%s",finalPdf->GetName()),*finalPdf,*finalNorm);
  extendPdf = new RooExtendPdf(Form("extend%sThisLumi",finalPdf->GetName()),Form("extend%sThisLumi",finalPdf->GetName()),*finalPdf,*finalNormThisLum);
  // do secondary models
  if (doSecondaryModels){
    assert(secondaryModelVarsSet);
    // sm higgs as bkg
    RooSpline1D *eaSpline_SM = graphToSpline(Form("fea_%s_cat%d_SM",proc_.c_str(),cat_),eaGraph,MH_SM);
    RooSpline1D *xs_SM = xsSplines_SM[proc_];
    finalNorm_SM = new RooFormulaVar(Form("%s_norm",finalPdf_SM->GetName()),Form("%s_norm",finalPdf_SM->GetName()),"@0*@1*@2",RooArgList(*xs_SM,*brSpline_SM,*eaSpline_SM));
    // second degen higgs
    RooSpline1D *eaSpline_2 = graphToSpline(Form("fea_%s_cat%d_2",proc_.c_str(),cat_),eaGraph,MH_2);
    RooSpline1D *xs_2 = xsSplines_2[proc_];
    finalNorm_2 = new RooFormulaVar(Form("%s_norm",finalPdf_2->GetName()),Form("%s_norm",finalPdf_2->GetName()),"@0*@1*@2",RooArgList(*xs_2,*brSpline_2,*eaSpline_2));
    // natural width
    RooSpline1D *eaSpline_NW = graphToSpline(Form("fea_%s_cat%d_NW",proc_.c_str(),cat_),eaGraph,MH);
    RooSpline1D *xs_NW = xsSplines_NW[proc_];
    finalNorm_NW = new RooFormulaVar(Form("%s_norm",finalPdf_NW->GetName()),Form("%s_norm",finalPdf_NW->GetName()),"@0*@1*@2",RooArgList(*xs_NW,*brSpline_NW,*eaSpline_NW));
  }

}

void FinalModelConstruction::save(RooWorkspace *work){
  work->import(*finalPdf,RecycleConflictNodes());
  work->import(*finalNorm,RecycleConflictNodes());
  work->import(*finalNormThisLum,RecycleConflictNodes());
  work->import(*extendPdf,RecycleConflictNodes());
  work->import(*extendPdfRel,RecycleConflictNodes());
  for (map<int,RooDataSet*>::iterator it=stdDatasets.begin(); it!=stdDatasets.end(); it++){
    work->import(*(it->second));
  }

  // do secondary models
  work->import(*finalPdf_SM,RecycleConflictNodes());
  work->import(*finalPdf_2,RecycleConflictNodes());
  work->import(*finalPdf_NW,RecycleConflictNodes());
  work->import(*finalNorm_SM,RecycleConflictNodes());
  work->import(*finalNorm_2,RecycleConflictNodes());
  work->import(*finalNorm_NW,RecycleConflictNodes());
}
