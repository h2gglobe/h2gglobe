#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooMinimizer.h"
#include "RooExtendPdf.h"

using namespace std;
using namespace RooFit;

int binning_=70;

// get effective sigma from culmalative distribution function
pair<double,double> getEffSigma(RooRealVar *mass, RooAbsPdf *pdf, double wmin=110., double wmax=130., double step=0.002, double epsilon=1.e-4){

  RooAbsReal *cdf = pdf->createCdf(RooArgList(*mass));
  cout << "Computing effSigma...." << endl;
  TStopwatch sw;
  sw.Start();
  double point=wmin;
  vector<pair<double,double> > points;
  
  while (point <= wmax){
    mass->setVal(point);
    if (pdf->getVal() > epsilon){
      points.push_back(pair<double,double>(point,cdf->getVal())); 
    }
    point+=step;
  }
  double low = wmin;
  double high = wmax;
  double width = wmax-wmin;
  for (unsigned int i=0; i<points.size(); i++){
    for (unsigned int j=i; j<points.size(); j++){
      double wy = points[j].second - points[i].second;
      if (TMath::Abs(wy-0.683) < epsilon){
        double wx = points[j].first - points[i].first;
        if (wx < width){
          low = points[i].first;
          high = points[j].first;
          width=wx;
        }
      }
    }
  }
  sw.Stop();
  cout << "effSigma: [" << low << "-" << high << "] = " << width/2. << endl;
  cout << "\tTook: "; sw.Print();
  pair<double,double> result(low,high);
  return result;
}

// get effective sigma from finely binned histogram
pair<double,double> getEffSigBinned(RooRealVar *mass, RooAbsPdf *pdf, double wmin=110., double wmax=130.,int stepsize=1 ){
  
  int nbins = int((wmax-wmin)/0.001/double(stepsize));
  TH1F *h = new TH1F("h","h",nbins,wmin,wmax);
  pdf->fillHistogram(h,RooArgList(*mass));

  double narrowest=1000.;
  double bestInt;
  int lowbin;
  int highbin;

  double oneSigma=1.-TMath::Prob(1,1);
  
  TStopwatch sw;
  sw.Start();
  // get first guess
  cout << "Getting first guess info. stepsize (MeV) = " << stepsize*100 << endl;
  for (int i=0; i<h->GetNbinsX(); i+=(stepsize*100)){
    for (int j=i; j<h->GetNbinsX(); j+=(stepsize*100)){
      double integral = h->Integral(i,j)/h->Integral();
      if (integral<oneSigma) continue;
      double width = h->GetBinCenter(j)-h->GetBinCenter(i);
      if (width<narrowest){
        narrowest=width;
        bestInt=integral;
        lowbin=i;
        highbin=j;
        i++;
      }
    }
  }
  cout << "Took: "; sw.Print(); 
  // narrow down result
  int thisStepSize=32;
  cout << "Narrowing....." << endl;
  while (thisStepSize>stepsize) {
    cout << "\tstepsize (MeV) = " << thisStepSize << endl;
    for (int i=(lowbin-10*thisStepSize); i<(lowbin+10*thisStepSize); i+=thisStepSize){
      for (int j=(highbin-10*thisStepSize); j<(highbin+10*thisStepSize); j+=thisStepSize){
        double integral = h->Integral(i,j)/h->Integral();
        if (integral<oneSigma) continue;
        double width = h->GetBinCenter(j)-h->GetBinCenter(i);
        if (width<narrowest){
          narrowest=width;
          bestInt=integral;
          lowbin=i;
          highbin=j;
          i++;
        }
      }
    }
    thisStepSize/=2;
  }

  sw.Stop();
  cout << narrowest/2. << " " << bestInt << " [" << h->GetBinCenter(lowbin) << "," << h->GetBinCenter(highbin) << "]" << endl;
  cout << "Took:"; sw.Print();
  pair<double,double> result(h->GetBinCenter(lowbin),h->GetBinCenter(highbin));
  delete h;
  return result;
}

// get FWHHM
vector<double> getFWHM(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, double wmin=110., double wmax=130., double step=0.0004) {
 
  cout << "Computing FWHM...." << endl;
  double nbins = (wmax-wmin)/step;
  TH1F *h = new TH1F("h","h",int(floor(nbins+0.5)),wmin,wmax);
  if (data){
    pdf->fillHistogram(h,RooArgList(*mass),data->sumEntries());
  }
  else {
    pdf->fillHistogram(h,RooArgList(*mass));
  }
  
  double hm = h->GetMaximum()*0.5;
  double low = h->GetBinCenter(h->FindFirstBinAbove(hm));
  double high = h->GetBinCenter(h->FindLastBinAbove(hm));

  cout << "FWHM: [" << low << "-" << high << "] Max = " << hm << endl;
  vector<double> result;
  result.push_back(low);
  result.push_back(high);
  result.push_back(hm);
  result.push_back(h->GetBinWidth(1));
  
  delete h;
  return result;
}

void performClosure(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string closurename, double wmin=110., double wmax=130., double slow=110., double shigh=130., double step=0.002) {
  
  // plot to perform closure test
  cout << "Performing closure test..." << endl; 
  double nbins = (wmax-wmin)/step;
  TH1F *h = new TH1F("h","h",int(floor(nbins+0.5)),wmin,wmax);
  if (data){
    pdf->fillHistogram(h,RooArgList(*mass),data->sumEntries());
    h->Scale(2*h->GetNbinsX()/double(binning_));
  }
  else {
    pdf->fillHistogram(h,RooArgList(*mass));
  }
  int binLow = h->FindBin(slow);
  int binHigh = h->FindBin(shigh)-1;
  TH1F *copy = new TH1F("copy","c",binHigh-binLow,h->GetBinLowEdge(binLow),h->GetBinLowEdge(binHigh+1));
  for (int b=0; b<copy->GetNbinsX(); b++) copy->SetBinContent(b+1,h->GetBinContent(b+1+binLow));
  double areaCov = 100*h->Integral(binLow,binHigh)/h->Integral();
 
  // style
  h->SetLineColor(kBlue);
  h->SetLineWidth(3);
  h->SetLineStyle(7);
  copy->SetLineWidth(3);
  copy->SetFillColor(kGray);
  
  TCanvas *c = new TCanvas();
  if (data){
    RooPlot *plot = mass->frame(Bins(binning_),Range("higgsRange"));
    plot->addTH1(h,"hist");
    plot->addTH1(copy,"same f");
    if (data) data->plotOn(plot);
    pdf->plotOn(plot,Normalization(h->Integral(),RooAbsReal::NumEvent),NormRange("higgsRange"),Range("higgsRange"),LineWidth(1),LineColor(kRed),LineStyle(kDashed));
    plot->Draw();
    c->Print(closurename.c_str());
  }
  else {
    RooPlot *plot = mass->frame(Bins(binning_),Range("higgsRange"));
    h->Scale(plot->getFitRangeBinW()/h->GetBinWidth(1));
    copy->Scale(plot->getFitRangeBinW()/h->GetBinWidth(1));
    pdf->plotOn(plot,LineColor(kRed),LineWidth(3));
    plot->Draw();
    h->Draw("hist same");
    copy->Draw("same f");
    c->Print(closurename.c_str());
  }
  cout << "IntH: [" << h->GetBinLowEdge(binLow) << "-" << h->GetBinLowEdge(binHigh+1) << "] Area = " << areaCov << endl;
  delete c;
  delete copy;
  delete h;
}

void Plot(RooRealVar *mass, RooDataSet *data, RooAbsPdf *pdf, pair<double,double> sigRange, vector<double> fwhmRange, string title, string savename){

  double semin=sigRange.first;
  double semax=sigRange.second;
  double fwmin=fwhmRange[0];
  double fwmax=fwhmRange[1];
  double halfmax=fwhmRange[2];
  double binwidth=fwhmRange[3];

  RooPlot *plot = mass->frame(Bins(binning_),Range("higgsRange"));
  if (data) data->plotOn(plot,Invisible());
  pdf->plotOn(plot,NormRange("higgsRange"),Range(semin,semax),FillColor(19),DrawOption("F"),LineWidth(2),FillStyle(1001),VLines(),LineColor(15));
  TObject *seffLeg = plot->getObject(int(plot->numItems()-1));
  pdf->plotOn(plot,NormRange("higgsRange"),Range(semin,semax),LineColor(15),LineWidth(2),FillStyle(1001),VLines());
  pdf->plotOn(plot,NormRange("higgsRange"),Range("higgsRange"),LineColor(kBlue),LineWidth(2),FillStyle(0));
  TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
  if (data) data->plotOn(plot,MarkerStyle(kOpenSquare));
  TObject *dataLeg = plot->getObject(int(plot->numItems()-1));
  TLegend *leg = new TLegend(0.15,0.89,0.5,0.55);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.03);
  if (data) leg->AddEntry(dataLeg,"Simulation","lep");
  leg->AddEntry(pdfLeg,"Parametric model","l");
  leg->AddEntry(seffLeg,Form("#sigma_{eff} = %1.2f GeV",0.5*(semax-semin)),"fl");

  plot->GetXaxis()->SetNdivisions(509);
  halfmax*=(plot->getFitRangeBinW()/binwidth);
  TArrow *fwhmArrow = new TArrow(fwmin,halfmax,fwmax,halfmax,0.02,"<>");
  fwhmArrow->SetLineWidth(2.);
  TPaveText *fwhmText = new TPaveText(0.15,0.45,0.45,0.58,"brNDC");
  fwhmText->SetFillColor(0);
  fwhmText->SetLineColor(kWhite);
  fwhmText->SetTextSize(0.03);
  fwhmText->AddText(Form("FWHM = %1.2f GeV",(fwmax-fwmin)));

  TLatex lat1(0.65,0.85,"#splitline{CMS Preliminary}{Simulation}");
  lat1.SetNDC(1);
  lat1.SetTextSize(0.03);
  TLatex lat2(0.65,0.75,title.c_str());
  lat2.SetNDC(1);
  lat2.SetTextSize(0.025);

  TCanvas *canv = new TCanvas("c","c",600,600);
  plot->SetTitle("");
  plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  plot->Draw();
  leg->Draw("same");
  fwhmArrow->Draw("same <>");
  fwhmText->Draw("same");
  lat1.Draw("same");
  lat2.Draw("same");
  canv->Print(Form("%s.pdf",savename.c_str()));
  canv->Print(Form("%s.png",savename.c_str()));
  string path = savename.substr(0,savename.find('/'));
  canv->Print(Form("%s/animation.gif+100",path.c_str()));
  delete canv;

}

string runQuickRejig(string oldfilename, int ncats){
  
  string newfilename="TestWS.root";
  TFile *oldFile = TFile::Open(oldfilename.c_str());
  TFile *newFile = new TFile(newfilename.c_str(),"RECREATE");

  RooWorkspace *oldWS = (RooWorkspace*)oldFile->Get("cms_hgg_workspace");
  RooWorkspace *newWS = new RooWorkspace("wsig_8TeV");
  
  RooRealVar *mass = (RooRealVar*)oldWS->var("CMS_hgg_mass");

  RooDataSet *tot;

  RooRealVar *normT1 = new RooRealVar("nT1","nT1",500,0,550);
  RooRealVar *meanT1 = new RooRealVar("mT1","mT1",125,122,128);
  RooRealVar *sigmaT1 = new RooRealVar("sT1","sT1",3.1,1.,10.);
  RooGaussian *gausT1 = new RooGaussian("gT1","gT1",*mass,*meanT1,*sigmaT1);
  RooRealVar *normT2 = new RooRealVar("nT2","nT2",2,0,500);
  RooRealVar *meanT2 = new RooRealVar("mT2","mT2",125,122,128);
  RooRealVar *sigmaT2 = new RooRealVar("sT2","sT2",1.7,0.,10.);
  RooGaussian *gausT2 = new RooGaussian("gT2","gT2",*mass,*meanT2,*sigmaT2);
  RooRealVar *normT3 = new RooRealVar("nT3","nT3",2,0,500);
  RooRealVar *meanT3 = new RooRealVar("mT3","mT3",125,122,128);
  RooRealVar *sigmaT3 = new RooRealVar("sT3","sT3",1.7,0.,10.);
  RooGaussian *gausT3 = new RooGaussian("gT3","gT3",*mass,*meanT3,*sigmaT3);
  RooAddPdf *gausT = new RooAddPdf("gT","gT",RooArgList(*gausT1,*gausT2,*gausT3),RooArgList(*normT1,*normT2,*normT3));

  for (int cat=0; cat<ncats; cat++){
    RooDataSet *thisData = (RooDataSet*)oldWS->data(Form("sig_ggh_mass_m125_cat%d",cat));
    newWS->import(*thisData);
    RooRealVar *norm1 = new RooRealVar(Form("n1%d",cat),Form("n1%d",cat),500,0,550);
    RooRealVar *mean1 = new RooRealVar(Form("m1%d",cat),Form("m1%d",cat),125,122,128);
    RooRealVar *sigma1 = new RooRealVar(Form("s1%d",cat),Form("s1%d",cat),3.1,1.,10.);
    RooGaussian *gaus1 = new RooGaussian(Form("g1%d",cat),Form("g1%d",cat),*mass,*mean1,*sigma1);
    RooRealVar *norm2 = new RooRealVar(Form("n2%d",cat),Form("n2%d",cat),2,0,500);
    RooRealVar *mean2 = new RooRealVar(Form("m2%d",cat),Form("m2%d",cat),125,122,128);
    RooRealVar *sigma2 = new RooRealVar(Form("s2%d",cat),Form("s2%d",cat),1.7,0.,10.);
    RooGaussian *gaus2 = new RooGaussian(Form("g2%d",cat),Form("g2%d",cat),*mass,*mean2,*sigma2);
    RooRealVar *norm3 = new RooRealVar(Form("n3%d",cat),Form("n3%d",cat),2,0,500);
    RooRealVar *mean3 = new RooRealVar(Form("m3%d",cat),Form("m3%d",cat),125,122,128);
    RooRealVar *sigma3 = new RooRealVar(Form("s3%d",cat),Form("s3%d",cat),1.7,0.,10.);
    RooGaussian *gaus3 = new RooGaussian(Form("g3%d",cat),Form("g3%d",cat),*mass,*mean3,*sigma3);
    RooAddPdf *gaus = new RooAddPdf(Form("g%d",cat),"g",RooArgList(*gaus1,*gaus2,*gaus3),RooArgList(*norm1,*norm2,*norm3));
    gaus->fitTo(*thisData,SumW2Error(kTRUE));
    newWS->import(*gaus);
    if (cat==0) {
      tot = thisData;
      tot->SetName("sig_ggh_m125");
    }
    else tot->append(*thisData);
  }
  newWS->import(*tot);
  gausT->fitTo(*tot,SumW2Error(kTRUE));
  newWS->import(*gausT);
  newWS->Write();
  oldFile->Close();
  newFile->Close();
  delete newFile;
  delete newWS;

  return newfilename;

}

void printInfo(map<string,RooDataSet*> data, map<string,RooAddPdf*> pdfs){
  
  for (map<string,RooDataSet*>::iterator dat=data.begin(); dat!=data.end(); dat++){
    if (!dat->second) {
      cout << "Dataset for " << dat->first << " not found" << endl;
      exit(1);
    }
    cout << dat->first << " : ";
    dat->second->Print();
  }
  for (map<string,RooAddPdf*>::iterator pdf=pdfs.begin(); pdf!=pdfs.end(); pdf++){
    if (!pdf->second) {
      cout << "Pdf for " << pdf->first << " not found" << endl;
      exit(1);
    }
    cout << pdf->first << " : ";
    pdf->second->Print();
  }

}

map<string,RooDataSet*> getGlobeData(RooWorkspace *work, int ncats, int m_hyp){
  
  map<string,RooDataSet*> result;
  for (int cat=0; cat<ncats; cat++){
    result.insert(pair<string,RooDataSet*>(Form("cat%d",cat),(RooDataSet*)work->data(Form("sig_mass_m%3d_cat%d",m_hyp,cat))));
  }
  result.insert(pair<string,RooDataSet*>("all",(RooDataSet*)work->data(Form("sig_mass_m%3d_AllCats",m_hyp))));

  return result;
}

map<string,RooAddPdf*> getGlobePdfs(RooWorkspace *work, int ncats){

  map<string,RooAddPdf*> result;
  for (int cat=0; cat<ncats; cat++){
    result.insert(pair<string,RooAddPdf*>(Form("cat%d",cat),(RooAddPdf*)work->pdf(Form("sigpdfrelcat%d_allProcs",cat))));
  }
  result.insert(pair<string,RooAddPdf*>("all",(RooAddPdf*)work->pdf("sigpdfrelAllCats_allProcs")));

  return result;
}

map<string,RooAddPdf*> getMITPdfs(RooWorkspace *work, int ncats, bool is2011){

  map<string,RooAddPdf*> result;

  RooArgList pdfs;
  RooArgList norms;

  string sqrts;
  string catstring;
  if (is2011){
    sqrts="7TeV";
    catstring="cat";
  }
  else {
    sqrts="8TeV";
    catstring="mvacat";
  }

  for (int cat=0; cat<ncats; cat++){
    RooAddPdf *ggh = (RooAddPdf*)work->pdf(Form("hggpdfsmrel_%s%d_%s_ggh",catstring.c_str(),cat,sqrts.c_str())); 
    RooAddPdf *vbf = (RooAddPdf*)work->pdf(Form("hggpdfsmrel_%s%d_%s_vbf",catstring.c_str(),cat,sqrts.c_str())); 
    RooAddPdf *wzh = (RooAddPdf*)work->pdf(Form("hggpdfsmrel_%s%d_%s_wzh",catstring.c_str(),cat,sqrts.c_str())); 
    RooAddPdf *tth = (RooAddPdf*)work->pdf(Form("hggpdfsmrel_%s%d_%s_tth",catstring.c_str(),cat,sqrts.c_str()));
    RooFormulaVar *ggh_norm = (RooFormulaVar*)work->function(Form("hggpdfsmrel_%s%d_%s_ggh_norm",catstring.c_str(),cat,sqrts.c_str()));
    RooFormulaVar *vbf_norm = (RooFormulaVar*)work->function(Form("hggpdfsmrel_%s%d_%s_vbf_norm",catstring.c_str(),cat,sqrts.c_str()));
    RooFormulaVar *wzh_norm = (RooFormulaVar*)work->function(Form("hggpdfsmrel_%s%d_%s_wzh_norm",catstring.c_str(),cat,sqrts.c_str()));
    RooFormulaVar *tth_norm = (RooFormulaVar*)work->function(Form("hggpdfsmrel_%s%d_%s_tth_norm",catstring.c_str(),cat,sqrts.c_str()));
    pdfs.add(*ggh); pdfs.add(*vbf); pdfs.add(*wzh); pdfs.add(*tth);
    norms.add(*ggh_norm); norms.add(*vbf_norm); norms.add(*wzh_norm); norms.add(*tth_norm);
    result.insert(pair<string,RooAddPdf*>(Form("cat%d",cat),new RooAddPdf(Form("sig_cat%d",cat),"sig",RooArgList(*ggh,*vbf,*wzh,*tth),RooArgList(*ggh_norm,*vbf_norm,*wzh_norm,*tth_norm))));
  }
  result.insert(pair<string,RooAddPdf*>("all",new RooAddPdf("sig_all","sig",pdfs,norms)));
  return result;

}

pair<double,double> bkgEvPerGeV(RooWorkspace *work, int m_hyp, int cat, pair<double,double> &bkgTotal){
  
  RooRealVar *mass = (RooRealVar*)work->var("CMS_hgg_mass");
  mass->setRange(100,180);
  RooAbsPdf *pdf = (RooAbsPdf*)work->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
  RooAbsData *data = (RooDataSet*)work->data(Form("data_mass_cat%d",cat));
  RooPlot *tempFrame = mass->frame();
  data->plotOn(tempFrame,Binning(80));
  pdf->plotOn(tempFrame);
  RooCurve *curve = (RooCurve*)tempFrame->getObject(tempFrame->numItems()-1);
  double nombkg = curve->Eval(double(m_hyp));
 
  RooRealVar *nlim = new RooRealVar(Form("nlim%d",cat),"",0.,0.,1.e5);
  //double lowedge = tempFrame->GetXaxis()->GetBinLowEdge(FindBin(double(m_hyp)));
  //double upedge  = tempFrame->GetXaxis()->GetBinUpEdge(FindBin(double(m_hyp)));
  //double center  = tempFrame->GetXaxis()->GetBinUpCenter(FindBin(double(m_hyp)));

  nlim->setVal(nombkg);
  mass->setRange("errRange",m_hyp-0.5,m_hyp+0.5);
  RooAbsPdf *epdf = 0;
  epdf = new RooExtendPdf("epdf","",*pdf,*nlim,"errRange");
		
  RooAbsReal *nll = epdf->createNLL(*data,Extended(),NumCPU(4));
  RooMinimizer minim(*nll);
  minim.setStrategy(0);
  minim.setPrintLevel(-1);
  minim.migrad();
  minim.minos(*nlim);
  
  double error = (nlim->getErrorLo(),nlim->getErrorHi())/2.;
  data->Print();
  bkgTotal.first += nombkg;
  bkgTotal.second += error*error;
  return pair<double,double>(nombkg,error); 
}

vector<double> sigEvents(RooWorkspace *work, int m_hyp, int cat, string altSigFileName, vector<double> &sigTotal, string spinProc=""){

  RooWorkspace *tempWork;
  if (altSigFileName!=""){
    TFile *temp = TFile::Open(altSigFileName.c_str());
    tempWork = (RooWorkspace*)temp->Get("cms_hgg_workspace");
  }
  else {
    tempWork = work;
  }
  vector<double> result;
  
  if (spinProc==""){
    
    RooDataSet *ggh = (RooDataSet*)tempWork->data(Form("sig_ggh_mass_m%d_cat%d",m_hyp,cat));
    RooDataSet *vbf = (RooDataSet*)tempWork->data(Form("sig_vbf_mass_m%d_cat%d",m_hyp,cat));
    RooDataSet *wzh = (RooDataSet*)tempWork->data(Form("sig_wzh_mass_m%d_cat%d",m_hyp,cat));
    RooDataSet *tth = (RooDataSet*)tempWork->data(Form("sig_tth_mass_m%d_cat%d",m_hyp,cat));
    
    double total = ggh->sumEntries()+vbf->sumEntries()+wzh->sumEntries()+tth->sumEntries();
    result.push_back(total);
    result.push_back(ggh->sumEntries());
    result.push_back(vbf->sumEntries());
    result.push_back(wzh->sumEntries());
    result.push_back(tth->sumEntries());

    sigTotal[0] += total;
    sigTotal[1] += ggh->sumEntries();
    sigTotal[2] += vbf->sumEntries();
    sigTotal[3] += wzh->sumEntries();
    sigTotal[4] += tth->sumEntries();
  }
  else if (spinProc=="gg_grav"){
    
    RooDataHist *ggh = (RooDataHist*)tempWork->data(Form("roohist_sig_gg_grav_mass_m%d_cat%d",m_hyp,cat));

    double total = ggh->sumEntries();
    result.push_back(total);
    result.push_back(ggh->sumEntries());
    result.push_back(0.);
    result.push_back(0.);
    result.push_back(0.);

    sigTotal[0] += total;
    sigTotal[1] += ggh->sumEntries();
  }
  else if (spinProc=="qq_grav"){
    RooDataHist *vbf = (RooDataHist*)tempWork->data(Form("roohist_sig_qq_grav_mass_m%d_cat%d",m_hyp,cat));
    
    double total = vbf->sumEntries();
    result.push_back(total);
    result.push_back(0.);
    result.push_back(vbf->sumEntries());
    result.push_back(0.);
    result.push_back(0.);

    sigTotal[0] += total;
    sigTotal[1] += vbf->sumEntries();
  }
  else {
    cout << "WARNING -- spinProc " << spinProc << " is no recognised." << endl;
  }

  return result;
}

pair<double,double> datEvents(RooWorkspace *work, int m_hyp, int cat, pair<double,double> &runningTotal){
  
  vector<double> result;
  RooDataSet *data = (RooDataSet*)work->data(Form("data_mass_cat%d",cat));
  double evs = data->numEntries();
  double evsPerGev;
  evsPerGev = data->sumEntries(Form("CMS_hgg_mass>=%4.1f && CMS_hgg_mass<%4.1f",double(m_hyp)-0.5,double(m_hyp)+0.5));
  runningTotal.first += evs;
  runningTotal.second += evsPerGev;
  return pair<double,double>(evs,evsPerGev);
}

void makeParametricSignalModelPlots(string sigFitFileName, string pathName, int ncats=9, bool is2011=false, int m_hyp=120, string bkgdatFileName="0", bool isMassFac = true, bool blind=true, bool doTable=true, string altSigFileName="", bool spin=false, string spinProc="", bool doCrossCheck=false, bool doMIT=false, bool rejig=false){

  gROOT->SetBatch();
  gStyle->SetTextFont(42);

  string newFileName;
  if (rejig) newFileName = runQuickRejig(sigFitFileName,ncats);
  else newFileName = sigFitFileName;

  TFile *hggFile = TFile::Open(newFileName.c_str());
  
  string sqrts;
  if (is2011) sqrts="7TeV";
  else sqrts="8TeV";

  RooWorkspace *hggWS;
  if (is2011) hggWS = (RooWorkspace*)hggFile->Get(Form("wsig_%s",sqrts.c_str()));
  else        hggWS = (RooWorkspace*)hggFile->Get(Form("wsig_%s",sqrts.c_str()));
 
  if (!hggWS) {
    cerr << "Workspace is null" << endl;
    exit(1);
  }

  RooRealVar *mass= (RooRealVar*)hggWS->var("CMS_hgg_mass");
  
  RooRealVar *mh = (RooRealVar*)hggWS->var("MH");
  mh->setVal(m_hyp);
  mass->setRange("higgsRange",m_hyp-20.,m_hyp+15.);

  map<string,string> labels;
  if (is2011) {
    if (isMassFac){
      labels.insert(pair<string,string>("cat0","BDT >= 0.89"));
      labels.insert(pair<string,string>("cat1","0.72 <= BDT <= 0.89"));
      labels.insert(pair<string,string>("cat2","0.55 <= BDT <= 0.72"));
      labels.insert(pair<string,string>("cat3","0.05 <= BDT <= 0.55"));
    }
    else {
      labels.insert(pair<string,string>("cat0","EBEB, min(R9) > 0.94"));
      labels.insert(pair<string,string>("cat1","EBEB, min(R9) < 0.94"));
      labels.insert(pair<string,string>("cat2","!EBEB, min(R9) > 0.94"));
      labels.insert(pair<string,string>("cat3","!EBEB, min(R9) < 0.94"));
    }
    labels.insert(pair<string,string>("cat4","BDT >= 0.05 VBF Tag"));
    labels.insert(pair<string,string>("all","All Categories Combined"));
  }
  else {
    if (isMassFac){
      labels.insert(pair<string,string>("cat0","BDT_{#gamma#gamma} >= 0.91"));
      labels.insert(pair<string,string>("cat1","0.79 <= BDT_{#gamma#gamma} <= 0.91"));
      labels.insert(pair<string,string>("cat2","0.49 <= BDT_{#gamma#gamma} <= 0.79"));
      labels.insert(pair<string,string>("cat3","-0.05 <= BDT_{#gamma#gamma} <= 0.49"));
    }
    else {
      labels.insert(pair<string,string>("cat0","EBEB, min(R9) > 0.94"));
      labels.insert(pair<string,string>("cat1","EBEB, min(R9) < 0.94"));
      labels.insert(pair<string,string>("cat2","!EBEB, min(R9) > 0.94"));
      labels.insert(pair<string,string>("cat3","!EBEB, min(R9) < 0.94"));
    }
    labels.insert(pair<string,string>("cat4","BDT_{jj} >= 0.985 Dijet Tag"));
    labels.insert(pair<string,string>("cat5","BDT_{jj} >= 0.93 Dijet Tag")); 
    labels.insert(pair<string,string>("cat6","BDT_{#gamma#gamma} >= -0.05 Muon Tag")); 
    labels.insert(pair<string,string>("cat7","BDT_{#gamma#gamma} >= -0.05 Electron Tag")); 
    labels.insert(pair<string,string>("cat8","BDT_{#gamma#gamma} >= -0.05 MET Tag")); 
    labels.insert(pair<string,string>("all","All Categories Combined"));
  }
  if (spin){
    labels.clear();
    labels.insert(pair<string,string>("cat0","#splitline{|#eta|_{max} < 1.44, R_{9min} > 0.94}{|cos(#theta*)| < 0.2}"));
    labels.insert(pair<string,string>("cat1","#splitline{|#eta|_{max} < 1.44, R_{9min} > 0.94}{0.2 < |cos(#theta*)| < 0.375}"));
    labels.insert(pair<string,string>("cat2","#splitline{|#eta|_{max} < 1.44, R_{9min} > 0.94}{0.375 < |cos(#theta*)| < 0.55}"));
    labels.insert(pair<string,string>("cat3","#splitline{|#eta|_{max} < 1.44, R_{9min} > 0.94}{0.55 < |cos(#theta*)| < 0.75}"));
    labels.insert(pair<string,string>("cat4","#splitline{|#eta|_{max} < 1.44, R_{9min} > 0.94}{0.75 < |cos(#theta*)| < 0.1}"));
    labels.insert(pair<string,string>("cat5","#splitline{|#eta|_{max} < 1.44, R_{9min} < 0.94}{|cos(#theta*)| < 0.2}"));
    labels.insert(pair<string,string>("cat6","#splitline{|#eta|_{max} < 1.44, R_{9min} < 0.94}{0.2 < |cos(#theta*)| < 0.375}"));
    labels.insert(pair<string,string>("cat7","#splitline{|#eta|_{max} < 1.44, R_{9min} < 0.94}{0.375 < |cos(#theta*)| < 0.55}"));
    labels.insert(pair<string,string>("cat8","#splitline{|#eta|_{max} < 1.44, R_{9min} < 0.94}{0.55 < |cos(#theta*)| < 0.75}"));
    labels.insert(pair<string,string>("cat9","#splitline{|#eta|_{max} < 1.44, R_{9min} < 0.94}{0.75 < |cos(#theta*)| < 0.1}"));
    labels.insert(pair<string,string>("cat10","#splitline{|#eta|_{max} > 1.44, R_{9min} > 0.94}{|cos(#theta*)| < 0.2}"));
    labels.insert(pair<string,string>("cat11","#splitline{|#eta|_{max} > 1.44, R_{9min} > 0.94}{0.2 < |cos(#theta*)| < 0.375}"));
    labels.insert(pair<string,string>("cat12","#splitline{|#eta|_{max} > 1.44, R_{9min} > 0.94}{0.375 < |cos(#theta*)| < 0.55}"));
    labels.insert(pair<string,string>("cat13","#splitline{|#eta|_{max} > 1.44, R_{9min} > 0.94}{0.55 < |cos(#theta*)| < 0.75}"));
    labels.insert(pair<string,string>("cat14","#splitline{|#eta|_{max} > 1.44, R_{9min} > 0.94}{0.75 < |cos(#theta*)| < 0.1}"));
    labels.insert(pair<string,string>("cat15","#splitline{|#eta|_{max} > 1.44, R_{9min} < 0.94}{|cos(#theta*)| < 0.2}"));
    labels.insert(pair<string,string>("cat16","#splitline{|#eta|_{max} > 1.44, R_{9min} < 0.94}{0.2 < |cos(#theta*)| < 0.375}"));
    labels.insert(pair<string,string>("cat17","#splitline{|#eta|_{max} > 1.44, R_{9min} < 0.94}{0.375 < |cos(#theta*)| < 0.55}"));
    labels.insert(pair<string,string>("cat18","#splitline{|#eta|_{max} > 1.44, R_{9min} < 0.94}{0.55 < |cos(#theta*)| < 0.75}"));
    labels.insert(pair<string,string>("cat19","#splitline{|#eta|_{max} > 1.44, R_{9min} < 0.94}{0.75 < |cos(#theta*)| < 0.1}"));
    labels.insert(pair<string,string>("all","All Categories Combined"));
  }
  /*
  else {
    labels.insert(pair<string,string>("cat0","BDT >= 0.88"));
    labels.insert(pair<string,string>("cat1","0.71 <= BDT <= 0.88"));
    labels.insert(pair<string,string>("cat2","0.50 <= BDT <= 0.71"));
    labels.insert(pair<string,string>("cat3","-0.05 <= BDT <= 0.50"));
    labels.insert(pair<string,string>("cat4","BDT >= -0.05 Tight VBF"));
    labels.insert(pair<string,string>("cat5","BDT >= -0.05 Loose VBF")); 
    labels.insert(pair<string,string>("all","All Categories Combined"));
  }
  */

  map<string,RooDataSet*> dataSets;
  map<string,RooAddPdf*> pdfs;

  if (doMIT){
    dataSets = getGlobeData(hggWS,ncats,m_hyp);
    pdfs = getMITPdfs(hggWS,ncats,is2011);
  }
  else {
    dataSets = getGlobeData(hggWS,ncats,m_hyp);
    pdfs = getGlobePdfs(hggWS,ncats); 
  }
  
  printInfo(dataSets,pdfs);

  map<string,double> sigEffs;
  map<string,double> fwhms;
  
  system(Form("mkdir -p %s",pathName.c_str()));
  system(Form("rm %s/animation.gif",pathName.c_str()));
  for (map<string,RooDataSet*>::iterator dataIt=dataSets.begin(); dataIt!=dataSets.end(); dataIt++){
    pair<double,double> thisSigRange = getEffSigma(mass,pdfs[dataIt->first],m_hyp-10.,m_hyp+10.);
    //pair<double,double> thisSigRange = getEffSigBinned(mass,pdf[dataIt->first],m_hyp-10.,m_hyp+10);
    vector<double> thisFWHMRange = getFWHM(mass,pdfs[dataIt->first],dataIt->second,m_hyp-10.,m_hyp+10.);
    sigEffs.insert(pair<string,double>(dataIt->first,(thisSigRange.second-thisSigRange.first)/2.));
    fwhms.insert(pair<string,double>(dataIt->first,thisFWHMRange[1]-thisFWHMRange[0]));
    if (doCrossCheck) performClosure(mass,pdfs[dataIt->first],dataIt->second,Form("%s/closure_%s.pdf",pathName.c_str(),dataIt->first.c_str()),m_hyp-10.,m_hyp+10.,thisSigRange.first,thisSigRange.second);
    Plot(mass,dataIt->second,pdfs[dataIt->first],thisSigRange,thisFWHMRange,labels[dataIt->first],Form("%s/%s",pathName.c_str(),dataIt->first.c_str()));
  }
  
  map<string,pair<double,double> > bkgVals;
  map<string,vector<double> > sigVals;
  map<string,pair<double,double> > datVals;

  // make PAS table
  if (doTable) {
    TFile *bkgFile;
    if (bkgdatFileName!="0"){
      bkgFile = TFile::Open(bkgdatFileName.c_str());
    }
    else {
      bkgFile = hggFile; 
    }
    RooWorkspace *bkgWS = (RooWorkspace*)bkgFile->Get("cms_hgg_workspace");
    // keep track of sums
    pair<double,double> bkgTotal(0.,0.);
    pair<double,double> datTotal(0.,0.);
    vector<double> sigTotal;
    for (int i=0; i<5; i++) sigTotal.push_back(0.);
    for (int cat=0; cat<ncats; cat++){
      bkgVals.insert(pair<string,pair<double,double> >(Form("cat%d",cat),bkgEvPerGeV(bkgWS,m_hyp,cat,bkgTotal)));
      sigVals.insert(pair<string,vector<double> >(Form("cat%d",cat),sigEvents(bkgWS,m_hyp,cat,altSigFileName,sigTotal,spinProc)));
      datVals.insert(pair<string,pair<double,double> >(Form("cat%d",cat),datEvents(bkgWS,m_hyp,cat,datTotal)));
    }
    bkgTotal.second = sqrt(bkgTotal.second);
    bkgVals.insert(pair<string,pair<double,double> >("all",bkgTotal));
    sigVals.insert(pair<string,vector<double> > ("all",sigTotal));
    datVals.insert(pair<string,pair<double,double> >("all",datTotal));
    bkgFile->Close();
    
    FILE *file = fopen(Form("%s/table.tex",pathName.c_str()),"w");
    FILE *nfile = fopen(Form("%s/table.txt",pathName.c_str()),"w");
    fprintf(nfile,"----------------------------------------------------------------------------------------------\n");
    fprintf(nfile,"Cat    SigY    ggh    vbf    wzh    tth   sEff  FWHM  FWHM/2.35  BkgEv/GeV    Data  DataEv/GeV\n");
    fprintf(nfile,"----------------------------------------------------------------------------------------------\n");
    for (int cat=0; cat<=ncats; cat++){
      string thisCatName;
      if (cat==ncats) thisCatName = "all";
      else thisCatName = Form("cat%d",cat);
      pair<double,double> bkg = bkgVals[thisCatName];
      vector<double> sigs = sigVals[thisCatName];
      pair<double,double> dat = datVals[thisCatName];
      // print to file
      fprintf(nfile,"%5s  ",thisCatName.c_str());
      fprintf(nfile,"%5.1f  ",sigs[0]);
      fprintf(nfile,"%4.1f%%  ",100.*sigs[1]/sigs[0]);
      fprintf(nfile,"%4.1f%%  ",100.*sigs[2]/sigs[0]);
      fprintf(nfile,"%4.1f%%  ",100.*sigs[3]/sigs[0]);
      fprintf(nfile,"%4.1f%%  ",100.*sigs[4]/sigs[0]);
      fprintf(nfile,"%4.2f  ",sigEffs[thisCatName]);
      fprintf(nfile,"%4.2f  ",fwhms[thisCatName]);
      fprintf(nfile,"%4.2f   ",fwhms[thisCatName]/2.35);
      fprintf(nfile,"%6.1f +/- %3.1f  ",bkg.first,bkg.second);
      fprintf(nfile,"%6.0f    ",dat.first);
      if (blind) fprintf(nfile,"%7s","----");
      else fprintf(nfile,"%7.1f  ",dat.second);
      fprintf(nfile,"\n");
      // print to file
      fprintf(file,"&  cat%d  ",cat);
      fprintf(file,"&  %5.1f  ",sigs[0]);
      fprintf(file,"&  %4.1f\\%%  ",100.*sigs[1]/sigs[0]);
      fprintf(file,"&  %4.1f\\%%  ",100.*sigs[2]/sigs[0]);
      fprintf(file,"&  %4.1f\\%%  ",100.*sigs[3]/sigs[0]);
      fprintf(file,"&  %4.1f\\%%  ",100.*sigs[4]/sigs[0]);
      fprintf(file,"&  %4.2f  ",sigEffs[thisCatName]);
      fprintf(file,"&  %4.2f  ",fwhms[thisCatName]);
      fprintf(file,"&  %4.2f  ",fwhms[thisCatName]/2.35);
      fprintf(file,"&  %6.1f & $\\pm$ %3.1f \\tabularnewline ",bkg.first,bkg.second);
      fprintf(file,"&  %7.1f  ",dat.first);
      if (blind) fprintf(file,"& %7s","----");
      else fprintf(file,"&  %7.1f  ",dat.second);
      fprintf(file,"\n");
    }
    fclose(nfile);
    fclose(file);
    system(Form("cat %s/table.txt",pathName.c_str()));
    cout << "-->" << endl;
    cout << Form("--> LaTeX version of this table has been written to %s/table.tex",pathName.c_str()) << endl;
  }

  hggFile->Close();
}
