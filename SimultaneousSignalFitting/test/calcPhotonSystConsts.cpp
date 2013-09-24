#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

string filename_;
string outfilename_;
string photonCatStr_;
vector<string> photonCats_;
string procStr_;
vector<string> procs_;
int mh_;
int nCats_;

void OptionParser(int argc, char *argv[]){
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                					"Show help")
    ("infilename,i", po::value<string>(&filename_),                                           					"Input file name")
    ("outfilename,o", po::value<string>(&outfilename_)->default_value("dat/photonCatSyst.dat"), 				"Output file name")
    ("mh,m", po::value<int>(&mh_)->default_value(125),                                  								"Mass point")
    ("nCats,n", po::value<int>(&nCats_)->default_value(9),                                    					"Number of total categories")
		("phoCats,c",po::value<string>(&photonCatStr_)->default_value("EBlowR9,EBhighR9,EElowR9,EEhighR9"),	"Photon cats (comma sep)")
		("procs,p",po::value<string>(&procStr_)->default_value("ggh,vbf,wh,zh,tth"),												"Processes (comma sep)")
  ;                                                                                             		
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")){ cout << desc << endl; exit(1);}
	split(photonCats_,photonCatStr_,boost::is_any_of(","));
	split(procs_,procStr_,boost::is_any_of(","));
}

//effsigma function from Chris
Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
 
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  //Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
//   if(total < 100.) {
//     cout << "effsigma: Too few entries " << total << endl;
//     return 0.;
//   }
  Int_t ierr=0;
  Int_t ismin=999;
 
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }  
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
 
  return widmin;
 
}

void plotVariation(TH1F *nom, TH1F *up, TH1F *down, string phoCat, string name){

	TCanvas *canv = new TCanvas();
	canv->SetWindowSize(750,750);
	canv->SetCanvasSize(750,750);
	gStyle->SetOptStat(0);
	nom->SetLineWidth(3);
	nom->SetLineColor(kBlack);
	up->SetLineWidth(2);
	up->SetLineColor(kBlue);
	down->SetLineWidth(2);
	down->SetLineColor(kRed);

	double max=0.;
	max = TMath::Max(max,nom->GetMaximum());
	max = TMath::Max(max,up->GetMaximum());
	max = TMath::Max(max,down->GetMaximum());

	nom->GetYaxis()->SetRangeUser(0,max*1.1);
	up->GetYaxis()->SetRangeUser(0,max*1.1);
	down->GetYaxis()->SetRangeUser(0,max*1.1);

	nom->GetXaxis()->SetRangeUser(mh_-10,mh_+10);
	up->GetXaxis()->SetRangeUser(mh_-10,mh_+10);
	down->GetXaxis()->SetRangeUser(mh_-10,mh_+10);

	nom->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	up->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	down->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

	nom->SetTitle(Form("%s_%s",name.c_str(),phoCat.c_str()));
	up->SetTitle(Form("%s_%s",name.c_str(),phoCat.c_str()));
	down->SetTitle(Form("%s_%s",name.c_str(),phoCat.c_str()));

	TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->AddEntry(nom,"Nominal","L");
	leg->AddEntry(up,"+3#sigma","L");
	leg->AddEntry(down,"-3#sigma","L");

	down->Draw("HIST");
	up->Draw("HISTsame");
	nom->Draw("HISTsame");
	leg->Draw();
	
	canv->Print(Form("plots/systematics/%s_%s.pdf",name.c_str(),phoCat.c_str()));
	canv->Print(Form("plots/systematics/%s_%s.png",name.c_str(),phoCat.c_str()));
}

double getMeanVar(TH1F* nom, TH1F *up, TH1F* down){
	double u = (up->GetMean()-nom->GetMean())/nom->GetMean();
	double d = (down->GetMean()-nom->GetMean())/nom->GetMean();
	return (TMath::Abs(u)+TMath::Abs(d))/2.;
}

double getSigmaVar(TH1F* nom, TH1F *up, TH1F* down){
	double effSignom = effSigma(nom);
	double effSigup = effSigma(up);
	double effSigdown = effSigma(down);
	double u = (effSigup-effSignom)/effSignom;
	double d = (effSigdown-effSignom)/effSignom;
	return (TMath::Abs(u)+TMath::Abs(d))/2.;
}

double getRateVar(TH1F* nom, TH1F *up, TH1F* down){
	double u = (up->Integral()-nom->Integral())/nom->Integral();
	double d = (down->Integral()-nom->Integral())/nom->Integral();
	return (TMath::Abs(u)+TMath::Abs(d))/2.;
}

int main(int argc, char *argv[]){
 
  OptionParser(argc,argv);

  TStopwatch sw;
  sw.Start();

	system("mkdir -p plots/systematics");
  
	TFile *inFile = TFile::Open(filename_.c_str());
	inFile->ls();
	ofstream outfile;
	outfile.open(outfilename_.c_str());
	outfile << "# this file has been autogenerated by calcPhotonSystConsts.cpp" << endl;
	outfile << endl;
	outfile << "photonCats=" << photonCatStr_ << endl;
	outfile << endl;
	outfile << "# photonCat       mean_change    sigma_change    rate_change" << endl;

	for (int cat=0; cat<nCats_; cat++){
		for (vector<string>::iterator proc=procs_.begin(); proc!=procs_.end(); proc++){
		
			outfile << Form("diphotonCat=%d",cat) << endl;
			outfile << Form("proc=%s",proc->c_str()) << endl;

			TH1F *nominal = (TH1F*)inFile->Get(Form("th1f_sig_%s_mass_m%d_cat%d",proc->c_str(),mh_,cat));
			for (vector<string>::iterator phoCat=photonCats_.begin(); phoCat!=photonCats_.end(); phoCat++){

				TH1F *scaleUp = (TH1F*)inFile->Get(Form("th1f_sig_%s_mass_m%d_cat%d_E_scale_%sUp01_sigma",proc->c_str(),mh_,cat,phoCat->c_str()));
				TH1F *scaleDown = (TH1F*)inFile->Get(Form("th1f_sig_%s_mass_m%d_cat%d_E_scale_%sDown01_sigma",proc->c_str(),mh_,cat,phoCat->c_str()));
			
				plotVariation(nominal,scaleUp,scaleDown,*phoCat,Form("%s_cat%d_scale",proc->c_str(),cat));
				
				outfile << *phoCat+"scale";
				for (unsigned int i=0; i<(15-phoCat->size()); i++) outfile << " ";
				outfile << Form("%4.4f     %4.4f     %4.4f    ",getMeanVar(nominal,scaleUp,scaleDown),getSigmaVar(nominal,scaleUp,scaleDown),getRateVar(nominal,scaleUp,scaleDown)) << endl;

				TH1F *smearUp = (TH1F*)inFile->Get(Form("th1f_sig_%s_mass_m%d_cat%d_E_res_%sUp01_sigma",proc->c_str(),mh_,cat,phoCat->c_str()));
				TH1F *smearDown = (TH1F*)inFile->Get(Form("th1f_sig_%s_mass_m%d_cat%d_E_res_%sDown01_sigma",proc->c_str(),mh_,cat,phoCat->c_str()));
				
				plotVariation(nominal,smearUp,smearDown,*phoCat,Form("%s_cat%d_smear",proc->c_str(),cat));
			
				outfile << *phoCat+"smear";
				for (unsigned int i=0; i<(15-phoCat->size()); i++) outfile << " ";
				outfile << Form("%4.4f     %4.4f     %4.4f    ",getMeanVar(nominal,smearUp,smearDown),getSigmaVar(nominal,smearUp,smearDown),getRateVar(nominal,smearUp,smearDown)) << endl;

			}
			outfile << endl;
		}
	}
	outfile.close();

	return 0;
}
