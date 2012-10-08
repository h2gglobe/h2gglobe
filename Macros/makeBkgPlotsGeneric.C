// Makes partially blinded mass distribution + fit plots for Mass-fac MVA analysis

#ifndef __CINT__
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include <iostream>
#endif 

void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname);

void makeBkgPlotsGeneric(std::string filebkg, bool blind=true, bool doBands=true, bool baseline=false, bool useBinnedData=false){

	// Globals
	gROOT->SetStyle("Plain");
	gROOT->SetBatch(1);
	gStyle->SetOptStat(0);
	const int ncats = 9;

	RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::WARNING));

	std::string * labels;
	
	std::string baselinelabels[6] = { "Both photons in barrel, R_{9}^{min} > 0.94"
					  ,"Both photons in barrel, R_{9}^{min} < 0.94"
					  ,"One or more photons in endcap, R_{9}^{min} > 0.94"
					  ,"One or more photons in endcap, R_{9}^{min} < 0.94"
					  ,"Dijet-tagged class m_{jj} > 500 GeV"
					  ,"Dijet-tagged class 250 < m_{jj} < 500 GeV"
	};
	std::string massfactlabels[ncats] = { 
					"BDT_{#gamma#gamma} >= 0.88"
					,"0.71  <= BDT_{#gamma#gamma} < 0.88"
					,"0.5 <= BDT_{#gamma#gamma} < 0.71"
					,"-0.05  <= BDT_{#gamma#gamma} < 0.5"
					,"Dijet-tagged class BDT_{#gamma#gamma}_{jj} > 0.93"
					,"Dijet-tagged class BDT_{jj} > 0.98"
          ,"Muon-tagged class"
          ,"Electron-tagged class"
          ,"MET-tagged class"
	};
	
	if( baseline ) { labels = baselinelabels; }
	else { labels = massfactlabels; }

	TFile *fb = TFile::Open(filebkg.c_str());
	
	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
//	w_bkg->Print();

	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
  RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");

	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.025);
	latex->SetNDC();
	
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.04);
	cmslatex->SetNDC();

	double totalGGHinDIJET = 0;
	double totalVBFinDIJET = 0;
	double totalGGHinINCL = 0;
	double totalVBFinINCL = 0;

	double totalTTHinDIJET = 0;
	double totalWZHinDIJET = 0;
	double totalTTHinINCL = 0;
	double totalWZHinINCL = 0;

	for (int cat=0;cat<ncats;cat++){
		
		TCanvas *can = new TCanvas("c","",800,800);
		TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);

		// Get Dataset ->
		RooAbsData *data;
		if (useBinnedData) data =  (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
		else  data =  (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
		data->Print();

		// Background Pdf ->
		/// RooExtendPdf *bkg =  (RooExtendPdf*)w_bkg->pdf(Form("data_pol_model_8TeV_cat%d",cat));
		RooAbsPdf *bkg =  (RooAbsPdf*)w_bkg->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
		bkg->Print();
		bkg->fitTo(*data);
		RooFitResult *r = bkg->fitTo(*data,RooFit::Save(1));
		
		// Get Signal pdf norms
		std::cout << "Getting Signal Components" << std::endl;
		TH1F *gghnorm = (TH1F*)fb->Get(Form("th1f_sig_ggh_mass_m125_cat%d",cat));
		TH1F *vbfnorm = (TH1F*)fb->Get(Form("th1f_sig_vbf_mass_m125_cat%d",cat));
		TH1F *wzhnorm = (TH1F*)fb->Get(Form("th1f_sig_wzh_mass_m125_cat%d",cat));
		TH1F *tthnorm = (TH1F*)fb->Get(Form("th1f_sig_tth_mass_m125_cat%d",cat));
		
		if (cat<=3){
			totalGGHinINCL+=gghnorm->Integral();
			totalVBFinINCL+=vbfnorm->Integral();
			totalTTHinINCL+=tthnorm->Integral();
			totalWZHinINCL+=wzhnorm->Integral();
		}else{
			totalGGHinDIJET+=gghnorm->Integral();
			totalVBFinDIJET+=vbfnorm->Integral();
			totalTTHinDIJET+=tthnorm->Integral();
			totalWZHinDIJET+=wzhnorm->Integral();
		}
		
		std::cout << "Rescaling Signal Components" << std::endl;
		gghnorm->Add(vbfnorm);
		gghnorm->Add(wzhnorm);
		gghnorm->Add(tthnorm);


		TH1F *allsig = (TH1F*)gghnorm->Clone();
		allsig->Rebin(2);
		allsig->SetLineColor(4);allsig->SetFillColor(38);allsig->SetFillStyle(3001) ;allsig->SetLineWidth(2);
		/// allsig->SetLineColor(1);
		/// allsig->SetFillColor(38);
		TH1F dumData("d","",80,100,180); dumData.Sumw2();dumData.SetMarkerSize(1.0);dumData.SetMarkerStyle(20);dumData.SetLineWidth(3);
		dumData.Fill(101);
		// TH1F dumSignal("s","",80,100,180); dumSignal.SetLineColor(4);dumSignal.SetFillColor(38);dumSignal.SetFillStyle(3001) ;dumSignal.SetLineWidth(2);
		TH1F dum1Sig("1s","",80,100,180); dum1Sig.SetFillColor(kYellow);dum1Sig.SetFillStyle(1001);
		TH1F dum2Sig("2s","",80,100,180); dum2Sig.SetFillColor(kGreen);dum2Sig.SetFillStyle(1001);
		TH1F dumBkg("b","",80,100,180); dumBkg.SetLineColor(kRed);dumBkg.SetLineWidth(3);
		dumBkg.Draw("P");// dumSignal.Draw("LFsame");
		dumBkg.Draw("Fsame");dum1Sig.Draw("Fsame");dum2Sig.Draw("Lsame");

		// Plot background
		RooPlot *frame = x->frame();

		std::cout << "Plotting Components" << std::endl;
		data->plotOn(frame,RooFit::Binning(80),RooFit::Invisible());
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,2,1),RooFit::FillColor(kGreen));
		/// bkg->plotOn(frame,RooFit::VisualizeError(*r,1,1),RooFit::FillColor(kYellow));
		bkg->plotOn(frame,RooFit::LineColor(kRed));
		TGraphAsymmErrors *onesigma = 0, *twosigma = 0;
		if( doBands ) {
			onesigma = new TGraphAsymmErrors();
			twosigma = new TGraphAsymmErrors();
			doBandsFit(onesigma, twosigma, x, bkg, dynamic_cast<RooCurve*>(frame->getObject(frame->numItems()-1)), 
				   data, frame, Form("cat%d",cat) );
		}
		if( blind ) {
			x->setRange("unblind_up",150,180);
			data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_up"));
			x->setRange("unblind_down",100,110);
			data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_down"));
		} else {
			data->plotOn(frame,RooFit::Binning(80));
		}
		///// x->setRange("unblind_1",100,110);
		///// x->setRange("unblind_2",150,180);
		///// data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_1"));
		///// data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_2"));

		frame->SetTitle("");
		frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
		frame->GetXaxis()->SetNdivisions(5,5,0);
		frame->GetYaxis()->SetTitle("Events / (1 GeV)");
		frame->GetYaxis()->SetTitleOffset(1.2);
		
		leg->AddEntry(&dumData,"Data","PEL");
		leg->AddEntry(&dumBkg,"Bkg Model","L");
		leg->AddEntry(&dum1Sig,"#pm 1#sigma","F");
		leg->AddEntry(&dum2Sig,"#pm 2#sigma","F");
		/// leg->AddEntry(&dumSignal,"1xSM m_{H} = 125 GeV","F");
		leg->AddEntry(allsig,"1xSM m_{H} = 125 GeV","F");
		
		frame->Draw();
		frame->SetMinimum(0.0001);
 		if( doBands ) {
 			twosigma->SetLineColor(kGreen);
 			twosigma->SetFillColor(kGreen);
 			twosigma->SetMarkerColor(kGreen);
 			twosigma->Draw("L3 SAME");     
 			
 			onesigma->SetLineColor(kYellow);
 			onesigma->SetFillColor(kYellow);
 			onesigma->SetMarkerColor(kYellow);
 			onesigma->Draw("L3 SAME");
 			frame->Draw("same");
 		}
		allsig->Draw("samehistF");
		leg->Draw();
		cmslatex->DrawLatex(Form(0.15,0.8,"#splitline{CMS Preliminary}{#sqrt{s} = 8TeV L = %2.1fb^{-1}}",intL->getVal()));
		latex->DrawLatex(0.1,0.92,labels[cat].c_str());
		can->SaveAs(Form( (baseline ? "baselinecat%d.pdf" : "massfacmvacat%d.pdf"),cat));
		can->SaveAs(Form( (baseline ? "baselinecat%d.png" : "massfacmvacat%d.png"),cat));
	}

	// JET ID Systematics
	std::cout << "The following can be used (nick knows what they mean) as the JET Migration systematics" <<std::endl;
	std::cout << "XXX " << 1-(0.7*totalGGHinDIJET/(totalGGHinINCL))<<std::endl;
	std::cout << "YYY " << 1-(0.7*totalTTHinDIJET/(totalTTHinINCL))<<std::endl;
	std::cout << "MMM " << 1-(0.1*totalVBFinDIJET/(totalVBFinINCL))<<std::endl;
	std::cout << "NNN " << 1-(0.1*totalWZHinDIJET/(totalWZHinINCL))<<std::endl;
}


using namespace RooFit;


void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname)
{
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catname.Data()),"",0.0,0.0,1e+5.0);

	for (int i=1; i<(plot->GetXaxis()->GetNbins()+1); ++i) {

		double lowedge = plot->GetXaxis()->GetBinLowEdge(i);
		double upedge = plot->GetXaxis()->GetBinUpEdge(i);
		double center = plot->GetXaxis()->GetBinCenter(i);
        
		double nombkg = nomcurve->interpolate(center);

		nlim->setVal(nombkg);
		hmass->setRange("errRange",lowedge,upedge);
		RooAbsPdf *epdf = 0;
		epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        
		RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(4));
		RooMinimizer minim(*nll);
		minim.setStrategy(0);
		minim.setPrintLevel(-1);
		double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
		double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        
		minim.migrad();
		minim.minos(*nlim);
		
		printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
		
		onesigma->SetPoint(i-1,center,nombkg);
		onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        
		minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
		// eventually if cl = 0.95 this is the usual 1.92!      
        	minim.migrad();
		minim.minos(*nlim);
		
		twosigma->SetPoint(i-1,center,nombkg);
		twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());      
        		
		delete nll;
		delete epdf;
	}

	onesigma->Print("V");

}
