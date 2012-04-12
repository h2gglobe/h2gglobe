#include <iostream>
#include <string>
#include <vector>

#include "tmvaglob.C"

#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
#include "TText.h"

using namespace std;

std::vector<double> optimizedReverseBinning(TH1F *hb,int nTargetBins,bool revise_target, bool use_n_target){
	// Return a set of bins which are "smoother" 

	if (revise_target) {
		if (use_n_target){
		   std::cerr << "WARNING -- RooContainer::OptimizedBinning -- Can't use number of Entries as target in revised binning algo " << std::endl; 
		   use_n_target = false;  // geometric algo always use revised number of bins, not number of entries
		
		}
	}

	int nBins = hb->GetNbinsX();
	std::vector<double> binEdges;

	double targetNumbers;
	if (use_n_target) targetNumbers = nTargetBins; 
	else targetNumbers = hb->Integral()/nTargetBins;

	if (hb->Integral() < 2*targetNumbers){
		std::cout << "RooContainer::OptimizedBinning -- Not enough entries in histogram for target numbers calculated - " 
			  << targetNumbers 
			  << ", Returning current bin boundaries "  << std::endl;
		//for (int j=2;j<=nBins+1;j++) binEdges.push_back(hb->GetBinLowEdge(j));
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(nBins+1));
		return binEdges;
	}
	binEdges.push_back(hb->GetBinLowEdge(nBins+1));

	std::cout << "RooContainer::optimizedBinning -- Performing Reverse Optimize Binning" <<std::endl;
	double sumBin = 0;
	int i=nBins;
	while (i>=1){
		if (revise_target) targetNumbers = hb->Integral(1,i)/nTargetBins;
		sumBin=hb->GetBinContent(i);
		double highEdge=hb->GetBinLowEdge(i);

		bool carryOn = sumBin <= targetNumbers;
		while ( carryOn){
			if (i>1){
			  sumBin+=hb->GetBinContent(i-1);
			  highEdge = hb->GetBinLowEdge(i-1);
			  carryOn =(sumBin <targetNumbers && i>=1);
			  i--;
			} else {
			  highEdge = hb->GetBinLowEdge(i);
			  carryOn=false;
			}
		}
	        binEdges.push_back(highEdge);
		i--;
	}
        if (sumBin < 10) binEdges.erase(binEdges.end()-2);
	reverse(binEdges.begin(),binEdges.end());
	return binEdges;

}

void histogramSmoothing(TH1F* h, int n){
   // Nothing too special, a function which will smooth a histogram but ignore the first and last
   // bins, useful for the "flat-binning" approach! 
   if (h->GetNbinsX()>3){
     int nbin = h->GetNbinsX();
     TH1F *h2 = new TH1F(Form("hn%s",h->GetName()),Form("hn%s",h->GetName()),nbin-2,0,1);
     for (int i=1;i<=nbin-2;i++){
           h2->SetBinContent(i,h->GetBinContent(i+1));
     }
     h2->Smooth(n);
     for (int i=2;i<=nbin-1;i++){
           h->SetBinContent(i,h2->GetBinContent(i-1));
     }

   }

   return;

   
}

TH1F* rebinBinnedDataset(std::string new_name,std::string name, TH1F *hb,std::vector<double> binEdges){

	double *arrBins = new double[binEdges.size()];
	int j=0;
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		arrBins[j]=*it;
		j++;
		
	}
	//const char *h_name = (const char *) hb->GetName;
	//const char *title  = (const char *) hb->GetTitle;
	
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,hb->GetName(),arrBins);
	hbnew->SetName(Form("th1f_%s",new_name.c_str()));
	//cout << "title for new re-binned histogream - " << hb->GetTitle()<<endl; 
	hbnew->SetTitle(hb->GetTitle());

	// Just a quick test, mask the last "channel"
	//hbnew->SetBinContent(hbnew->GetNbinsX(),0);
	//cout << "DONT DO THIS IN MAIN PROGRAM ----- LINE 1563 rebin setting last bin to 0" <<endl;
	delete [] arrBins;
  return hbnew;	
}

std::vector<double> soverBOptimizedBinning(TH1F *hs,TH1F *hb,int nTargetBins,double penaltyScale){

	// Performs Optimized Binning based on a Signal and Background (S/B) distributions
	// First runs the optimizedBinning on background and rebins S and B clones, note, always performs 
	// revise_target=false,direction=-1 and use_n_entries=true
	// nTargetBins is used for the flat binning, decision to merge based on penaltyScale
	int ninitBins = hb->GetNbinsX();
	if (hs->Integral()==0 ||  hb->Integral()==0 || ninitBins < 2) {
		std::vector<double> binEdges;
		binEdges.push_back(hb->GetBinLowEdge(1));
		binEdges.push_back(hb->GetBinLowEdge(ninitBins+1));
		return binEdges;
	}

	std::vector<double> binEdges = optimizedReverseBinning(hb,nTargetBins,false,true);

	int j =0;
	double *arrBins = new double[binEdges.size()];
	for (std::vector<double>::iterator it=binEdges.begin();it!=binEdges.end();it++){
		//cout << *it << endl;
		arrBins[j]=*it;
		j++;	
	}
	// Create new rebinned histograms (only temporary)
	TH1F *hbnew =(TH1F*) hb->Rebin(binEdges.size()-1,"hbnew",arrBins);
	TH1F *hsnew =(TH1F*) hs->Rebin(binEdges.size()-1,"hsnew",arrBins);
	

	// Better smoothing which doesn't use the first and last bins	
	if (hsnew->Integral()!=0 && hbnew->Integral()!=0 && binEdges.size()-1 > 10){
		histogramSmoothing(hsnew,1000);
		histogramSmoothing(hbnew,1000);
		//hsnew->Smooth(1000);
		//hbnew->Smooth(1000);
        }
	// Do we really need the background histogram ?  we will be assuming that the first step is nentries per bin

	// Smooth signal new binned histograms, the size of smoothing should be ~1% of the total bins	
	//int nSmooth = (int) 0.01*hsnew->GetNbinsX();
	//hsnew->Smooth(nSmooth);

	delete [] arrBins;

	if (hbnew->Integral()==0 || hsnew->Integral()==0) return binEdges;
	std::vector<double> newbinEdges;
	newbinEdges.push_back(hbnew->GetBinLowEdge(1));
	int nNewBins = hbnew->GetNbinsX();
	int i=1;
	double maxSoB=0;
        for (int j=1;j<=hbnew->GetNbinsX();j++){
                double newMaximum = hsnew->GetBinContent(j)/hbnew->GetBinContent(j) ;
                if (newMaximum>maxSoB) maxSoB=newMaximum;
        }

	while (i<=nNewBins){

		int k = i+1;
		double highEdge=hbnew->GetBinLowEdge(i+1);
		double S = hsnew->GetBinContent(i);
		double B = hbnew->GetBinContent(i);
                double Stot =S;
                double Btot =B;

		if (B!=0){ 
		  bool carryOn=true;

		  while ( carryOn){

		  	double SoB = S/B;

			if (k<=nNewBins){

			  double S1 = hsnew->GetBinContent(k);
			  double B1 = hbnew->GetBinContent(k);
			  if (B1==0) {
				carryOn=true;
			      	highEdge = hbnew->GetBinLowEdge(k+1);
				Stot+=S1;
				k++;
			  }
			  else{

			    double SoB1 = S1/B1;
			    double scaler=penaltyScale;
                            double prob   = TMath::Prob((scaler*(Stot+S1))*(scaler*(Stot+S1))/(Btot+B1),1);
                            double prob1  = TMath::Prob(((scaler*Stot)*(scaler*Stot)/(Btot)) + ((scaler*S1)*(scaler*S1)/(B1)),2);
                  	    double importance = SoB/maxSoB;

//			      if (fabs(SoB-SoB1)/SoB < penaltyScale/importance ){
			      if ( prob<prob1 ){

			      highEdge = hbnew->GetBinLowEdge(k+1);
			      Stot+=S1;
			      Btot+=B1;
			      carryOn = true;
			      k++;
			    } else {carryOn=false;}
			 }

			} else {
			  highEdge = hbnew->GetBinLowEdge(k);
			  carryOn=false;
			}
		  }
		}
	        newbinEdges.push_back(highEdge);
		i=k;
	}

	// now we have new Bin edges to return to the 
	return newbinEdges;

}

void FillHist(TH1F* hist, string sorb, TTree* tree, string method){

  int SorB=-1;
  if (sorb=="Signal") SorB=0;
  else if (sorb=="Background") SorB=1;
  else {
    cout << "ERROR: must specify Signal or Background" << endl;
    exit(1);
  }
  if (method!="Likelihood" && method!="LikelihoodD" && method!="BDTada" && method!="BDTgrad"){
    cout << "ERROR: " << method << " is not a valid method" << endl;
    exit(1);
  }
  if (SorB!=0 && SorB!=1){
    cout << "ERROR: " << endl;
    exit(1);
  }

  int classID;
  float BDTadaMIT, BDTgradMIT, LikelihoodMIT, LikelihoodDMIT,weight;
  tree->SetBranchAddress("weight",&weight);
  tree->SetBranchAddress("classID",&classID);
  tree->SetBranchAddress("BDTadaMIT",&BDTadaMIT);
  tree->SetBranchAddress("BDTgradMIT",&BDTgradMIT);
  tree->SetBranchAddress("LikelihoodMIT",&LikelihoodMIT);
  tree->SetBranchAddress("LikelihoodDMIT",&LikelihoodDMIT);

  for (int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if (classID==SorB){
      if (method=="BDTada") hist->Fill(BDTadaMIT,weight);
      if (method=="BDTgrad") hist->Fill(BDTgradMIT,weight);
      if (method=="Likelihood") hist->Fill(LikelihoodMIT,weight);
      if (method=="LikelihoodD") hist->Fill(LikelihoodDMIT,weight);
    }
  }

}

TH1F* neatBin(TH1F* hist){

  string name(hist->GetName());
  TH1F *newHist = new TH1F(name.c_str(),name.c_str(),hist->GetNbinsX(),0,hist->GetNbinsX());
  for (int bin=1; bin<=hist->GetNbinsX(); bin++){
    newHist->SetBinContent(bin,hist->GetBinContent(bin));
  }
  return newHist;
}

void Rebin(){
  
  TFile *inFile = TFile::Open("07_01_12_MIT_all.root");

  TTree *testTree = (TTree*)inFile->FindObjectAny("TestTree");
  TTree *trainTree = (TTree*)inFile->FindObjectAny("TrainTree");

  const int nMeths=4;
  string methods[nMeths] = {"Likelihood","LikelihoodD","BDTada","BDTgrad"};
  //string methods[nMeths] = {"Likelihood","LikelihoodD"};

  for (int m=0; m<nMeths; m++){
    
    TH1F *sig = new TH1F(Form("sig%s",methods[m].c_str()),Form("sig%s",methods[m].c_str()),5000,-1.,1.);
    TH1F *bgd = new TH1F(Form("bgd%s",methods[m].c_str()),Form("bgd%s",methods[m].c_str()),5000,-1.,1.);
    TH1F *sigOv = new TH1F(Form("sigOv%s",methods[m].c_str()),Form("sigOv%s",methods[m].c_str()),5000,-1.,1.);
    TH1F *bgdOv = new TH1F(Form("bgdOv%s",methods[m].c_str()),Form("bgdOv%s",methods[m].c_str()),5000,-1.,1.);

    FillHist(sig,"Signal",testTree,methods[m]);
    FillHist(bgd,"Background",testTree,methods[m]);
    FillHist(sigOv,"Signal",trainTree,methods[m]);
    FillHist(bgdOv,"Background",trainTree,methods[m]);

    TMVAGlob::SetSignalAndBackgroundStyle( sig, bgd );
    sig->Scale(72.39/sig->Integral());
    bgd->Scale(3938.7/bgd->Integral());
    sigOv->Scale(72.39/sigOv->Integral());
    bgdOv->Scale(3938.7/bgdOv->Integral());

    Int_t col = sig->GetLineColor();
    sigOv->SetMarkerColor( col );
    sigOv->SetMarkerSize( 0.7 );
    sigOv->SetMarkerStyle( 20 );
    sigOv->SetLineWidth( 1 );
    sigOv->SetLineColor( col );

    col = bgd->GetLineColor();
    bgdOv->SetMarkerColor( col );
    bgdOv->SetMarkerSize( 0.7 );
    bgdOv->SetMarkerStyle( 20 );
    bgdOv->SetLineWidth( 1 );
    bgdOv->SetLineColor( col );

    std::vector<double> binEdges = soverBOptimizedBinning(sigOv,bgdOv,20,50);
    TH1F *sigReb = rebinBinnedDataset(Form("sig_reb_%s",methods[m].c_str())," ",sig,binEdges);
    TH1F *bgdReb = rebinBinnedDataset(Form("bgd_reb_%s",methods[m].c_str())," ",bgd,binEdges);
    TH1F *sigOvReb = rebinBinnedDataset(Form("sig_train_reb_%s",methods[m].c_str())," ",sigOv,binEdges);
    TH1F *bgdOvReb = rebinBinnedDataset(Form("bgd_train_reb_%s",methods[m].c_str())," ",bgdOv,binEdges);

    sigReb->GetYaxis()->SetRangeUser(1,ceil(TMath::Max(sigReb->GetMaximum(),bgdReb->GetMaximum())*1.1));
    TCanvas *c = new TCanvas();
    c->cd();
    TMVAGlob::NormalizeHists( sig, bgd );
    TMVAGlob::NormalizeHists( sigOv, bgdOv );
    sig->Rebin(50);
    bgd->Rebin(50);
    sigOv->Rebin(50);
    bgdOv->Rebin(50);
    sig->GetYaxis()->SetRangeUser(0.,ceil(TMath::Max(sig->GetMaximum(),bgd->GetMaximum())*1.1));
    sig->Draw("hist");
    bgd->Draw("samehist");
    sigOv->Draw("e1same");
    bgdOv->Draw("e1same");
    Double_t kolS_ = sig->KolmogorovTest( sigOv );
    Double_t kolB_ = bgd->KolmogorovTest( bgdOv );
    TString probatext_ = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS_, kolB_ );
    TText* tt_ = new TText( 0.12, 0.84, probatext_ );
    tt_->SetNDC(); tt_->SetTextSize( 0.032 ); 
    tt_->Draw("same");
    c->Print(Form("plots/NormBin_%s.pdf",methods[m].c_str()),"pdf");
   
    TH1F *sigNice = neatBin(sigReb);
    TH1F *bgdNice = neatBin(bgdReb);
    TH1F *sigOvNice = neatBin(sigOvReb);
    TH1F *bgdOvNice = neatBin(bgdOvReb);

    Double_t kolS = sigNice->KolmogorovTest( sigOvNice );
    Double_t kolB = bgdNice->KolmogorovTest( bgdOvNice );
    TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
    TText* tt = new TText( 0.12, 0.84, probatext );
    tt->SetNDC(); tt->SetTextSize( 0.032 ); 
    
    sigNice->Scale(5.);
    sigOvNice->Scale(5.);

    //sigNice->GetYaxis()->SetRangeUser(1,ceil(TMath::Max(sigNice->GetMaximum(),bgdNice->GetMaximum())+8000));
    sigNice->GetYaxis()->SetRangeUser(1,10000);
    sigNice->GetXaxis()->SetTitle(Form("%s Output Bin Number",methods[m].c_str()));
    sigNice->SetTitle(Form("%s",methods[m].c_str()));

    TMVAGlob::SetSignalAndBackgroundStyle( sigNice, bgdNice );
    Int_t col = sig->GetLineColor();
    sigOvNice->SetMarkerColor( col );
    sigOvNice->SetMarkerSize( 0.7 );
    sigOvNice->SetMarkerStyle( 20 );
    sigOvNice->SetLineWidth( 1 );
    sigOvNice->SetLineColor( col );
    col = bgd->GetLineColor();
    bgdOvNice->SetMarkerColor( col );
    bgdOvNice->SetMarkerSize( 0.7 );
    bgdOvNice->SetMarkerStyle( 20 );
    bgdOvNice->SetLineWidth( 1 );
    bgdOvNice->SetLineColor( col );
    
    TLegend *legend= new TLegend( 0.6,0.5,0.9,0.8 ); 
    legend->SetFillStyle( 0 );
    legend->AddEntry(sig,TString("Signal") + "(test sample)", "F");
    legend->AddEntry(bgd,TString("Background") + "(test sample)", "F");
    legend->SetBorderSize(1);
    legend->AddEntry(sigOv,"Signal (training sample)","lep");
    legend->AddEntry(bgdOv,"Background (training sample)","lep");
    legend->Draw("same");


    TCanvas *c2 = new TCanvas();
    c2->cd();
    gStyle->SetOptStat(0);
    c2->SetLogy();
    sigNice->Draw("hist");
    bgdNice->Draw("samehist");
    sigOvNice->Draw("e1same");
    bgdOvNice->Draw("e1same");
    tt->Draw("same");
    legend->Draw("same");
    c2->Print(Form("plots/ICbinning_%s.pdf",methods[m].c_str()),"pdf");
  }

}
