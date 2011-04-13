/*
 * TestReweight.C
 *
 *  Created on: Apr 8, 2011
 *      Author: adavid
 */

TH1 *MCGen;
TH1 *WeightsGenPoisson;
TH1 *TargetGenPoisson;
TH1 *GenTrial;

TH2 *MC2DTruth;
TH2 *MC2DResponse;
TH1 *MC2DGen;
TH1 *Weights2DPoisson;
TH1 *Target2DPoisson;
TH1 *Trial2D;

TCanvas *c;
TCanvas *l;

void TestReweight(int n=50, int mean=15,int meantarget = 10){

	ReweightMC::target_t mode;
	int halfway = (n+mean)/2;

	//
	mode = 0;
	//
	ReweightMC *rPassMCGen = new ReweightMC(mode, n);
	for(int gen=0; gen<=n; ++gen){
		double weight = 1-TMath::Gaus(gen,halfway,0.5,true);
		if(gen<halfway/2) weight *= (float)(gen)/halfway/2.;
		if(gen<1) weight = 0;
		rPassMCGen->fillMC(gen, weight);
	}
	rPassMCGen->makeResponse();
	MCGen = (TH1*) rPassMCGen->getResponse()->Clone();
	delete rPassMCGen;
	//
	ReweightMC *rPassPoissonGen = new ReweightMC(MCGen, meantarget);
	WeightsGenPoisson = (TH1*) rPassPoissonGen->getWeights()->Clone();
	TargetGenPoisson = (TH1*) rPassPoissonGen->getTarget()->Clone();
	delete rPassPoissonGen;
		
	//
	mode = 1;
	//
	ReweightMC *rPassMCGenReco = new ReweightMC(mode, n);
	TH2 * MC2DModel = (TH2 *) rPassMCGenReco->getResponse()->Clone();
	for(int reco=0; reco<=n; ++reco){
		for(int gen=0; gen<=n; ++gen){

			double weight=1;
//			if (reco<mean) weight = 2;

//			double weight= (reco+1)*gen;

// 			double weight=1;
// 			if (gen<mean) weight = 2;

// 			double weight=1;
// 			if (reco<gen) weight = 2;

//			double weight=1;
//			if (reco<gen-3) weight = 2;
//			if (reco>gen+3) weight = 0;
//			if (reco<gen-6) weight = 0;

//			double weight=0;
//			if (reco<=gen+1){
//				weight = TMath::Gaus(reco,gen,(gen+1.)/2.,true);
//				if(gen>=halfway) weight *= ( (float)(n-gen)/(n-halfway));
//			}

//  			double weight=0;
//  			weight = TMath::Gaus(reco,gen,(gen+1.)/2.,true)
//  					* TMath::PoissonI(gen,mean);
//			weight = TMath::Gaus(reco,gen-0.001*gen*gen,2.,true) * TMath::PoissonI(gen,mean);
			weight = TMath::Gaus(reco,gen,5.,true) * TMath::PoissonI(gen,mean);
			
			/// MC2DModel->Fill(gen, reco, 1e+9*weight);
			rPassMCGenReco->fillMC(gen, reco, weight);
		}
	}
	/// rPassMCGenReco->getResponse()->FillRandom(MC2DModel,1000000000);
	MC2DTruth = (TH2 *) rPassMCGenReco->getResponse()->Clone();
	rPassMCGenReco->makeResponse();
	MC2DResponse = (TH2 *) rPassMCGenReco->getResponse()->Clone();
	MC2DGen = (TH1 *) rPassMCGenReco->getMCGen()->Clone();
	//
//	for(int reco=3; reco<=7; ++reco){
//		rPassMCGenReco->fillTarget(reco);
//	}
//	Target2DPoisson = (TH1 *) rPassMCGenReco->getTarget()->Clone();
	delete rPassMCGenReco;


//	ReweightMC *rPassPoissonGenReco = new ReweightMC(MC2DResponse, MC2DGen, Target2DPoisson);
	ReweightMC *rPassPoissonGenReco = new ReweightMC(MC2DResponse, MC2DGen, meantarget);
	Weights2DPoisson = (TH1*) rPassPoissonGenReco->getWeights()->Clone();
	Target2DPoisson = (TH1*) rPassPoissonGenReco->getTarget()->Clone();
	delete rPassPoissonGenReco;



	//
	//
	c = new TCanvas("checks", "Reweight checks", 900, 600 );
	c->Divide(3,2);

	c->cd(1);
	MCGen->Draw();

	c->cd(2);
	WeightsGenPoisson->Draw();

	c->cd(3);
	TargetGenPoisson->Draw();
	GenTrial = (TH1*) TargetGenPoisson->Clone("temp1");
	GenTrial->Reset("ICE");
	GenTrial->SetEntries(0);
	GenTrial->Multiply(MCGen, WeightsGenPoisson);
	GenTrial->SetMarkerColor(kBlue);
	GenTrial->Draw("same");

	c->cd(4);
	MC2DResponse->Draw("colz");
	/// MC2DGen->DrawNormalized("hist");
	/// ((TH2*)MC2DResponse)->ProjectionX("a")->DrawNormalized("samehistp");
	/// /// ((TH2*)MC2DTruth)->ProjectionY("b")->DrawNormalized("samehist");

	c->cd(5);
	Weights2DPoisson->Draw("histp");

	c->cd(6);

	//

//	TargetGenRecoPoisson->Draw("hist");
//	GenRecoTrial = (TH1*) TargetGenRecoPoisson->Clone("temp2");
//	GenRecoTrial->Reset("ICE");
//	GenRecoTrial->SetEntries(0);
//	for(int r=0; r<MC2DResponse->GetNbinsY(); ++r){
//		double sum = 0;
//		for(int g=0; g<MC2DResponse->GetNbinsX(); ++g){
//			sum += MC2DResponse->GetBinContent(MC2DResponse->FindBin(g,r))
//					* Weights2DPoisson->GetBinContent(Weights2DPoisson->FindBin(g));
//		}
//		Trial2D->Fill(r, sum);
//	}
//	Trial2D->SetMarkerColor(kBlue);
//	Trial2D->Draw("same");

	//

	TH2 *h2 = (TH2*) MC2DResponse->Clone("temp2");
	h2->Reset("ICE");
	h2->SetEntries(0);

	for(int r=0; r<MC2DTruth->GetNbinsY(); ++r){
		for(int g=0; g<MC2DTruth->GetNbinsX(); ++g){
		h2->Fill(g, r,
			 MC2DTruth->GetBinContent(MC2DResponse->FindBin(g,r))
			 *
			 Weights2DPoisson->GetBinContent(Weights2DPoisson->FindBin(g))
				);
		}
	}
	h2->Draw("colz");
	Trial2D = (TH1 *)h2->ProjectionY()->Clone();

	l = new TCanvas();
	Trial2D->Print();
	Trial2D->DrawNormalized("histp");
	Target2DPoisson->Print();
	Target2DPoisson->DrawNormalized("samehist");

	
	/// TH2 *A = MC2DTruth;
	/// TH1 *y = Target2DPoisson;
	/// TUnfold unfold(A,TUnfold::kHistMapOutputVert);
	/// unfold.SetInput(y);
	/// Int_t nScan=30;
	/// Int_t iBest;
	/// TSpline *logTauX,*logTauY;
	/// TGraph *lCurve;
	/// iBest=unfold.ScanLcurve(nScan,0.0,0.0,&lCurve);
	/// std::cout<<"tau="<<unfold.GetTau()<<"\n";
	/// TH1D *x=unfold.GetOutput("x","myVariable");
	/// TH2D *rhoij=unfold.GetRhoIJ("correlation","myVariable");
	/// 
	/// m = new TCanvas();
	/// x->Draw("histp");
}
