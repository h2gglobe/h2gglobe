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

void TestReweight(int n=20, int mean=10,int meantarget = 8){

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
	RooUnfoldResponse response(n+1, -0.5, n+0.5);
	ReweightMC *rPassMCGenReco = new ReweightMC(mode, n);
	TH2 * MC2DModel = (TH2 *) rPassMCGenReco->getResponse()->Clone();
//	for(int reco=0; reco<=n; ++reco){
		for(int gen=1; gen<=n; ++gen){

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
//			weight = TMath::Gaus(reco,gen-0.005*gen*gen,2.,true) * TMath::PoissonI(gen,mean);
//			weight = TMath::Gaus(reco,gen,2.,true) * TMath::PoissonI(gen,mean);

			for(int nev= 10000*TMath::PoissonI(gen,mean); nev--;){

				double reco = (int) gRandom->Gaus(gen-0.01*gen*gen,0.1);
				if(reco <= 0){
					response.Miss(gen);
				}
				else{
					response.Fill(reco, gen);
					MC2DModel->Fill(gen, reco);
				}
				//rPassMCGenReco->fillMC(gen, reco, weight);
			}
		}
//	}



	/// rPassMCGenReco->getResponse()->FillRandom(MC2DModel,1000000000);
//	MC2DTruth = (TH2 *) rPassMCGenReco->getResponse()->Clone();
//	rPassMCGenReco->makeResponse();
//	MC2DResponse = (TH2 *) rPassMCGenReco->getResponse()->Clone();
//	MC2DGen = (TH1 *) rPassMCGenReco->getMCGen()->Clone();
//	//
////	for(int reco=3; reco<=7; ++reco){
////		rPassMCGenReco->fillTarget(reco);
////	}
////	Target2DPoisson = (TH1 *) rPassMCGenReco->getTarget()->Clone();
//	delete rPassMCGenReco;
//
//
////	ReweightMC *rPassPoissonGenReco = new ReweightMC(MC2DResponse, MC2DGen, Target2DPoisson);
//	ReweightMC *rPassPoissonGenReco = new ReweightMC(MC2DResponse, MC2DGen, meantarget);
//	Weights2DPoisson = (TH1*) rPassPoissonGenReco->getWeights()->Clone();
//	Target2DPoisson = (TH1*) rPassPoissonGenReco->getTarget()->Clone();
//	delete rPassPoissonGenReco;

	//
	//
	c = new TCanvas("checks", "Reweight checks", 900, 600 );
	c->Divide(3,2);

	c->cd(1);
//	MCGen->Draw();

	c->cd(2);
//	WeightsGenPoisson->Draw();

	c->cd(3);
//	TargetGenPoisson->Draw();
	GenTrial = (TH1*) TargetGenPoisson->Clone("temp1");
	GenTrial->Reset("ICE");
	GenTrial->SetEntries(0);
	GenTrial->Multiply(MCGen, WeightsGenPoisson);
	GenTrial->SetMarkerColor(kBlue);
//	GenTrial->Draw("same");


	c->cd(4);
	//MC2DModel->Draw("colz");
	//response.Hresponse()->Draw("colz");
	response.Htruth()->Draw("h");
	response.Hmeasured()->Draw("pesame");

	c->cd(5);
	RooUnfoldBayes unfold (&response, Target2DPoisson, 5);
	//RooUnfoldTUnfold unfold (&response, Target2DPoisson);
	TH1* unfolded = unfold.Hreco();
	unfolded->Draw();
	c->cd(6);
	Target2DPoisson->Draw();
	TH1* transformed = response.ApplyToTruth(unfolded);
	transformed->SetMarkerColor(kBlue);
	transformed->Draw("pesame");


}

