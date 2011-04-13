#include "../interface/ReweightMC.h"

#include <assert.h>
#include "TMatrixD.h"
#include "TMath.h"
#include "TClass.h"
#include <iostream>
#include <algorithm>

#define assert_inrange(X) assert(0 <= X && X <= nVtx_)

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
const TString ReweightMC::UID("ReweightMC_");

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::ReweightMC(target_t target, int nVtx):
	target_(target),
	nVtx_(nVtx),
	hWeights_(0)
{
	TH1::SetDefaultSumw2(true);

	switch(target_){
	case Gen:
		hMCGenOrGenReco_ = h_.MCGen = new TH1D(UID+"hMCGen","Number of Gen vertices in MC; nGenMC; A.U.",
				nVtx_+1,-0.5,nVtx_+0.5);
		hTarget_ = new TH1D(UID+"nGenTarget","Number of Target gen vertices; nGenTarget",
				nVtx_+1,-0.5,nVtx_+0.5);
		break;
	case Reco:
		hMCGenOrGenReco_ = h_.MCGenReco = new TH2D(UID+"hMCGenReco","Smearing between number of Gen and Reco vertices in MC; nGenMC; nRecoMC",
				nVtx_+1,-0.5,nVtx_+0.5,nVtx_+1,-0.5,nVtx_+0.5);
		hTarget_ = new TH1D(UID+"nRecoTarget","Number of Target reco vertices; nRecoTarget",
				nVtx_+1,-0.5,nVtx_+0.5);
		break;
	}

	toDelete_.push_back(hMCGenOrGenReco_);
	toDelete_.push_back(hTarget_);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::ReweightMC(TH2 *hMCGenReco, TH1 *hMCGen, float RecoPoissonMean):
	target_(Reco),
	hMCGenOrGenReco_(hMCGenReco),
	hMCGen_(hMCGen)
{
	assert(hMCGenOrGenReco_ != 0);
	assert(hMCGen_ != 0);

	h_.MCGenReco = hMCGenReco;
	nVtx_ = h_.MCGenReco->GetNbinsY()-1; // minus one, since the number of bins is +1...

	hTarget_ = new TH1D(UID+"nRecoTarget","Number of Target reco vertices; nRecoTarget; A.U.",
			nVtx_+1,-0.5,nVtx_+0.5);
	toDelete_.push_back(hTarget_);

	for(int iVtx=0; iVtx<=nVtx_; ++iVtx){
		hTarget_->SetBinContent(iVtx+1, TMath::PoissonI(iVtx, RecoPoissonMean));
	}

	makeWeights();
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::ReweightMC(TH1* hMCGen, float GenPoissonMean):
	target_(Gen),
	hMCGenOrGenReco_(hMCGen),
	hMCGen_(0)
{
	assert(hMCGenOrGenReco_ != 0);

	h_.MCGen = hMCGen;
	nVtx_ = h_.MCGen->GetNbinsX()-1; // minus one, since the number of bins is +1...

	hTarget_ = new TH1D(UID+"nGenTarget","Number of Target gen vertices; nGenTarget; A.U.",
			nVtx_+1,-0.5,nVtx_+0.5);
	toDelete_.push_back(hTarget_);

	for(int iVtx=0; iVtx<=nVtx_; ++iVtx){
		hTarget_->SetBinContent(iVtx+1, TMath::PoissonI(iVtx, GenPoissonMean));
	}
	
	makeWeights();
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::ReweightMC(TH1* hMCGen, TH1* hGenTarget):
	hMCGenOrGenReco_(hMCGen),
	hMCGen_(0),
	hTarget_(hGenTarget)

{
	assert(hMCGenOrGenReco_ != 0);
	assert(hTarget_ != 0);

	makeWeights();
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::ReweightMC(TH2 *hMCGenReco, TH1 *hMCGen, TH1 *hRecoTarget):
	hMCGenOrGenReco_(hMCGenReco),
	hMCGen_(hMCGen),
	hTarget_(hRecoTarget)
{
	assert(hMCGenOrGenReco_ != 0);
	assert(hMCGen_ != 0);
	assert(hTarget_ != 0);

	makeWeights();
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void ReweightMC::makeWeights()
{
	TH1::SetDefaultSumw2(true);

	// First of all make sure that all histos are the right size
	if(hMCGenOrGenReco_->IsA()->InheritsFrom(TH2::Class())){
		h_.MCGenReco= (TH2*)hMCGenOrGenReco_;
		assert(hTarget_->GetNbinsX() == h_.MCGenReco->GetNbinsY());
		assert(h_.MCGenReco->GetNbinsX() == h_.MCGenReco->GetNbinsY());
		target_=Reco;
	}else{
		h_.MCGen= hMCGenOrGenReco_;
		assert(hTarget_->GetNbinsX() == h_.MCGen->GetNbinsX());
		target_=Gen;
	}

	// set the number of vertices to the dimension
	nVtx_=hTarget_->GetNbinsX()-1; // minus one, since the number of bins is +1...

	// the weights histograms has the same dimensions as the MC weights passed
	hWeights_ = (TH1*) hMCGenOrGenReco_->Clone(UID+"MCWeights");
	toDelete_.push_back(hWeights_);
	hWeights_->Reset("ICE");
	hWeights_->SetEntries(0);
	switch(target_){
	case Reco:
		hWeights_->SetTitle("Weights to apply to MC events; nGenMC; nRecoMC");
		break;
	case Gen:
		hWeights_->SetTitle("Weights to apply to MC events; nGenMC; A.U.");
		break;
	}

	// normalize the target distribution
	assert(hTarget_->Integral() != 0);

	switch(target_){
	case Gen:
		hWeights_->Divide( hTarget_, h_.MCGen );
		break;

	case Reco:
		// makeResponse();
		TMatrixD respMatrix(nVtx_+1,nVtx_+1);
		for(int nReco=0; nReco<=nVtx_; ++nReco){
			for(int nGen=0; nGen<=nVtx_; ++nGen){
				respMatrix(nGen,nReco) =  h_.MCGenReco->GetBinContent(h_.MCGenReco->FindBin(nGen, nReco));
			}
		}
		TMatrixD inverse = respMatrix.Invert();
		for(int nReco=0; nReco<=nVtx_; ++nReco){
			for(int nGen=0; nGen<=nVtx_; ++nGen){
				if(hMCGen_->GetBinContent(hMCGen_->FindBin(nGen))==0) continue;
				hWeights_->SetBinContent(hWeights_->FindBin(nGen, nReco), inverse(nGen,nReco) * hTarget_->GetBinContent(hTarget_->FindBin(nReco)));
			}
		}
		hWeights_ = ((TH2*)hWeights_)->ProjectionX();
		hWeights_->Divide(hMCGen_);

		//// for(int nReco=0; nReco<=nVtx_; ++nReco){
		//// 	for(int nGen=0; nGen<=nVtx_; ++nGen){
		//// 		if(hMCGen_->GetBinContent(hMCGen_->FindBin(nGen))==0) continue;
		//// 
		//// 		hWeights_->SetBinContent(hWeights_->FindBin(nGen, nReco),
		//// 				hTarget_->GetBinContent(hTarget_->FindBin(nReco))
		//// 				* h_.MCGenReco->GetBinContent(h_.MCGenReco->FindBin(nGen, nReco))
		//// 				/ hMCGen_->GetBinContent(hMCGen_->FindBin(nGen))
		//// 				);
		//// 	}
		//// }
		break;
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
void ReweightMC::makeResponse(){

	switch(target_){
	case Gen:
		assert(h_.MCGen->Integral() != 0);
		h_.MCGen->Scale( 1.0 / h_.MCGen->Integral() );
		hMCGen_ = h_.MCGen;
		break;

	case Reco:
		hMCGen_ = (TH1 *) h_.MCGenReco->ProjectionX()->Clone(UID+"MCGen");
		hMCGen_ -> Scale( 1. / hMCGen_->Integral() );
		toDelete_.push_back(hMCGen_);

		assert(h_.MCGenReco->Integral() != 0);

		TH1 *hGen = h_.MCGenReco->ProjectionX("_Gen");
		//flatten the gen projection
		for(int nGen=0; nGen<=nVtx_; ++nGen){
			// Get the Reco distribution for each Gen value
			int binGen = h_.MCGenReco->GetXaxis()->FindBin(nGen);
			TH1 *hReco = h_.MCGenReco->ProjectionY("_Reco",binGen,binGen);
			if(hGen->GetBinContent( hGen->FindBin(nGen))==0) continue;
			hReco->Scale( 1.0 / hGen->GetBinContent( hGen->FindBin(nGen)) );

			for(int nReco=0; nReco<=nVtx_; ++nReco){
				h_.MCGenReco->SetBinContent(h_.MCGenReco->FindBin(nGen, nReco),
						hReco->GetBinContent(hReco->FindBin(nReco))
						);
			}
		}

//		TH1 *hReco = h_.MCGenReco->ProjectionY("_Reco");
//		// flatten the reco projection
//		for(int nReco=0; nReco<=nVtx_; ++nReco){
//			// Get the Gen distribution for each Reco value
//			int binReco = h_.MCGenReco->GetYaxis()->FindBin(nReco);
//			TH1 *hGen = h_.MCGenReco->ProjectionX("_Gen",binReco,binReco);
//			if(hReco->GetBinContent( hReco->FindBin(nReco))==0) continue;
//			hGen->Scale( 1.0 / hReco->GetBinContent( hReco->FindBin(nReco)) );
//
//			for(int nGen=0; nGen<=nVtx_; ++nGen){
//				h_.MCGenReco->SetBinContent(h_.MCGenReco->FindBin(nGen, nReco),
//						hGen->GetBinContent(hGen->FindBin(nGen))
//						);
//			}
//		}

		break;
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
ReweightMC::~ReweightMC(){
	vector<TH1* >::iterator i = toDelete_.begin();
	vector<TH1* >::iterator end = toDelete_.end();
	for(; i != end; i++){
		delete(*i);
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
int ReweightMC::fillMC(int nGen, int nReco, double weight){
	assert_inrange(nGen);
	assert_inrange(nReco);
	assert(weight>=0);
	assert(target_==Reco);

	return h_.MCGenReco->Fill(nGen, nReco, weight);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
int ReweightMC::fillMC(int nGen, double weight){
	assert_inrange(nGen);
	assert(weight>=0);
	assert(target_==Gen);

	return h_.MCGen->Fill(nGen, weight);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
int ReweightMC::fillTarget(int n, double weight){
	assert_inrange(n);
	assert(weight>=0);

	return hTarget_->Fill(n, weight);
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------
double ReweightMC::getWeight(int nGenInMC, int nRecoInMC){
	assert_inrange(nGenInMC);
	assert_inrange(nRecoInMC);

	return hWeights_->GetBinContent(hWeights_->FindBin(nGenInMC, nRecoInMC));
}
