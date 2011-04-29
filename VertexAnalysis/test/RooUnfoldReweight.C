//=====================================================================-*-C++-*-
// Adapted from RooUnfold examples
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldBinByBin.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t gen)
{
//  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
//  Double_t x= gRandom->Rndm();
//  if (x>xeff) return cutdummy;
  Double_t smear= gRandom->Gaus( -0.5-(gen/5.), 0.5+(gen/20.) );     // bias and smear
  Double_t reco = gen+ smear;
  if (reco<1) return cutdummy;
  return reco;
}

TH1* Unfold(TH2*hgen_vs_reco, TH1*hgen, TH1*hdata){
	  RooUnfoldResponse response (hgen_vs_reco->ProjectionX(), hgen, hgen_vs_reco);
	  RooUnfoldBayes    unfold (&response, hdata, 4);
	  return (TH1D*) unfold.Hreco(3)->Clone("hunfolded");
}

/*

.L ../../../h2gglobe/VertexAnalysis/test/RooUnfoldReweight.C
TFile b("/afs/cern.ch/user/m/malberti/public/xLino/NDataRecoVertices.root");
TFile a("/afs/cern.ch/user/m/malberti/public/xLino/HistosForPUReweighting.root");
TH1D *data  = (TH1D*) b.Get("NvtAll");
TH1D *gen  = (TH1D*) a.Get("hgen");
TH2D *gen_vs_reco  = (TH2D*) a.Get("hnpugen_vs_nvtxreco");
TFile f("unfolded.root","RECREATE");
TH1D *unf = Unfold(gen_vs_reco, gen, data);
unf->Write();
f.Close();

*/


//==============================================================================
// Example Unfolding
//==============================================================================

void RooUnfoldReweight()
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif

  cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse MCresponse (30, 1, 30);
  TH2D *hResponse2D= new TH2D("Response2D","Response 2D; reco; gen",30,1,30,30,1,30);
  TH1D *hResponseGen= new TH1D("ResponseGen","Response Gen; gen",30,1,30);

  // Fill response
  for (Int_t i= 0; i<110000; i++) {
    Double_t gen= 0;

    if(i<100000){
    	gen= gRandom->Uniform(0, 16);
    }else{
        gen= 16+10*gRandom->Exp(0.2);
    }
    Double_t reco= smear (gen);
    if (reco!=cutdummy){
      MCresponse.Fill (reco, gen);
      hResponse2D->Fill(reco,gen);
      hResponseGen->Fill(gen);
    }else{
      MCresponse.Miss (gen);
      hResponseGen->Fill(gen);
    }
  }

  RooUnfoldResponse MChistoResponse (hResponse2D->ProjectionX(), hResponseGen, hResponse2D);

  cout << "==================================== TEST =====================================" << endl;
  TH1D* hGen= new TH1D ("Gen", "Generated; gen",    30, 1, 30);
  TH1D* hReco= new TH1D ("Reco", "Reconstructed; reco", 30, 1, 30);
  // Test with a Poisson
  for (Int_t i=0; i<1000; i++) {
    Double_t gen= 1+gRandom->Poisson(15);
	Double_t reco= smear (gen);
    hGen->Fill(gen);
    if (reco!=cutdummy) hReco->Fill(reco);
  }

  cout << "==================================== UNFOLD ===================================" << endl;

  RooUnfoldBayes    unfold (&MChistoResponse, hReco, 4);    // OR
//  RooUnfoldBayes    unfold (&MCresponse, hReco, 4);    // OR

  //  RooUnfoldSvd      unfold (&MCresponse, hReco, 6);   // OR
//  RooUnfoldBinByBin unfold (&MCresponse, hReco);
//  RooUnfoldTUnfold	unfold (&MCresponse, hReco);
//  RooUnfoldInvert	unfold (&MCresponse, hReco);

  TH1D* hUnfolded= (TH1D*) unfold.Hreco(2);

  unfold.PrintTable (cout, hGen);

  TCanvas *c = new TCanvas("unfold test", "unfold test", 1200,600);
  c->Divide(3,2);

  c->cd(1);
  TH2* r = (TH2*) MCresponse.Hresponse()->Clone("Response");
  r->SetTitle("Response; reco; gen");
  r->Draw("box");

  c->cd(2);
  TH1 *hMCReco = (TH1*)MCresponse.Hmeasured()->Clone("Reco");
  hMCReco->SetTitle("Reconstructed");
  hMCReco->SetLineColor(8);
  hMCReco->SetFillColor(8);
  hMCReco->SetFillStyle(3003);
  hMCReco->Draw();

  TH1 *hMCGen = (TH1*)MCresponse.Htruth()->Clone("Gen");
  hMCGen->SetTitle("Generated");
  hMCGen->Draw("same");

  leg2 = new TLegend(0.7,0.7,0.95,0.95);
  leg2->SetHeader("MC response");
  leg2->AddEntry(hMCReco,"Reco","l");
  leg2->AddEntry(hMCGen,"Gen","l");
  leg2->Draw();

  c->cd(3);
  hReco->SetLineColor(8);
  hReco->SetFillColor(8);
  hReco->SetFillStyle(3003);
  hReco->Draw("");

  hGen->Draw("SAME");

  hUnfolded->Draw("same");

  leg1 = new TLegend(0.7,0.7,0.95,0.95);
  leg1->SetHeader("Unfolding");
  leg1->AddEntry(hReco,"Reco","f");
  leg1->AddEntry(hGen,"Gen","l");
  leg1->AddEntry(hUnfolded,"Unfolded","ep");
  leg1->Draw();

  c->cd(4);
  TH1* hWeights = (TH1*) hUnfolded->Clone("Weights");
  hWeights->SetTitle("Weights;;weight");
  hWeights->Divide(hUnfolded, hMCGen, 1./hUnfolded->Integral(), 1./hMCGen->Integral() );
  hWeights->Draw("e");

  c->cd(5);
  TH1D* hSampled= new TH1D ("Sampled", "Sampled; reco",    30, 1, 30);
  TH2D* hPull= new TH2D ("Pulls", "Pulls; pull; ",    30, 1, 30, 100, -10, 10);

  for(Int_t toy=0; toy<100; ++toy){
	  hSampled->Reset("ICE");
	  hSampled->SetEntries(0);

	  Double_t oversampling = 10.;
	  for (Int_t i=0; i<1000*oversampling; ++i) {
		  Double_t gen= hMCGen->GetRandom();
		  Double_t w= gRandom->Gaus(
				  hWeights->GetBinContent(hWeights->FindBin(gen)),
				  hWeights->GetBinError(hWeights->FindBin(gen))
				  );
		  Double_t reco=smear(gen);
		  hSampled->Fill(reco, w/oversampling);
	  }

	  for(Int_t bin=1; bin<=hSampled->GetNbinsX(); ++bin){
		  Double_t data= hSampled->GetBinContent(bin);
		  Double_t err= hSampled->GetBinError(bin);
		  Double_t ref= hReco->GetBinContent(bin);
		  if(err==0) continue;

		  hPull->Fill(hSampled->GetBinCenter(bin), (data-ref)/err);
	  }
  }

  hSampled->Draw("e");
  hReco->Draw("same");

  leg5 = new TLegend(0.7,0.7,0.95,0.95);
  leg5->SetHeader("Test weights");
  leg5->AddEntry(hReco,"Reco","f");
  leg5->AddEntry(hSampled,"Sampled","ep");
  leg5->Draw();

  TH1D* hMean= new TH1D ("Mean", "Mean;;mean", 30, 1, 30);
  TH1D* hSigma= new TH1D ("Sigma", "Sigma;;sigma", 30, 1, 30);
  TF1 *norm= new TF1("normal","gaus",-10,10);

  for(Int_t bin=0; bin<=hPull->GetNbinsX()+1; ++bin){
	  TH1D *hSlice = hPull->ProjectionY("PullSlice",bin,bin,"e");
	  if (hSlice->Integral() == 0) continue;

	  norm->SetParameter(0, hSlice->Integral());
	  norm->SetParameter(1, 0.);
	  norm->SetParameter(2, 1.);

	  hSlice->Fit(norm, "MINE");

	  hMean->SetBinContent(hPull->GetBinCenter(bin), norm->GetParameter(1) );
	  hMean->SetBinError(hPull->GetBinCenter(bin), norm->GetParError(1) );

	  hSigma->SetBinContent(hPull->GetBinCenter(bin), norm->GetParameter(2) );
	  hSigma->SetBinError(hPull->GetBinCenter(bin), norm->GetParError(2) );
  }


  c->cd(6);
//  TH1* hRatio = (TH1*) hSampled->Clone("Ratio");
//  hRatio->SetTitle("Ratio;;sampled/reco");
//  hRatio->Divide(hReco);
//  hRatio->Fit("pol0");
//  hRatio->Draw("e");

  hMean->SetMarkerColor(8);
  hMean->Draw();
  hMean->SetMarkerStyle(24);
  hSigma->Draw("same");

  leg6 = new TLegend(0.7,0.7,0.95,0.95);
  leg6->SetHeader("Pulls gaus fit");
  leg6->AddEntry(hMean,"Mean","ep");
  leg6->AddEntry(hSigma,"Sigma","ep");
  leg6->Draw();


}

#ifndef __CINT__
int main () { RooUnfoldReweight(); return 0; }  // Main program when run stand-alone
#endif
