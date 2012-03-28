#include "../interface/MvaAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0 

using namespace std;

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::MvaAnalysis()  : 
    name_("MvaAnalysis"),
    vtxAna_(vtxAlgoParams), vtxConv_(vtxAlgoParams)
{

    systRange  = 3.; // in units of sigma
    nSystSteps = 1;    
    nMasses  = 9;
}

// ----------------------------------------------------------------------------------------------------
MvaAnalysis::~MvaAnalysis() 
{
    //Default
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Term(LoopAll& l) 
{
    eventListText.close();
    if (doTraining){
        for (int i = 0; i<nMasses;i++){
            cout << " ----------- Writing trees ------------ " << endl;
            mvaFile_->cd();
            if (0<signalTrainTree_[i]->GetEntries()    ){
                signalTrainTree_[i]->Write(("sig_train"+names[i]).c_str());
                cout << "Wrote " << signalTrainTree_[i]->GetName() << endl;
            }
            if (0<signalTestTree_[i]->GetEntries()    ){
                signalTestTree_[i]->Write(("sig_test"+names[i]).c_str());
                cout << "Wrote " << signalTestTree_[i]->GetName() << endl;
            }
            if (0<backgroundTrainTree_7pt_[i]->GetEntries()    ){
                backgroundTrainTree_7pt_[i]->Write(("bkg_train_7pt"+names[i]).c_str());
                cout << "Wrote " << backgroundTrainTree_7pt_[i]->GetName() << endl;
            }
            if (0<backgroundTestTree_7pt_[i]->GetEntries()    ){
                backgroundTestTree_7pt_[i]->Write(("bkg_test_7pt"+names[i]).c_str());
                cout << "Wrote " << backgroundTestTree_7pt_[i]->GetName() << endl;
            }
            if (0<backgroundTrainTree_2pt_[i]->GetEntries()    ){
                backgroundTrainTree_2pt_[i]->Write(("bkg_train_2pt"+names[i]).c_str());
                cout << "Wrote " << backgroundTrainTree_2pt_[i]->GetName() << endl;
            }
            if (0<backgroundTestTree_2pt_[i]->GetEntries()    ){
                backgroundTestTree_2pt_[i]->Write(("bkg_test_2pt"+names[i]).c_str());
                cout << "Wrote " << backgroundTestTree_2pt_[i]->GetName() << endl;
            }
        }
        mvaFile_->Close();
    }
    else{

        // -----------------------------
        // l.rooContainer->AddRealVar("r1",-8,-15.,0.);
        // l.rooContainer->AddRealVar("r2",-0.05,-15.,0.);
        // l.rooContainer->AddRealVar("f2",0.01,0.,1.);

        // Power law
        // std::vector<std::string> data_pow_pars(3,"p");   
        // data_pow_pars[0] = "r1";
        /// data_pow_pars[1] = "r2";
        //data_pow_pars[2] = "f2";
        //l.rooContainer->AddGenericPdf("data_pow_model", "(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)","CMS_hgg_mass",data_pow_pars,0);
        /*// 5th Order Polynomial
          l.rooContainer->AddRealVar("pol0",-0.5,-1.5,1.5);
          l.rooContainer->AddRealVar("pol1",0.4,-1.5,1.5);
          l.rooContainer->AddRealVar("pol2",0.04,-1.5,1.5);
          l.rooContainer->AddRealVar("pol3",-0.03,-1.5,1.5);
          l.rooContainer->AddRealVar("pol4",-0.03,-1.5,1.5);
          l.rooContainer->AddFormulaVar("modpol0","@0*@0","pol0");
          l.rooContainer->AddFormulaVar("modpol1","@0*@0","pol1");
          l.rooContainer->AddFormulaVar("modpol2","@0*@0","pol2");
          l.rooContainer->AddFormulaVar("modpol3","@0*@0","pol3");
          l.rooContainer->AddFormulaVar("modpol4","@0*@0","pol4");

          std::vector<std::string> data_pol_pars(5,"p");     
          data_pol_pars[0] = "modpol0";
          data_pol_pars[1] = "modpol1";
          data_pol_pars[2] = "modpol2";
          data_pol_pars[3] = "modpol3";
          data_pol_pars[4] = "modpol4";
          l.rooContainer->AddGenericPdf("data_pow_model",
          "0","CMS_hgg_mass",data_pol_pars,75); // >= 71 means RooBernstein of order >= 1
          // -----------------------------
          */

        // loop hypothesis  
        for (double mass=110.0; mass<150.2; mass+=0.5){
            if (mass < masses[1] || mass > masses[nMasses-1] ) continue;  

            // define hypothesis masses for the sidebands
            float mass_hypothesis = mass;

            // define the signal Region
            float sideband_boundaries[2];

            sideband_boundaries[0] = mass_hypothesis*(1.-sidebandWidth);
            sideband_boundaries[1] = mass_hypothesis*(1.+sidebandWidth);

            // -----------------------------
            l.rooContainer->AddRealVar(Form("r1_%3.1f",mass),-8.,-10.,0.);
            l.rooContainer->AddRealVar(Form("r2_%3.1f",mass),-0.05,-10.,0.);
            l.rooContainer->AddRealVar(Form("f2_%3.1f",mass),0.01,0.,1.);
            /*
            // 5th Order Polynomial
            l.rooContainer->AddRealVar(Form("pol0_%3.1f",mass),-0.05,-1.5,1.5);
            l.rooContainer->AddRealVar(Form("pol1_%3.1f",mass),0.05,-1.5,1.5);
            l.rooContainer->AddRealVar(Form("pol2_%3.1f",mass),0.05,-1.5,1.5);
            l.rooContainer->AddRealVar(Form("pol3_%3.1f",mass),-0.001,-1.5,1.5);
            l.rooContainer->AddRealVar(Form("pol4_%3.1f",mass),-0.001,-1.5,1.5);
            l.rooContainer->AddFormulaVar(Form("modpol0_%3.1f",mass),"@0*@0",Form("pol0_%3.1f",mass));
            l.rooContainer->AddFormulaVar(Form("modpol1_%3.1f",mass),"@0*@0",Form("pol1_%3.1f",mass));
            l.rooContainer->AddFormulaVar(Form("modpol2_%3.1f",mass),"@0*@0",Form("pol2_%3.1f",mass));
            l.rooContainer->AddFormulaVar(Form("modpol3_%3.1f",mass),"@0*@0",Form("pol3_%3.1f",mass));
            l.rooContainer->AddFormulaVar(Form("modpol4_%3.1f",mass),"@0*@0",Form("pol4_%3.1f",mass));

            std::vector<std::string> data_poly4_pars(4,"p");     
            data_poly4_pars[0] = Form("modpol0_%3.1f",mass);
            data_poly4_pars[1] = Form("modpol1_%3.1f",mass);
            data_poly4_pars[2] = Form("modpol2_%3.1f",mass);
            data_poly4_pars[3] = Form("modpol3_%3.1f",mass);

            l.rooContainer->AddGenericPdf(Form("data_poly4_model_%3.1f",mass),
            "0","CMS_hgg_mass",data_poly4_pars,64); // >= 71 means RooBernstein of order >= 1

            std::vector<std::string> data_pol_pars(5,"p");   
            data_pol_pars[0] = Form("modpol0_%3.1f",mass);
            data_pol_pars[1] = Form("modpol1_%3.1f",mass);
            data_pol_pars[2] = Form("modpol2_%3.1f",mass);
            data_pol_pars[3] = Form("modpol3_%3.1f",mass);
            data_pol_pars[4] = Form("modpol4_%3.1f",mass);

            l.rooContainer->AddGenericPdf(Form("data_pow_model_%3.1f",mass),
            "0","CMS_hgg_mass",data_pol_pars,65);   // >= 71 means RooBernstein of order >= 1
            // -----------------------------
            */

            // Power law
            std::vector<std::string> data_pow_pars(3,Form("p_%3.1f",mass));   
            data_pow_pars[0] = Form("r1_%3.1f",mass);
            data_pow_pars[1] = Form("r2_%3.1f",mass);
            data_pow_pars[2] = Form("f2_%3.1f",mass);
            l.rooContainer->AddGenericPdf(Form("data_pow_model_%3.1f",mass), "(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)","CMS_hgg_mass",data_pow_pars,0);

            //      l.rooContainer->FitToData(Form("data_poly4_model_%3.1f",mass), "data_mass",massMin,sideband_boundaries[0],sideband_boundaries[1],massMax); // try to fit with lower order first
            l.rooContainer->FitToData(Form("data_pow_model_%3.1f",mass), "data_mass",massMin,sideband_boundaries[0],sideband_boundaries[1],massMax);
            std::vector<std::pair<double,double> > N_sigErr = l.rooContainer->GetFitNormalisationsAndErrors(Form("data_pow_model_%3.1f",mass),"data_mass",sideband_boundaries[0],sideband_boundaries[1]);
            //l.rooContainer->FitToData("data_pow_model", "data_mass",massMin,sideband_boundaries[0],sideband_boundaries[1],massMax);
            //std::vector<std::pair<double,double> > N_sigErr = l.rooContainer->GetFitNormalisationsAndErrors("data_pow_model","data_mass",sideband_boundaries[0],sideband_boundaries[1]);

            //l.rooContainer->AddNormalisationSystematics(Form("bkg_norm_%3.1f",mass),N_sigErr, 1); // 1 means it effect the background only

            // Calculate weights to apply to the sidebands
            std::vector<double> N_sig = l.rooContainer->GetFitNormalisations(Form("data_pow_model_%3.1f",mass),"data_mass",sideband_boundaries[0],sideband_boundaries[1]);
            //std::vector<double> N_sig = l.rooContainer->GetFitNormalisations("data_pow_model","data_mass",sideband_boundaries[0],sideband_boundaries[1]);
            // Used to store this number in the histograms of original sideband model, now just store as a RooRealVar
            l.rooContainer->AddConstant(Form("NBkgInSignal_mH%3.1f",mass),N_sig[0]);
            std::cout << "Number Of Events In Background -- " << N_sig[0] << std::endl;

            bool scale = true;//true;
            std::vector<string> ada_bkgsets;
            std::vector<string> grad_bkgsets;
            std::vector<string> vbf_bkgsets;
            //std::vector<string> ada_datasets;
            //std::vector<string> grad_datasets;

            // For summing over sidebands, can skip some of them (defined in .dat) and only include those iside the range 100->180
            // Dont use all of the sidebands availablbe but rather only up to the number of sidebands set by numberOfSidebandsForAlgos
            // Actually, we only *need* these sums for 5GeV steps

            if ((int)mass - mass==0 && ((int)mass)%5==0){
                for (int sideband_i=numberOfSidebandGaps+1;sideband_i<=numberOfSidebandsForAlgos+numberOfSidebandGaps;sideband_i++) {

                    // Calculate mass hypothesis of sideband and check its boudaries
                    double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                    double mass_hypothesis_sb = (mass_hypothesis*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));               
                    double sideband_boundaries_high= mass_hypothesis_sb*(1.+sidebandWidth);

                    if (sideband_boundaries_high <= massSidebandMax){

                        ada_bkgsets.push_back(Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass));
                        //   ada_datasets.push_back(Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass));
                        grad_bkgsets.push_back(Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass));
                        vbf_bkgsets.push_back(Form("bkg_%dhigh_VBF_%3.1f",sideband_i,mass));
                        //   grad_datasets.push_back(Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass));
                    }

                    hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                    mass_hypothesis_sb = (mass_hypothesis*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                    double sideband_boundaries_low = mass_hypothesis_sb*(1.-sidebandWidth);

                    if (sideband_boundaries_low >= massSidebandMin){
                        ada_bkgsets.push_back(Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass));
                        //           ada_datasets.push_back(Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass));
                        grad_bkgsets.push_back(Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass));
                        vbf_bkgsets.push_back(Form("bkg_%dlow_VBF_%3.1f",sideband_i,mass));
                        //           grad_datasets.push_back(Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass));
                    }

                }

                // Also for the bkg MC we want to throw in the signal region for the binning algos
                ada_bkgsets.push_back(Form("bkg_BDT_ada_%3.1f",mass)) ;
                grad_bkgsets.push_back(Form("bkg_BDT_grad_%3.1f",mass)) ;
                vbf_bkgsets.push_back(Form("bkg_VBF_%3.1f",mass)) ;


                l.rooContainer->SumMultiBinnedDatasets(Form("bkg_BDT_ada_all_%3.1f",mass),ada_bkgsets, N_sig, scale);
                l.rooContainer->SumMultiBinnedDatasets(Form("bkg_BDT_grad_all_%3.1f",mass),grad_bkgsets, N_sig, scale);
                l.rooContainer->SumMultiBinnedDatasets(Form("bkg_VBF_all_%3.1f",mass),vbf_bkgsets, -1, true); // N<0, implies no rescaling when true flag is used

                //      l.rooContainer->SumMultiBinnedDatasets(Form("data_BDT_sideband_ada_%3.1f",mass),ada_datasets, N_sig, scale);
                //      l.rooContainer->SumMultiBinnedDatasets(Form("data_BDT_sideband_grad_%3.1f",mass),grad_datasets, N_sig, scale);
            }

        }
        for (int i=1; i<nMasses; i++){
            // Need to binning on MC!
            //        std::vector <std::vector<double> > optimizedGradBins =  l.rooContainer->OptimizedBinning("bkg_BDT_grad_all"+names[i],50,false,true,-1);
            //        std::vector<std::vector <double> > optimizedAdaBins =  l.rooContainer->OptimizedBinning("bkg_BDT_ada_all"+names[i],50,false,true,-1);
            //        std::vector <std::vector<double> > optimizedGradBins =  l.rooContainer->RebinConstantEdges("bkg_BDT_grad_all"+names[i],100);
            //        std::vector<std::vector <double> > optimizedAdaBins =  l.rooContainer->RebinConstantEdges("bkg_BDT_ada_all"+names[i],100);

            // use ggh as signal modell , should ideall sum them first

            std::vector<std::string> grad_sigsets;
            std::vector<std::string> ada_sigsets;
            std::vector<std::string> vbf_sigsets;
            ada_sigsets.push_back("sig_BDT_ada_ggh"+names[i]);
            ada_sigsets.push_back("sig_BDT_ada_vbf"+names[i]);
            ada_sigsets.push_back("sig_BDT_ada_wzh"+names[i]);
            ada_sigsets.push_back("sig_BDT_ada_tth"+names[i]);
            grad_sigsets.push_back("sig_BDT_grad_ggh"+names[i]);
            grad_sigsets.push_back("sig_BDT_grad_vbf"+names[i]);
            grad_sigsets.push_back("sig_BDT_grad_wzh"+names[i]);
            grad_sigsets.push_back("sig_BDT_grad_tth"+names[i]);
            vbf_sigsets.push_back("sig_VBF_ggh"+names[i]);
            vbf_sigsets.push_back("sig_VBF_vbf"+names[i]);
            vbf_sigsets.push_back("sig_VBF_wzh"+names[i]);
            vbf_sigsets.push_back("sig_VBF_tth"+names[i]);

            l.rooContainer->SumMultiBinnedDatasets("sig_BDT_ada_all"+names[i],ada_sigsets,-1.,true);
            l.rooContainer->SumMultiBinnedDatasets("sig_BDT_grad_all"+names[i],grad_sigsets,-1.,true);
            l.rooContainer->SumMultiBinnedDatasets("sig_VBF_all"+names[i],vbf_sigsets,-1.,true);

            std::vector<double> optimizedVbfBins;
            std::vector <std::vector<double> > optimizedGradBins;
            std::vector<std::vector <double> > optimizedAdaBins; 


            // Bin edges can be rederived if flag is on --> Requires background MC
            if (rederiveOptimizedBinEdges) {
                //optimizedGradBins =  l.rooContainer->SignificanceOptimizedBinning("sig_BDT_grad_all"+names[i],"bkg_BDT_grad_all"+names[i],10);
                //optimizedAdaBins =  l.rooContainer->SignificanceOptimizedBinning("sig_BDT_ada_all"+names[i],"bkg_BDT_ada_all"+names[i],10);
                optimizedGradBins =  l.rooContainer->SoverBOptimizedBinning("sig_BDT_grad_all"+names[i],"bkg_BDT_grad_all"+names[i],20,50);
                optimizedAdaBins =  l.rooContainer->SoverBOptimizedBinning("sig_BDT_ada_all"+names[i],"bkg_BDT_ada_all"+names[i],20,50);
                //      optimizedVbfBins =  l.rooContainer->SoverBOptimizedBinning("sig_VBF_all"+names[i],"bkg_VBF_all"+names[i],10,50);
                //      optimizedVbfBins =  l.rooContainer->OptimizedBinning("bkg_VBF_all"+names[i],3,false,false-1);
                optimizedVbfBins.push_back(1.);
                //  optimizedVbfBins.push_back(1.013333333333333);
                //  optimizedVbfBins.push_back(1.026666666666666);
                optimizedVbfBins.push_back(1.04);
            } else {
                if (i==1){
                    optimizedVbfBins=VbfBinEdges_110;
                    optimizedGradBins.push_back(GradBinEdges_110);
                    optimizedAdaBins.push_back(AdaBinEdges_110);
                }
                else if (i==2){
                    optimizedVbfBins=VbfBinEdges_115;
                    optimizedGradBins.push_back(GradBinEdges_115);
                    optimizedAdaBins.push_back(AdaBinEdges_115);
                }
                else if (i==3){
                    optimizedVbfBins=VbfBinEdges_120;
                    optimizedGradBins.push_back(GradBinEdges_120);
                    optimizedAdaBins.push_back(AdaBinEdges_120);
                }
                else if (i==4){
                    optimizedVbfBins=VbfBinEdges_125;
                    optimizedGradBins.push_back(GradBinEdges_125);
                    optimizedAdaBins.push_back(AdaBinEdges_125);
                }
                else if (i==5){
                    optimizedVbfBins=VbfBinEdges_130;
                    optimizedGradBins.push_back(GradBinEdges_130);
                    optimizedAdaBins.push_back(AdaBinEdges_130);
                }
                else if (i==6){
                    optimizedVbfBins=VbfBinEdges_135;
                    optimizedGradBins.push_back(GradBinEdges_135);
                    optimizedAdaBins.push_back(AdaBinEdges_135);
                }
                else if (i==7){
                    optimizedVbfBins=VbfBinEdges_140;
                    optimizedGradBins.push_back(GradBinEdges_140);
                    optimizedAdaBins.push_back(AdaBinEdges_140);
                }
                else if (i==8){
                    optimizedVbfBins=VbfBinEdges_150;
                    optimizedGradBins.push_back(GradBinEdges_150);
                    optimizedAdaBins.push_back(AdaBinEdges_150);
                }

            }


            double mass_h_low;      
            double mass_h_high;

            if (i==7){ // the 140 case
                mass_h_low  =-2.5;
                mass_h_high =4.6;
            } else if (i==8){  // the 150 case
                mass_h_low  =-5.0;
                mass_h_high =0.1;
            } else {
                mass_h_low  =-2.5;
                mass_h_high =2.1;
            }

            for (double mass=mass_h_low; mass<mass_h_high; mass+=0.5){

                double mass_hypothesis = masses[i]+mass;
                if (mass_hypothesis < 110 || mass_hypothesis > 150) continue;

                //l.rooContainer->RebinBinnedDataset(Form("bkg_grad_%3.1f",mass_hypothesis),Form("data_BDT_sideband_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);
                l.rooContainer->RebinBinnedDataset(Form("bkg_mc_grad_%3.1f",mass_hypothesis),Form("bkg_BDT_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);
                // l.rooContainer->RebinBinnedDataset(Form("zee_mc_grad_%3.1f",mass_hypothesis),Form("zee_BDT_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);
                //l.rooContainer->RebinBinnedDataset(Form("bkg_mc_balg_grad_%3.1f",mass_hypothesis),Form("bkg_BDT_grad_all_%3.1f",mass_hypothesis),optimizedGradBins,false); 
                l.rooContainer->RebinBinnedDataset(Form("data_grad_%3.1f",mass_hypothesis),Form("data_BDT_grad_%3.1f",mass_hypothesis),optimizedGradBins,false);

                //l.rooContainer->RebinBinnedDataset(Form("bkg_ada_%3.1f",mass_hypothesis),Form("data_BDT_sideband_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);
                l.rooContainer->RebinBinnedDataset(Form("bkg_mc_ada_%3.1f",mass_hypothesis),Form("bkg_BDT_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);
                // l.rooContainer->RebinBinnedDataset(Form("zee_mc_ada_%3.1f",mass_hypothesis),Form("zee_BDT_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);
                // l.rooContainer->RebinBinnedDataset(Form("bkg_mc_balg_ada_%3.1f",mass_hypothesis),Form("bkg_BDT_ada_all_%3.1f",mass_hypothesis),optimizedAdaBins,false);
                l.rooContainer->RebinBinnedDataset(Form("data_ada_%3.1f",mass_hypothesis),Form("data_BDT_ada_%3.1f",mass_hypothesis),optimizedAdaBins,false);

                l.rooContainer->RebinBinnedDataset(Form("bkg_mc_vbf_%3.1f",mass_hypothesis),Form("bkg_VBF_%3.1f",mass_hypothesis),optimizedVbfBins,false);
                l.rooContainer->RebinBinnedDataset(Form("data_vbf_%3.1f",mass_hypothesis),Form("data_VBF_%3.1f",mass_hypothesis),optimizedVbfBins,false);
    
                if (includeVBF){
                    // After Rebinning, now we want to include the VBF dataset
                    l.rooContainer->MergeHistograms(Form("bkg_mc_grad_%3.1f",mass_hypothesis),Form("bkg_mc_vbf_%3.1f",mass_hypothesis));
                    l.rooContainer->MergeHistograms(Form("data_grad_%3.1f",mass_hypothesis),Form("data_vbf_%3.1f",mass_hypothesis));
                    l.rooContainer->MergeHistograms(Form("bkg_mc_ada_%3.1f",mass_hypothesis),Form("bkg_mc_vbf_%3.1f",mass_hypothesis));
                    l.rooContainer->MergeHistograms(Form("data_ada_%3.1f",mass_hypothesis),Form("data_vbf_%3.1f",mass_hypothesis));
                }
                                                                                     

                for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){
                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
                    // l.rooContainer->RebinBinnedDataset(Form("zee_mc_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("zee_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);
                    // l.rooContainer->RebinBinnedDataset(Form("zee_mc_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("zee_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),optimizedGradBins,false);

                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
                    // l.rooContainer->RebinBinnedDataset(Form("zee_mc_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("zee_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);
                    // l.rooContainer->RebinBinnedDataset(Form("zee_mc_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("zee_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),optimizedAdaBins,false);

                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis),Form("data_%dhigh_VBF_%3.1f",sideband_i,mass_hypothesis),optimizedVbfBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis),Form("data_%dlow_VBF_%3.1f",sideband_i,mass_hypothesis),optimizedVbfBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_VBF_%3.1f",sideband_i,mass_hypothesis),optimizedVbfBins,false);
                    l.rooContainer->RebinBinnedDataset(Form("bkg_mc_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_VBF_%3.1f",sideband_i,mass_hypothesis),optimizedVbfBins,false);

                    if (includeVBF){
                        // Also add sidbeband from VBF
                        l.rooContainer->MergeHistograms(Form("bkg_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_mc_%dhigh_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_mc_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_mc_%dlow_grad_%3.1f",sideband_i,mass_hypothesis),Form("bkg_mc_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis));

                        l.rooContainer->MergeHistograms(Form("bkg_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_mc_%dhigh_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_mc_%dhigh_vbf_%3.1f",sideband_i,mass_hypothesis));
                        l.rooContainer->MergeHistograms(Form("bkg_mc_%dlow_ada_%3.1f",sideband_i,mass_hypothesis),Form("bkg_mc_%dlow_vbf_%3.1f",sideband_i,mass_hypothesis));
                    }
                }

            }

            for (int j=-1; j<2; j++){
                if ((i==1 && j==-1) || (i==nMasses-1 && j==1)) continue;
                l.rooContainer->RebinBinnedDataset("sig_grad_ggh"+names[i]+names[i+j],"sig_BDT_grad_ggh"+names[i+j],optimizedGradBins,true);
                l.rooContainer->RebinBinnedDataset("sig_grad_vbf"+names[i]+names[i+j],"sig_BDT_grad_vbf"+names[i+j],optimizedGradBins,true);
                l.rooContainer->RebinBinnedDataset("sig_grad_wzh"+names[i]+names[i+j],"sig_BDT_grad_wzh"+names[i+j],optimizedGradBins,true);
                l.rooContainer->RebinBinnedDataset("sig_grad_tth"+names[i]+names[i+j],"sig_BDT_grad_tth"+names[i+j],optimizedGradBins,true);
                l.rooContainer->RebinBinnedDataset("sig_ada_ggh"+names[i]+names[i+j],"sig_BDT_ada_ggh"  +names[i+j],optimizedAdaBins,true);
                l.rooContainer->RebinBinnedDataset("sig_ada_vbf"+names[i]+names[i+j],"sig_BDT_ada_vbf"  +names[i+j],optimizedAdaBins,true);
                l.rooContainer->RebinBinnedDataset("sig_ada_wzh"+names[i]+names[i+j],"sig_BDT_ada_wzh"  +names[i+j],optimizedAdaBins,true);
                l.rooContainer->RebinBinnedDataset("sig_ada_tth"+names[i]+names[i+j],"sig_BDT_ada_tth"  +names[i+j],optimizedAdaBins,true);
                l.rooContainer->RebinBinnedDataset("sig_vbf_ggh"+names[i]+names[i+j],"sig_VBF_ggh"  +names[i+j],optimizedVbfBins,true);
                l.rooContainer->RebinBinnedDataset("sig_vbf_vbf"+names[i]+names[i+j],"sig_VBF_vbf"  +names[i+j],optimizedVbfBins,true);
                l.rooContainer->RebinBinnedDataset("sig_vbf_wzh"+names[i]+names[i+j],"sig_VBF_wzh"  +names[i+j],optimizedVbfBins,true);
                l.rooContainer->RebinBinnedDataset("sig_vbf_tth"+names[i]+names[i+j],"sig_VBF_tth"  +names[i+j],optimizedVbfBins,true);

                if (includeVBF){
                    l.rooContainer->MergeHistograms("sig_grad_ggh"+names[i]+names[i+j],"sig_vbf_ggh"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_grad_vbf"+names[i]+names[i+j],"sig_vbf_vbf"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_grad_wzh"+names[i]+names[i+j],"sig_vbf_wzh"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_grad_tth"+names[i]+names[i+j],"sig_vbf_tth"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_ada_ggh"+names[i]+names[i+j],"sig_vbf_ggh"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_ada_vbf"+names[i]+names[i+j],"sig_vbf_vbf"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_ada_wzh"+names[i]+names[i+j],"sig_vbf_wzh"+names[i]+names[i+j],true);
                    l.rooContainer->MergeHistograms("sig_ada_tth"+names[i]+names[i+j],"sig_vbf_tth"+names[i]+names[i+j],true);
                }
        
            } 
        }

    }
    PhotonAnalysis::Term(l);
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Init(LoopAll& l) 
{
    if(PADEBUG) 
        cout << "InitRealMvaAnalysis START"<<endl;

    nevents=0., sumwei=0.; 
    sumaccept=0., sumsmear=0., sumev=0.;

    eventListText.open(Form("%s",l.outputTextFileName.c_str()));

    if (doTraining) {
        nMasses  = 2;
        names[0]="_121.0";
        BDTnames[0]="_121";
        masses[0] = 121.;

        names[1]="_123.0";
        BDTnames[1]="_123";
        masses[1] = 123.;
    }
    else {
        names[0]="_105.0";
        BDTnames[0]="_105";
        masses[0] = 105.;

        names[1]="_110.0";
        BDTnames[1]="_110";
        masses[1] = 110.;

        names[2]="_115.0";
        BDTnames[2]="_115";
        masses[2] = 115.;

        names[3]="_120.0";
        BDTnames[3]="_120";
        masses[3] = 120.;

        names[4]="_125.0";
        BDTnames[4]="_125";
        masses[4] = 125.;

        names[5]="_130.0";
        BDTnames[5]="_130";
        masses[5] = 130.;

        names[6]="_135.0";
        BDTnames[6]="_135";
        masses[6] = 135.;

        names[7]="_140.0";
        BDTnames[7]="_140";
        masses[7] = 140.;

        names[8]="_150.0";
        BDTnames[8]="_150";
        masses[8] = 150.;
        /*
          names[9]="_121.0";
          BDTnames[9]="_121";
          masses[9] = 121.;

          names[10]="_123.0";
          BDTnames[10]="_123";
          masses[10] = 123.;
        */
    }

    // Make sure that we wont try to use more sidebands than available
    assert(numberOfSidebandsForAlgos+numberOfSidebandGaps <= numberOfSidebands);

    //   for (float mass = 90.;
    //              mass<150.;
    //              mass*= (1+signalRegionWidth)/(1-signalRegionWidth)){
    //       bkg_masses.push_back(mass); 
    //       cout<<"background mass hypothesis: "<< mass << endl;
    //   }
    //cout << "end weird loop"<<endl;
    std::string outputfilename = (std::string) l.histFileName;
    //
    // These parameters are set in the configuration file
    std::cout
        << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "MvaAnalysis " << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << "leadEtCut "<< leadEtCut << "\n"
        << "subleadEtCut "<< subleadEtCut << "\n"
        << "doTriggerSelection "<< doTriggerSelection << "\n"
        << "nEtaCategories "<< nEtaCategories << "\n"
        << "nR9Categories "<< nR9Categories << "\n"    
        << "nPtCategories "<< nPtCategories << "\n"    
        << "doEscaleSyst "<< doEscaleSyst << "\n"
        << "doEresolSyst "<< doEresolSyst << "\n"
        << "doEcorrectionSyst "<< doEcorrectionSyst << "\n"
        << "efficiencyFile " << efficiencyFile << "\n"
        << "doPhotonIdEffSyst "<< doPhotonIdEffSyst << "\n"
        << "doR9Syst "<< doR9Syst << "\n"
        << "doVtxEffSyst "<< doVtxEffSyst << "\n"
        << "doTriggerEffSyst "<< doTriggerEffSyst << "\n"
        << "doKFactorSyst "<< doKFactorSyst << "\n"
        << "-------------------------------------------------------------------------------------- \n"
        << std::endl;

    cout << "using bdt philosophy:       " << bdtTrainingPhilosophy << endl;
    if (doTraining) cout << "Training ON:    saving trees" << endl;
    else cout << "Training OFF:     running analysis" << endl;
    if (!doTraining) cout << "Reading MVA of type:       " << MVAtype << endl; 
    // call the base class initializer
    PhotonAnalysis::Init(l);

    // Avoid reweighing from histo conainer
    for(size_t ind=0; ind<l.histoContainer.size(); ind++) {
        l.histoContainer[ind].setScale(1.);
    }

    diPhoCounter_ = l.countersred.size();
    l.countersred.resize(diPhoCounter_+1);

    // initialize the analysis variables
    nCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nCategories_ *= nR9Categories;
    if( nPtCategories != 0 ) nCategories_ *= nPtCategories;

    nPhotonCategories_ = nEtaCategories;
    if( nR9Categories != 0 ) nPhotonCategories_ *= nR9Categories;

    effSmearPars.categoryType = "2CatR9_EBEE";
    effSmearPars.n_categories = 4;
    effSmearPars.efficiency_file = efficiencyFile;

    diPhoEffSmearPars.n_categories = 8;
    diPhoEffSmearPars.efficiency_file = efficiencyFile;

    if( doEcorrectionSmear ) {
        // instance of this smearer done in PhotonAnalysis
        photonSmearers_.push_back(eCorrSmearer);
    }
    if( doEscaleSmear ) {
        photonSmearers_.push_back(eScaleSmearer);
    }
    if( doEresolSmear ) {
        // energy resolution smearing
        std::cerr << __LINE__ << std::endl; 
        eResolSmearer = new EnergySmearer( eSmearPars );
        eResolSmearer->name("E_res");
        eResolSmearer->doEnergy(false);
        eResolSmearer->scaleOrSmear(false);
        photonSmearers_.push_back(eResolSmearer);
    }
    if( doRegressionSmear ) {
        // energy regression. smearing
        std::cerr << __LINE__ << std::endl; 
        eRegressionSmearer = new EnergySmearer( eSmearPars );
        eRegressionSmearer->name("regSig");
        eRegressionSmearer->doEnergy(false);// allows for future reweighting also
        eRegressionSmearer->doRegressionSigma(true);
        photonSmearers_.push_back(eRegressionSmearer);
    }
    if( doPhotonIdEffSmear ) {
        // photon ID efficiency 
        std::cerr << __LINE__ << std::endl; 
        idEffSmearer = new EfficiencySmearer( effSmearPars );
        idEffSmearer->name("idEff");
        idEffSmearer->setEffName("ratioTP");
        idEffSmearer->init();
        idEffSmearer->doPhoId(true);
        photonSmearers_.push_back(idEffSmearer);
    }
    if( doR9Smear ) {
        // R9 re-weighting
        r9Smearer = new EfficiencySmearer( effSmearPars );
        r9Smearer->name("r9Eff");
        r9Smearer->setEffName("ratioR9");
        r9Smearer->init();
        r9Smearer->doR9(true);
        photonSmearers_.push_back(r9Smearer);
    }
    if( doVtxEffSmear ) {
        // Vertex ID
        std::cerr << __LINE__ << std::endl; 
        vtxEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );   // triplicate TF1's here
        vtxEffSmearer->name("vtxEff");
        vtxEffSmearer->setEffName("ratioVertex");
        vtxEffSmearer->doVtxEff(true);
        vtxEffSmearer->init();
        diPhotonSmearers_.push_back(vtxEffSmearer);
    }
    if( doTriggerEffSmear ) {
        // trigger efficiency
        std::cerr << __LINE__ << std::endl; 
        triggerEffSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
        triggerEffSmearer->name("triggerEff");
        triggerEffSmearer->setEffName("effL1HLT");
        triggerEffSmearer->doVtxEff(false);
        triggerEffSmearer->init();
        diPhotonSmearers_.push_back(triggerEffSmearer);
    }
    if( doPhotonMvaIdSmear ) {
        // trigger efficiency
        std::cerr << __LINE__ << std::endl; 
        photonMvaIdSmearer = new DiPhoEfficiencySmearer( diPhoEffSmearPars );
        photonMvaIdSmearer->name("phoIdMva");
        photonMvaIdSmearer->setEffName("effL1HLT");
        photonMvaIdSmearer->doVtxEff(false);
        photonMvaIdSmearer->doMvaIdEff(true);
        photonMvaIdSmearer->init();
        diPhotonSmearers_.push_back(photonMvaIdSmearer);
    }
    if(doKFactorSmear) {
        // kFactor efficiency
        std::cerr << __LINE__ << std::endl; 
        kFactorSmearer = new KFactorSmearer( kfacHist );
        kFactorSmearer->name("kFactor");
        kFactorSmearer->init();
        genLevelSmearers_.push_back(kFactorSmearer);
    }
    /*    if(doKFactorSmear2D) {
    // kFactor efficiency
    std::cerr << __LINE__ << std::endl; 
    kFactorSmearer2D = new KFactorSmearer2D( kfacHist );
    kFactorSmearer2D->name("kFactor2D");
    kFactorSmearer2D->init();
    genLevelSmearers_.push_back(kFactorSmearer2D);
    }
    */

    // Define the number of categories for the statistical analysis and
    // the systematic sets to be formed

    // FIXME move these params to config file
    l.rooContainer->SetNCategories(1);
    l.rooContainer->nsigmas = nSystSteps;
    l.rooContainer->sigmaRange = systRange;
    l.rooContainer->SaveRooDataHists(false);
    l.rooContainer->Verbose(false);

    if( doEcorrectionSmear && doEcorrectionSyst ) {
        // instance of this smearer done in PhotonAnalysis
        systPhotonSmearers_.push_back(eCorrSmearer);
        std::vector<std::string> sys(1,eCorrSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEscaleSmear && doEscaleSyst ) {
        systPhotonSmearers_.push_back( eScaleSmearer );
        std::vector<std::string> sys(1,eScaleSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doEresolSmear && doEresolSyst ) {
        systPhotonSmearers_.push_back( eResolSmearer );
        std::vector<std::string> sys(1,eResolSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doRegressionSmear && doRegressionSyst ) {
        systPhotonSmearers_.push_back( eRegressionSmearer );
        std::vector<std::string> sys(1,eRegressionSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonIdEffSmear && doPhotonIdEffSyst ) {
        systPhotonSmearers_.push_back( idEffSmearer );
        std::vector<std::string> sys(1,idEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doR9Smear && doR9Syst ) {
        systPhotonSmearers_.push_back( r9Smearer );
        std::vector<std::string> sys(1,r9Smearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doVtxEffSmear && doVtxEffSyst ) {
        systDiPhotonSmearers_.push_back( vtxEffSmearer );
        std::vector<std::string> sys(1,vtxEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doTriggerEffSmear && doTriggerEffSyst ) {
        systDiPhotonSmearers_.push_back( triggerEffSmearer );
        std::vector<std::string> sys(1,triggerEffSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if( doPhotonMvaIdSmear && doPhotonMvaIdSyst ) {
        systDiPhotonSmearers_.push_back( photonMvaIdSmearer );
        std::vector<std::string> sys(1,photonMvaIdSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    if(doKFactorSmear && doKFactorSyst) {
        systGenLevelSmearers_.push_back(kFactorSmearer);
        std::vector<std::string> sys(1,kFactorSmearer->name());
        std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
        l.rooContainer->MakeSystematicStudy(sys,sys_t);
    }
    /*    if(doKFactorSmear2D && doKFactorSyst2D) {
          systGenLevelSmearers_.push_back(kFactorSmearer2D);
          std::vector<std::string> sys(1,kFactorSmearer2D->name());
          std::vector<int> sys_t(1,-1);  // -1 for signal, 1 for background 0 for both
          l.rooContainer->MakeSystematicStudy(sys,sys_t);
          }
    */

    // ----------------------------------------------------
    // ----------------------------------------------------
    // Global systematics - Lumi
    l.rooContainer->AddGlobalSystematic("lumi",1.045,1.00);
    // ----------------------------------------------------
    l.rooContainer->AddObservable("CMS_hgg_mass",massMin,massMax);

    l.rooContainer->AddConstant("IntLumi",l.intlumi_);

    //  l.rooContainer->AddRealVar("r1",-4,-20.,0.);
    //  l.rooContainer->AddRealVar("r2",-5,-20.,0.);
    //  l.rooContainer->AddRealVar("f2",0.2,0.,1.);

    // Power law
    //std::vector<std::string> data_pow_pars(3,"p");   
    //data_pow_pars[0] = "r1";
    //data_pow_pars[1] = "r2";
    //data_pow_pars[2] = "f2";
    //l.rooContainer->AddGenericPdf("data_pow_model", "(1-@3)*TMath::Power(@0,@1) + @3*TMath::Power(@0,@2)","CMS_hgg_mass",data_pow_pars,0);

    // Initialize all MVA ---------------------------------------------------//
    l.SetAllMVA();
    // UCSD
    l.tmvaReaderID_UCSD->BookMVA("Gradient"      ,photonLevelMvaUCSD.c_str()  );
    l.tmvaReader_dipho_UCSD->BookMVA("Gradient"  ,eventLevelMvaUCSD.c_str()   );
    // MIT 
    l.tmvaReaderID_MIT_Barrel->BookMVA("AdaBoost",photonLevelMvaMIT_EB.c_str());
    l.tmvaReaderID_MIT_Endcap->BookMVA("AdaBoost",photonLevelMvaMIT_EE.c_str());
    l.tmvaReader_dipho_MIT->BookMVA("Gradient"   ,eventLevelMvaMIT.c_str()    );
    // ----------------------------------------------------------------------//

    if (doTraining){
        TString outfileName( "TMVA_input_" + (std::string) l.histFileName);
        mvaFile_ = TFile::Open( outfileName, "RECREATE" );
        mvaFile_->cd();
        for (int i = 0; i<nMasses;i++){
            cout << " --------- Creating trees --------- " << endl;
            signalTrainTree_[i] = new TTree(("sig_train"+names[i]).c_str(),("SignalTrainTree"+names[i]).c_str());
            SetBDTInputTree(signalTrainTree_[i]);
            cout << "Made " << signalTrainTree_[i]->GetName() << endl;
            signalTestTree_[i]  = new TTree(("sig_test"+names[i]).c_str(), ("SignalTestTree"+names[i]).c_str());
            SetBDTInputTree(signalTestTree_[i]);
            cout << "Made " << signalTestTree_[i]->GetName() << endl;

            backgroundTrainTree_2pt_[i] = new TTree(("bkg_train_2pt"+names[i]).c_str(),("BackgroundTrainTree_2pt"+names[i]).c_str());
            SetBDTInputTree(backgroundTrainTree_2pt_[i]);
            backgroundTestTree_2pt_[i]  = new TTree(("bkg_test_2pt"+names[i]).c_str(), ("BackgroundTestTree_2pt"+names[i]).c_str());
            SetBDTInputTree(backgroundTestTree_2pt_[i]);
            backgroundTrainTree_7pt_[i] = new TTree(("bkg_train_7pt"+names[i]).c_str(),("BackgroundTrainTree_7pt"+names[i]).c_str());
            SetBDTInputTree(backgroundTrainTree_2pt_[i]);
            backgroundTestTree_7pt_[i]  = new TTree(("bkg_test_7pt"+names[i]).c_str(), ("BackgroundTestTree_7pt"+names[i]).c_str());
            SetBDTInputTree(backgroundTestTree_2pt_[i]);
        }
    }
    else{

        l.rooContainer->AddObservable("BDT" ,-1.,1.);
        l.rooContainer->AddObservable("VBF" ,1,1+2*sidebandWidth);  // Basically a single bin just for VBF

        //Set up TMVA reader (only two variables)
        tmvaReader_= new TMVA::Reader();

        tmvaReader_->AddVariable("bdtoutput",&_bdtoutput);
        tmvaReader_->AddVariable("deltaMOverM", &_deltaMOverM);

        //Invariant Mass Spectra
        l.rooContainer->CreateDataSet("CMS_hgg_mass","data_mass",nDataBins);
        l.rooContainer->CreateDataSet("CMS_hgg_mass","bkg_mass" ,nDataBins);       
        l.rooContainer->CreateDataSet("CMS_hgg_mass","zee_mass" ,nDataBins);       

        int nBDTbins = 5000;

        // TESTING 

        for (std::vector<double>::iterator b=VbfBinEdges_110.begin();b!=VbfBinEdges_110.end();b++){

            std::cout << "Bin edge " << *b<<std::endl ;
        }
        //

        // Usual datasets
        for (double mass=110.0; mass<150.5; mass+=0.5){

            //Adaptive Boost
            l.rooContainer->CreateDataSet("BDT",Form("data_BDT_ada_%3.1f",mass)       ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("bkg_BDT_ada_%3.1f",mass)       ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("zee_BDT_ada_%3.1f",mass)       ,nBDTbins);

            //Gradient Boost
            l.rooContainer->CreateDataSet("BDT",Form("data_BDT_grad_%3.1f",mass)     ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("bkg_BDT_grad_%3.1f",mass)      ,nBDTbins);
            l.rooContainer->CreateDataSet("BDT",Form("zee_BDT_grad_%3.1f",mass)      ,nBDTbins);

            if (includeVBF){
                // VBF
                l.rooContainer->CreateDataSet("VBF",Form("data_VBF_%3.1f",mass)     ,nBDTbins);
                l.rooContainer->CreateDataSet("VBF",Form("bkg_VBF_%3.1f",mass)      ,nBDTbins);
                l.rooContainer->CreateDataSet("VBF",Form("zee_VBF_%3.1f",mass)      ,nBDTbins);
            }
      

            for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){
                // Always create all of the sidebands, even if we skip some of them
                // later on for the sums
                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass)   ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass)  ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass)  ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass) ,nBDTbins);

                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass)  ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass) ,nBDTbins);
                l.rooContainer->CreateDataSet("BDT",Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass),nBDTbins);

                if (includeVBF){
                    l.rooContainer->CreateDataSet("VBF",Form("bkg_%dlow_VBF_%3.1f",sideband_i,mass)  ,nBDTbins);
                    l.rooContainer->CreateDataSet("VBF",Form("bkg_%dhigh_VBF_%3.1f",sideband_i,mass) ,nBDTbins);
                    l.rooContainer->CreateDataSet("VBF",Form("data_%dlow_VBF_%3.1f",sideband_i,mass) ,nBDTbins);
                    l.rooContainer->CreateDataSet("VBF",Form("data_%dhigh_VBF_%3.1f",sideband_i,mass),nBDTbins);
                }

                //  l.rooContainer->CreateDataSet("BDT",Form("zee_%dlow_BDT_ada_%3.1f",sideband_i,mass)   ,nBDTbins);
                //  l.rooContainer->CreateDataSet("BDT",Form("zee_%dhigh_BDT_ada_%3.1f",sideband_i,mass)  ,nBDTbins);
                //  l.rooContainer->CreateDataSet("BDT",Form("zee_%dlow_BDT_grad_%3.1f",sideband_i,mass)   ,nBDTbins);
                //  l.rooContainer->CreateDataSet("BDT",Form("zee_%dhigh_BDT_grad_%3.1f",sideband_i,mass)  ,nBDTbins);

            }

        }

        // loop signal mass points signal datasets
        for (double sig=110.0;sig<=150;sig+=5.0){  // We are ignoring mass 105 for now
            if (sig==145.0) continue;
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_ggh_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_vbf_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_wzh_%3.1f",sig)      ,nBDTbins); 
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_grad_tth_%3.1f",sig)      ,nBDTbins); 

            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_ada_ggh_%3.1f",sig)       ,nBDTbins);    
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_ada_vbf_%3.1f",sig)       ,nBDTbins);    
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_ada_wzh_%3.1f",sig)       ,nBDTbins);    
            l.rooContainer->CreateDataSet("BDT",Form("sig_BDT_ada_tth_%3.1f",sig)       ,nBDTbins);   

            if (includeVBF){
                l.rooContainer->CreateDataSet("VBF",Form("sig_VBF_ggh_%3.1f",sig)       ,nBDTbins);    
                l.rooContainer->CreateDataSet("VBF",Form("sig_VBF_vbf_%3.1f",sig)       ,nBDTbins);    
                l.rooContainer->CreateDataSet("VBF",Form("sig_VBF_wzh_%3.1f",sig)       ,nBDTbins);    
                l.rooContainer->CreateDataSet("VBF",Form("sig_VBF_tth_%3.1f",sig)       ,nBDTbins);   
            }

            /*
              for (int sideband_i=1;sideband_i<=numberOfSidebands;sideband_i++){
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_ada_ggh_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_ada_ggh_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_ada_vbf_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_ada_vbf_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_ada_wzh_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_ada_wzh_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_ada_tth_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_ada_tth_%3.1f",sideband_i,sig)  ,nBDTbins);

              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_grad_ggh_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_grad_ggh_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_grad_vbf_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_grad_vbf_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_grad_wzh_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_grad_wzh_%3.1f",sideband_i,sig)  ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dlow_BDT_grad_tth_%3.1f",sideband_i,sig)   ,nBDTbins);
              l.rooContainer->CreateDataSet("BDT",Form("sig_%dhigh_BDT_grad_tth_%3.1f",sideband_i,sig)  ,nBDTbins);
              }
            */

            // Make the signal Systematic Sets
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_ggh_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_vbf_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_wzh_%3.1f",sig),-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_grad_tth_%3.1f",sig),-1);

            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_ada_ggh_%3.1f",sig) ,-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_ada_vbf_%3.1f",sig) ,-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_ada_wzh_%3.1f",sig) ,-1);
            l.rooContainer->MakeSystematics("BDT",Form("sig_BDT_ada_tth_%3.1f",sig) ,-1);

            if (includeVBF){
                l.rooContainer->MakeSystematics("VBF",Form("sig_VBF_ggh_%3.1f",sig) ,-1);
                l.rooContainer->MakeSystematics("VBF",Form("sig_VBF_vbf_%3.1f",sig) ,-1);
                l.rooContainer->MakeSystematics("VBF",Form("sig_VBF_wzh_%3.1f",sig) ,-1);
                l.rooContainer->MakeSystematics("VBF",Form("sig_VBF_tth_%3.1f",sig) ,-1);
            }


        }
        //TMVA Reader
        if (MVAtype=="Likelihood"){
            if (bdtTrainingPhilosophy=="MIT") {
                tmvaReader_->BookMVA("BDT_ada_123", mvaWeightsFolder+"/TMVAClassification_LikelihoodMIT.weights.xml");
                tmvaReader_->BookMVA("BDT_grad_123",mvaWeightsFolder+"/TMVAClassification_LikelihoodDMIT.weights.xml");
            }
            else {
                tmvaReader_->BookMVA("BDT_ada_123", mvaWeightsFolder+"/TMVAClassification_LikelihoodUCSD.weights.xml");
                tmvaReader_->BookMVA("BDT_grad_123",mvaWeightsFolder+"/TMVAClassification_LikelihoodDUCSD.weights.xml");
            }
        }
        else {
            if (bdtTrainingPhilosophy=="MIT") {
                tmvaReader_->BookMVA("BDT_ada_123", mvaWeightsFolder+"/TMVAClassification_BDTadaMIT.weights.xml");
                tmvaReader_->BookMVA("BDT_grad_123",mvaWeightsFolder+"/TMVAClassification_BDTgradMIT.weights.xml");
            }
            else {
                tmvaReader_->BookMVA("BDT_ada_123", mvaWeightsFolder+"/TMVAClassification_BDTadaUCSD.weights.xml");
                tmvaReader_->BookMVA("BDT_grad_123",mvaWeightsFolder+"/TMVAClassification_BDTgradUCSD.weights.xml");
            }
        }
    }

    FillSignalLabelMap();

    if(PADEBUG) 
        cout << "InitRealMvaAnalysis END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::Analysis(LoopAll& l, Int_t jentry) 
{

    if(PADEBUG) 
        cout << "Analysis START; cur_type is: " << l.itype[l.current] <<endl;
    int cur_type = l.itype[l.current];
    float weight = l.sampleContainer[l.current_sample_index].weight;

    if (!doTraining){
        if (splitSignalSample  && cur_type < 0) {
            if (jentry%2!=0) return;
            else weight*=2;
        }

        if (splitBackgroundSample && cur_type > 0) {
            if (jentry%2!=0) return;
            else weight*=2;
        }
    }


    // Get Signal Label for the current type
    std::string currentTypeSignalLabel = "";
    if  (cur_type<0) currentTypeSignalLabel = GetSignalLabel(cur_type);
    // -----------------------------------------------------------------------------------------------

    l.FillCounter( "Processed", 1. );
    assert( weight > 0. );  
    l.FillCounter( "XSWeighted", weight );
    nevents+=1.;

    // Set reRunCiC Only if this is an MC event since scaling of R9 and Energy isn't done at reduction
    if (cur_type==0) {
        l.runCiC=reRunCiCForData;
    } else {
        l.runCiC = true;
    }

    // -----------------------------------------------------------------------------------------------

    //PU reweighting
    unsigned int n_pu = l.pu_n;
    if ( cur_type !=0 && puHist != "" && cur_type < 100) {
        bool hasSpecificWeight = weights.find( cur_type ) != weights.end() ; 
        if( cur_type < 0 && !hasSpecificWeight && jentry == 1 ) {
            std::cerr  << "WARNING no pu weights specific for sample " << cur_type << std::endl;
        }
        std::vector<double> & puweights = hasSpecificWeight ? weights[ cur_type ] : weights[0]; 
        if(n_pu<puweights.size()){
            weight *= puweights[n_pu];
            sumwei+=puweights[n_pu]; 
        }    
        else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
            cout <<"n_pu ("<< n_pu<<") too big ("<<puweights.size()<<") ["<< l.itype[l.current]<<"], event will not be reweighted for pileup"<<endl;
        }
    }

    assert( weight >= 0. );  
    l.FillCounter( "PUWeighted", weight );

    if( jentry % 10000 ==  0 ) {
        std::cout << " nevents " <<  nevents << " sumpuweights " << sumwei << " ratio " << sumwei / nevents 
                  << " equiv events " << sumev << " accepted " << sumaccept << " smeared " << sumsmear << " "  
                  <<  sumaccept / sumev << " " << sumsmear / sumaccept
                  << std::endl;
    }
    // ------------------------------------------------------------
    //PT-H K-factors
    double gPT = 0;
    double gY  = 0;
    TLorentzVector gP4(0,0,0,0);
    if (cur_type<0){            // if background sample, gP4 remains 4vect(0)
        for (int gi=0;gi<l.gp_n;gi++){
            if (l.gp_pdgid[gi]==25){
                gP4 = *((TLorentzVector*)l.gp_p4->At(gi));
                gPT = gP4.Pt();
                gY  = fabs(gP4.Rapidity());
                break;
            }
        }
    }

    // ------------------------------------------------------------
    // smear all of the photons!
    std::pair<int,int> diphoton_index;

    // do gen-level dependent first (e.g. k-factor); only for signal
    double genLevWeight=1; 
    if(cur_type<0){
        for(std::vector<BaseGenLevelSmearer*>::iterator si=genLevelSmearers_.begin(); si!=genLevelSmearers_.end(); si++){
            float genWeight=1;
            (*si)->smearEvent( genWeight,gP4, l.pu_n, cur_type, 0. );
            if( genWeight < 0. ) {
                std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
                assert(0);
            }
            genLevWeight*=genWeight;
        }
    }

    // Nominal smearing
    std::vector<float> smeared_pho_energy(l.pho_n,0.); 
    std::vector<float> smeared_pho_r9(l.pho_n,0.); 
    std::vector<float> smeared_pho_weight(l.pho_n,1.);

    // TEMPORARY FIX -------------------------------------------------------------------------------------------------------//
    // Scale all the r9 of the photons in the MC
    // For now we just let it use the index but we specifically Change the r9 in the branch AFTER Energy regression smearing
    // Ideally we want to pass a smeared r9 too and apply after energy corrections, currently the smeared_pho_r9 isnt used!
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    if (cur_type !=0){
        for (int ipho=0;ipho<l.pho_n;ipho++){
            l.pho_isEB[ipho]=fabs((*((TVector3*)l.sc_xyz->At(l.pho_scind[ipho]))).Eta())<1.5;
            l.pho_r9[ipho]*=1.0035;
            if (l.pho_isEB[ipho]){ l.pho_sieie[ipho] = (0.87*l.pho_sieie[ipho]) + 0.0011 ;}
            else {l.pho_sieie[ipho]*=0.99;}
            l.sc_seta[l.pho_scind[ipho]]*=0.99;  
            l.sc_sphi[l.pho_scind[ipho]]*=0.99;  
            energyCorrectedError[ipho] *=(l.pho_isEB[ipho]) ? 1.07 : 1.045 ;
        }
    }
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    // ---------------------------------------------------------------------------------------------------------------------//
    photonInfoCollection.clear();
    for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
        std::vector<std::vector<bool> > p;
        PhotonReducedInfo phoInfo (// *((TVector3*)l.pho_calopos->At(ipho)), 
                                   *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
                                   ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
                                   energyCorrected[ipho],
                                   l.pho_isEB[ipho], l.pho_r9[ipho],
                                   l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
                                   (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
                                   );
        if (l.CheckSphericalPhoton(ipho)) phoInfo.setSphericalPhoton(true);
        float pweight = 1.;
        // smear MC. But apply energy shift to data 
        if( cur_type != 0 && doMCSmearing ) { // if it's MC
            for(std::vector<BaseSmearer *>::iterator si=photonSmearers_.begin(); si!= photonSmearers_.end(); ++si ) {
                float sweight = 1.;
                (*si)->smearPhoton(phoInfo,sweight,l.run,0.);     
                if( sweight < 0. ) {
                    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
                    assert(0);
                }
                pweight *= sweight;
            }
        } else if(cur_type == 0 ) {          // if it's data
            float sweight = 1.; 
            if (doEcorrectionSmear){
                eCorrSmearer->smearPhoton(phoInfo,sweight,l.run,0.);     // This Smearer is the same as for MC so can just re-use it
            }
            eScaleDataSmearer->smearPhoton(phoInfo,sweight,l.run,0.);
            pweight *= sweight;  
        }

        smeared_pho_energy[ipho] = phoInfo.energy();
        smeared_pho_r9[ipho] = phoInfo.r9();
        smeared_pho_weight[ipho] = pweight;
        photonInfoCollection.push_back(phoInfo);
    }

    sumev += weight;

    // Get Ready for VBF Tagging
    bool VBFevent = false;
    double leadEtCutVBF = 55.;
    // FIXME pass smeared R9
    int diphoton_id=-1;

    // VBF-TAGGING -------------------------------------------------------------------------- //
    // CP // NW Use Same pre-selected events but tag as VBF if pass Jet Requirements
    // Make sure VBF event is set to false before checking for the Jets
    VBFevent = false; 

    if(includeVBF) {
            
        // JET MET Corrections
        l.RescaleJetEnergy();

        int diphoton_id_vbf = l.DiphotonMITPreSelection(leadEtCutVBF,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 

        if (diphoton_id_vbf > -1){
            diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id_vbf],  l.dipho_subleadind[diphoton_id_vbf] );

            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
            float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id_vbf]];
            float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id_vbf]];
            TLorentzVector Higgs = lead_p4 + sublead_p4;   

            float jet1ptcut =0.0;
            float jet2ptcut =0.0;
            bool crosscheck = false;
            std::pair<int,int> highestPtJets(-1,-1);

            highestPtJets = l.Select2HighestPtJets(lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
            bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.second!=-1);

            if(VBFpresel){

                TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
                TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
                TLorentzVector dijet = (*jet1) + (*jet2);

                myVBFLeadJPt = jet1->Pt();
                myVBFSubJPt = jet2->Pt();
                myVBF_Mjj = dijet.M();
                myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
                myVBFZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
                myVBFdPhi = fabs(Higgs.DeltaPhi(dijet));
                myVBF_Mgg =Higgs.M();


                // Cannot Get Apply cuts to work -> Need to discuss with C.Palmer, for now, jyst apply the cuts
                //VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);
                VBFevent = l.ApplyCuts(0,1);
                if(VBFevent) diphoton_id = diphoton_id_vbf;

            }
        }
    }
    // CP // NW VBF Tagging
    // --------------------- END VBF-TAGGING --------------------------------------------------------//


    if (diphoton_id<0 ){ // then the VBF selection failed at some point and so we try for inclusive
        if (bdtTrainingPhilosophy=="MIT"){
            diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
        } else if (bdtTrainingPhilosophy=="UCSD"){
            diphoton_id = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
        }
    }
    /// std::cerr << "Selected pair " << l.dipho_n << " " << diphoton_id << std::endl;
    if (diphoton_id > -1 ) {

        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
        // bring all the weights together: lumi & Xsection, single gammas, pt kfactor
        float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

        l.countersred[diPhoCounter_]++;

        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
        float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
        float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
        TLorentzVector Higgs = lead_p4 + sublead_p4;   
        TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

        int category = 0;
        int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
        bool CorrectVertex;
        // FIXME pass smeared R9
        if( cur_type != 0 && doMCSmearing && cur_type < 100) {
            float pth = Higgs.Pt();
            for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
                float rewei=1.;
                (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), zero_ ,zero_,0.);
                if( rewei < 0. ) {
                    std::cerr << "Negative weight from smearer " << (*si)->name() << std::endl;
                    assert(0);
                }
                evweight *= rewei;
            }
            CorrectVertex=(*vtx- *((TVector3*)l.gv_pos->At(0))).Mag() < 1.;
        }

        // Keep some variables for plotting later (already have r9 from above);

        float mass        = Higgs.M();
        float ptHiggs     = Higgs.Pt();
        float etaHiggs    = fabs(Higgs.Eta());
        float pt_lead     = lead_p4.Pt();
        float pt_sublead  = sublead_p4.Pt();
        float eta_lead    = fabs(lead_p4.Eta());
        float eta_sublead = fabs(sublead_p4.Eta());
        float delta_phi   = fabs(lead_p4.DeltaPhi(sublead_p4));
        int histoplace    = 0;      // index for histocontaner, corresponds to sideband

        assert( evweight >= 0. ); 

        // Mass Resolution of the Event
        //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
        massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
        //    massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
        //    massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

        double massResolution = massResolutionCalculator->massResolutionCorrVtx();
        //  double massResolution = massResolutionCalculator->massResolution();  //no longer use one or other
        double vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
        float sigmaMrv = massResolutionCalculator->massResolutionEonly();
        float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
        float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
        float vtxProb = 1.-0.49*(vtx_mva+1.0);

        float bdtoutput,phoid_mvaout_lead,phoid_mvaout_sublead;
        if (bdtTrainingPhilosophy=="MIT"){
            bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"MIT");
            if (bdtoutput < 0.05) category = -1;
            phoid_mvaout_lead = l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,"MIT");
            phoid_mvaout_sublead = l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,"MIT");
        } 
        else if (bdtTrainingPhilosophy=="UCSD"){
            bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"UCSD");
            phoid_mvaout_lead = l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,"UCSD");
            phoid_mvaout_sublead = l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,"UCSD");
        }



        float bdt_grad,bdt_ada;

        if (doTraining ){//
            if (cur_type > 0 ){
                // Background:
                //
                // 2 Percent mass window N upper and N sidebands
                for (int i = 0; i<nMasses;i++) {
                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[i];
                    float mass_hypothesis_low = 0.;//mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth) -sidebandShift;
                    float mass_hypothesis_high = 0.;//mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth)+sidebandShift;

                    float sideband_boundaries[2];
                    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                    //Signal Window
                    if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight,category);
                        //std::cout<<_mgg<<", "<<mass<<std::endl;
                        _sideband = 0;
                        if (jentry%2==0) backgroundTrainTree_2pt_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_2pt_[i]->Fill();
                    }
                    // Loop over N lower sidebands
                    for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                        double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                        mass_hypothesis_low = (mass_hypothesis*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)
                            *(TMath::Power(hypothesisModifier,sideband_i-1));
                        double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
                        double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

                        if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                            SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_low,bdtoutput,evweight,category);
                            //std::cout<<_mgg<<", "<<mass<<std::endl;
                            _sideband = 1*sideband_i;
                            if (jentry%2==0) backgroundTrainTree_2pt_[i]->Fill();
                            else if (jentry%2==1) backgroundTestTree_2pt_[i]->Fill();
                        }
                    }
                    // Loop over N higher sidebands
                    for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                        double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                        mass_hypothesis_high = (mass_hypothesis*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)
                            *(TMath::Power(hypothesisModifier,sideband_i-1));
                        double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
                        double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);
                        //std::cout<<sideband_boundaries_low<<", "<<sideband_boundaries_high<<std::endl;

                        if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                            SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_high,bdtoutput,evweight,category);
                            //std::cout<<_mgg<<", "<<mass<<std::endl;
                            _sideband = -1*sideband_i;
                            if (jentry%2==0) backgroundTrainTree_2pt_[i]->Fill();
                            else if (jentry%2==1) backgroundTestTree_2pt_[i]->Fill();
                        }
                    }
                }


                // 7 Percent mass window with single upper and lower sidebands
                for (int i = 0; i<nMasses;i++) {
                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[i];
                    float mass_hypothesis_low = mass_hypothesis*(1-signalRegionWidth)/(1+signalRegionWidth);
                    float mass_hypothesis_high = mass_hypothesis*(1+signalRegionWidth)/(1-signalRegionWidth);
                    // define the sidebands
                    float sideband_boundaries[4];
                    sideband_boundaries[0] = mass_hypothesis_low*(1-signalRegionWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1-signalRegionWidth);
                    sideband_boundaries[2] = mass_hypothesis*(1+signalRegionWidth);
                    sideband_boundaries[3] = mass_hypothesis_high*(1+signalRegionWidth);

                    if( mass>sideband_boundaries[1] && mass<sideband_boundaries[2] ){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight,category);
                        _sideband = 0;
                        if (jentry%2==0) backgroundTrainTree_7pt_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_7pt_[i]->Fill();
                    }
                    else if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1] ){//lower sideband mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_low,bdtoutput,evweight,category);
                        _sideband = -1;
                        if (jentry%2==0) backgroundTrainTree_7pt_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_7pt_[i]->Fill();
                    }
                    else if( mass>sideband_boundaries[2] && mass<sideband_boundaries[3] ){//upper sideband mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_high,bdtoutput,evweight,category);
                        _sideband = 1;
                        if (jentry%2==0) backgroundTrainTree_7pt_[i]->Fill();
                        else if (jentry%2==1) backgroundTestTree_7pt_[i]->Fill();
                    }
                } 
            }
            else { //Signal 
                int i0 = SignalType(cur_type); 
                if (i0<0){
                    cout<<"CAN'T TRAIN ON DATA!\n";
                    return;
                }
                float mass_hypothesis = masses[i0];
                float sideband_boundaries[2];
                sideband_boundaries[0] = mass_hypothesis*(1-signalRegionWidth);
                sideband_boundaries[1] = mass_hypothesis*(1+signalRegionWidth);

                if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1] ){//Signal mass window cut
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,masses[i0],bdtoutput,evweight,category);
                    _sideband = 0;
                    if (jentry%2==0) signalTrainTree_[i0]->Fill();
                    else if (jentry%2==1) signalTestTree_[i0]->Fill();
                }
            }
        }
        else{

            // --- Fill invariant mass spectrum -------
            if (cur_type==0){  // Data
                l.rooContainer->InputDataPoint("data_mass",category,mass);
            } else if (cur_type>0){ // Background MC
                if (cur_type==6){l.rooContainer->InputBinnedDataPoint("zee_mass",category,mass,evweight);}
                else {l.rooContainer->InputBinnedDataPoint("bkg_mass",category,mass,evweight);}
            }

            if (bdtoutput>=0.05) l.FillHist("all_mass",0, mass, evweight);

            // ------ Deal with Signal MC first
            if (cur_type<0){ // signal MC
                // define hypothesis masses for the sidebands
                float mass_hypothesis = masses[SignalType(cur_type)];

                // define the sidebands
                float sideband_boundaries[2];
                sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                //Signal Window
                if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                    histoplace=1;
                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight);

                    bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                    bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                    if (bdtoutput>=0.05) {
                        l.FillHist("signal_pt_msig",0, Higgs.Pt(), evweight);
                        l.FillHist("signal_logpt_msig",0, log10(Higgs.Pt()), evweight);
                        l.FillHist("signal_ptOverM_msig",0, Higgs.Pt()/mass, evweight);
                        l.FillHist("signal_ptOverMH_msig",0, _H_ptOverM, evweight);
                        l.FillHist("signal_eta_msig",0, Higgs.Eta(), evweight);
                        l.FillHist("signal_deltaPhi_msig",0, _d_phi, evweight);
                        l.FillHist("signal_cosDeltaPhi_msig",0, cos(_d_phi), evweight);
                        l.FillHist("signal_deltaEta_msig",0, _d_eta, evweight);
                        l.FillHist("signal_helicityAngle_msig",0, _cos_theta_star, evweight);
                        l.FillHist("signal_pho1_pt_msig",0,lead_p4.Pt(), evweight);
                        l.FillHist("signal_pho2_pt_msig",0,sublead_p4.Pt(), evweight);
                        l.FillHist("signal_pho1_ptOverM_msig",0,lead_p4.Pt()/mass, evweight);
                        l.FillHist("signal_pho2_ptOverM_msig",0,sublead_p4.Pt()/mass, evweight);
                        l.FillHist("signal_pho1_ptOverMH_msig",0,_pho1_ptOverM, evweight);
                        l.FillHist("signal_pho2_ptOverMH_msig",0,_pho2_ptOverM, evweight);
                        l.FillHist("signal_pho1_eta_msig",0,lead_p4.Eta(), evweight);
                        l.FillHist("signal_pho2_eta_msig",0,sublead_p4.Eta(), evweight);
                        l.FillHist("signal_pho_minr9_msig",0,min(l.pho_r9[diphoton_index.first],l.pho_r9[diphoton_index.second]), evweight);
                        l.FillHist("signal_pho1_r9_msig",0,l.pho_r9[diphoton_index.first], evweight);
                        l.FillHist("signal_pho2_r9_msig",0,l.pho_r9[diphoton_index.second], evweight);
                        l.FillHist("signal_maxeta_msig",0,_max_eta, evweight);
                        l.FillHist("signal_deltaMOverMH_msig",0, _deltaMOverM, evweight);
                        l.FillHist("signal_sigmaM_msig",0,sigmaMrv, evweight);
                        l.FillHist("signal_sigmaM_wrongVtx_msig",0,sigmaMwv, evweight);
                        l.FillHist("signal_sigmaMOverM_msig",0,sigmaMrv/mass, evweight);
                        l.FillHist("signal_sigmaMOverM_wrongVtx_msig",0,sigmaMwv/mass, evweight);
                        l.FillHist("signal_sigmaMOverMH_msig",0,sigmaMrv/mass_hypothesis, evweight);
                        l.FillHist("signal_sigmaMOverMH_wrongVtx_msig",0,sigmaMwv/mass_hypothesis, evweight);
                        l.FillHist("signal_deltaMOverSigmaM_msig",0,_deltaMOverSigmaM, evweight);
                        l.FillHist("signal_deltaMOverSigmaM_wrongVtx_msig",0,(mass-mass_hypothesis)/sigmaMwv, evweight);
                        l.FillHist("signal_pho1_phoidMva_msig",0,phoid_mvaout_lead, evweight);
                        l.FillHist("signal_pho2_phoidMva_msig",0,phoid_mvaout_sublead, evweight);
                        l.FillHist("signal_vtxProb_msig",0,vtxProb, evweight);
                    }
                    l.FillHist("signal_bdtoutput_msig",0,bdtoutput, evweight);

                    if (VBFevent){ // note if includeVBF=false, VBFevent will always be false so no need to check both
                        l.rooContainer->InputBinnedDataPoint("sig_VBF_"+currentTypeSignalLabel  ,category,1.+sidebandWidth+_deltaMOverM,evweight);
                    } else {
                        l.rooContainer->InputBinnedDataPoint("sig_BDT_ada_"+currentTypeSignalLabel  ,category,bdt_ada,evweight);
                        l.rooContainer->InputBinnedDataPoint("sig_BDT_grad_"+currentTypeSignalLabel ,category,bdt_grad,evweight);
                    }
                }
                else {
                    /*  // Loop over N lower sidebands
                        for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                        double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                        double mass_hypothesis_low     = (mass_hypothesis*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                        double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
                        double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

                        if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_low,bdtoutput,evweight);
                        bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                        bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                        if (!VBFevent){ // note if includeVBF=false, VBFevent will always be false so no need to check both
                        l.rooContainer->InputBinnedDataPoint(Form("sig_%dlow_BDT_ada_",sideband_i)+currentTypeSignalLabel ,category,bdt_ada,evweight);
                        l.rooContainer->InputBinnedDataPoint(Form("sig_%dlow_BDT_grad_",sideband_i)+currentTypeSignalLabel ,category,bdt_grad,evweight);
                        }
                        }
                        }
                    */

                    // Loop over N higher sidebands
                    /*  for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){

                        double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                        double mass_hypothesis_high     = (mass_hypothesis*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                        double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
                        double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);

                        if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_high,bdtoutput,evweight);
                        bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                        bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
        
                        if (!VBFevent){ // note if includeVBF=false, VBFevent will always be false so no need to check both
                        l.rooContainer->InputBinnedDataPoint(Form("sig_%dhigh_BDT_ada_",sideband_i)+currentTypeSignalLabel ,category,bdt_ada,evweight);
                        l.rooContainer->InputBinnedDataPoint(Form("sig_%dhigh_BDT_grad_",sideband_i)+currentTypeSignalLabel ,category,bdt_grad,evweight);
                        }
                        }
                        }
                    */
                }


            }
            // ---- Now deal with background MC and data
            else {

                // Iterate over each BDT mass. 
                for (int i = 1; i<nMasses;i++){ //ignoring mass 105 for now

                    double mass_h_low;      
                    double mass_h_high;

                    if (i==7){ // the 140 case because currently 145 is missing
                        mass_h_low  =-2.5;
                        mass_h_high =4.6;
                    } else if (i==8){  // the 150 case
                        mass_h_low  =-5.0;
                        mass_h_high =0.1;
                    } else {
                        mass_h_low  =-2.5;
                        mass_h_high =2.1;
                    }

                    for (double mass_h=mass_h_low; mass_h<mass_h_high; mass_h+=0.5){

                        float mass_hypothesis = masses[i]+mass_h;
                        if (mass_hypothesis < masses[1] || mass_hypothesis > masses[nMasses-1]) continue;

                        float sideband_boundaries[2];
                        sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                        sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                        //Signal Window
                        if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                            //std::cout << "sig region" << std::endl;
                            histoplace=1;
                            SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight);

                            bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                            bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                            if (mass_hypothesis == masses[i]) {

                                if (bdtoutput>=0.05) {
                                    l.FillHist("pt_msig",i-1, Higgs.Pt(), evweight);
                                    l.FillHist("logpt_msig",i-1, log10(Higgs.Pt()), evweight);
                                    l.FillHist("ptOverM_msig",i-1, Higgs.Pt()/mass, evweight);
                                    l.FillHist("ptOverMH_msig",i-1, _H_ptOverM, evweight);
                                    l.FillHist("eta_msig",i-1, Higgs.Eta(), evweight);
                                    l.FillHist("deltaPhi_msig",i-1, _d_phi, evweight);
                                    l.FillHist("cosDeltaPhi_msig",i-1, cos(_d_phi), evweight);
                                    l.FillHist("deltaEta_msig",i-1, _d_eta, evweight);
                                    l.FillHist("helicityAngle_msig",i-1, _cos_theta_star, evweight);
                                    l.FillHist("pho1_pt_msig",i-1,lead_p4.Pt(), evweight);
                                    l.FillHist("pho2_pt_msig",i-1,sublead_p4.Pt(), evweight);
                                    l.FillHist("pho1_ptOverM_msig",i-1,lead_p4.Pt()/mass, evweight);
                                    l.FillHist("pho2_ptOverM_msig",i-1,sublead_p4.Pt()/mass, evweight);
                                    l.FillHist("pho1_ptOverMH_msig",i-1,_pho1_ptOverM, evweight);
                                    l.FillHist("pho2_ptOverMH_msig",i-1,_pho2_ptOverM, evweight);
                                    l.FillHist("pho1_eta_msig",i-1,lead_p4.Eta(), evweight);
                                    l.FillHist("pho2_eta_msig",i-1,sublead_p4.Eta(), evweight);
                                    l.FillHist("pho_minr9_msig",i-1,min(l.pho_r9[diphoton_index.first],l.pho_r9[diphoton_index.second]), evweight);
                                    l.FillHist("pho1_r9_msig",i-1,l.pho_r9[diphoton_index.first], evweight);
                                    l.FillHist("pho2_r9_msig",i-1,l.pho_r9[diphoton_index.second], evweight);
                                    l.FillHist("maxeta_msig",i-1,_max_eta, evweight);
                                    l.FillHist("deltaMOverMH_msig",i-1, _deltaMOverM, evweight);
                                    l.FillHist("sigmaM_msig",i-1,sigmaMrv, evweight);
                                    l.FillHist("sigmaM_wrongVtx_msig",i-1,sigmaMwv, evweight);
                                    l.FillHist("sigmaMOverM_msig",i-1,sigmaMrv/mass, evweight);
                                    l.FillHist("sigmaMOverM_wrongVtx_msig",i-1,sigmaMwv/mass, evweight);
                                    l.FillHist("sigmaMOverMH_msig",i-1,sigmaMrv/mass_hypothesis, evweight);
                                    l.FillHist("sigmaMOverMH_wrongVtx_msig",i-1,sigmaMwv/mass_hypothesis, evweight);
                                    l.FillHist("deltaMOverSigmaM_msig",i-1,_deltaMOverSigmaM, evweight);
                                    l.FillHist("deltaMOverSigmaM_wrongVtx_msig",i-1,(mass-mass_hypothesis)/sigmaMwv, evweight);
                                    l.FillHist("pho1_phoidMva_msig",i-1,phoid_mvaout_lead, evweight);
                                    l.FillHist("pho2_phoidMva_msig",i-1,phoid_mvaout_sublead, evweight);
                                    l.FillHist("vtxProb_msig",i-1,vtxProb, evweight);
                                }
                                l.FillHist("bdtoutput_msig",i-1,bdtoutput, evweight);

                            }

                            //std::cout << "ada: " << bdt_ada << "  grad: " << bdt_grad << std::endl;

                            if (cur_type == 0 ){//data
    
                                if (cur_type==0 && mass_hypothesis==124.0 && bdtoutput>=0.05){
                                    eventListText <<"Type="<< cur_type << " diphotonBDT="<<bdtoutput<<" mgg="<<mass<<" bdt="<<bdt_grad<<" vbf=" << VBFevent << endl;
                                }
                                if (VBFevent){ 
                                    l.rooContainer->InputBinnedDataPoint(Form("data_VBF_%3.1f",mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                } else {
                                    l.rooContainer->InputBinnedDataPoint(Form("data_BDT_ada_%3.1f",mass_hypothesis),category,bdt_ada,evweight);
                                    l.rooContainer->InputBinnedDataPoint(Form("data_BDT_grad_%3.1f",mass_hypothesis),category,bdt_grad,evweight);
                                }
                            }
                            if (cur_type > 0 ){// background MC
                                if (cur_type==6){
                                    //l.rooContainer->InputBinnedDataPoint(Form("zee_BDT_ada_%3.1f",mass_hypothesis),category,bdt_ada,evweight);
                                    //l.rooContainer->InputBinnedDataPoint(Form("zee_BDT_grad_%3.1f",mass_hypothesis) ,category,bdt_grad,evweight);
                                }
                                else{
                                    if (VBFevent){
                                        l.rooContainer->InputBinnedDataPoint(Form("bkg_VBF_%3.1f",mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                    } else {
                                        l.rooContainer->InputBinnedDataPoint(Form("bkg_BDT_ada_%3.1f",mass_hypothesis),category,bdt_ada,evweight);
                                        l.rooContainer->InputBinnedDataPoint(Form("bkg_BDT_grad_%3.1f",mass_hypothesis) ,category,bdt_grad,evweight);
                                    }
                                }
                            }
                        }

                        else {
                            // Loop over N lower sidebands
                            std::string sideband_str;
                            std::stringstream sideband_stream;
                            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){
                                histoplace=sideband_i+1;
                                double hypothesisModifier = (1.-sidebandWidth)/(1+sidebandWidth);
                                double mass_hypothesis_low     = (mass_hypothesis*(1.-signalRegionWidth)/(1.+sidebandWidth)-sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                                double sideband_boundaries_low = mass_hypothesis_low*(1.-sidebandWidth);
                                double sideband_boundaries_high= mass_hypothesis_low*(1.+sidebandWidth);

                                //cout << "sideband, "<< sideband_i << " mH, "<<mass_hypothesis_low<< " bands, " <<sideband_boundaries_low << " " << sideband_boundaries_high <<endl;
                                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_low,bdtoutput,evweight);
                                    bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                                    bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                                    if (mass_hypothesis == masses[i]) {

                                        if (sideband_i==1) {
                                            sideband_str="";
                                        } else {
                                            sideband_stream << sideband_i;
                                            sideband_str = sideband_stream.str();
                                        }

                                        if (bdtoutput>=0.05) {
                                            l.FillHist("pt_mlow"+sideband_str,i-1, Higgs.Pt(), evweight);
                                            l.FillHist("logpt_mlow"+sideband_str,i-1, log10(Higgs.Pt()), evweight);
                                            l.FillHist("ptOverM_mlow"+sideband_str,i-1, Higgs.Pt()/mass, evweight);
                                            l.FillHist("ptOverMH_mlow"+sideband_str,i-1, _H_ptOverM, evweight);
                                            l.FillHist("eta_mlow"+sideband_str,i-1, Higgs.Eta(), evweight);
                                            l.FillHist("deltaPhi_mlow"+sideband_str,i-1, _d_phi, evweight);
                                            l.FillHist("cosDeltaPhi_mlow"+sideband_str,i-1, cos(_d_phi), evweight);
                                            l.FillHist("deltaEta_mlow"+sideband_str,i-1, _d_eta, evweight);
                                            l.FillHist("helicityAngle_mlow"+sideband_str,i-1, _cos_theta_star, evweight);
                                            l.FillHist("pho1_pt_mlow"+sideband_str,i-1,lead_p4.Pt(), evweight);
                                            l.FillHist("pho2_pt_mlow"+sideband_str,i-1,sublead_p4.Pt(), evweight);
                                            l.FillHist("pho1_ptOverM_mlow"+sideband_str,i-1,lead_p4.Pt()/mass, evweight);
                                            l.FillHist("pho2_ptOverM_mlow"+sideband_str,i-1,sublead_p4.Pt()/mass, evweight);
                                            l.FillHist("pho1_ptOverMH_mlow"+sideband_str,i-1,_pho1_ptOverM, evweight);
                                            l.FillHist("pho2_ptOverMH_mlow"+sideband_str,i-1,_pho2_ptOverM, evweight);
                                            l.FillHist("pho1_eta_mlow"+sideband_str,i-1,lead_p4.Eta(), evweight);
                                            l.FillHist("pho2_eta_mlow"+sideband_str,i-1,sublead_p4.Eta(), evweight);
                                            l.FillHist("pho_minr9_mlow"+sideband_str,i-1,_min_r9, evweight);
                                            l.FillHist("pho1_r9_mlow"+sideband_str,i-1,l.pho_r9[diphoton_index.first], evweight);
                                            l.FillHist("pho2_r9_mlow"+sideband_str,i-1,l.pho_r9[diphoton_index.second], evweight);
                                            l.FillHist("maxeta_mlow"+sideband_str,i-1,_max_eta, evweight);
                                            l.FillHist("deltaMOverMH_mlow"+sideband_str,i-1, _deltaMOverM, evweight);
                                            l.FillHist("sigmaM_mlow"+sideband_str,i-1,sigmaMrv, evweight);
                                            l.FillHist("sigmaM_wrongVtx_mlow"+sideband_str,i-1,sigmaMwv, evweight);
                                            l.FillHist("sigmaMOverM_mlow"+sideband_str,i-1,sigmaMrv/mass, evweight);
                                            l.FillHist("sigmaMOverM_wrongVtx_mlow"+sideband_str,i-1,sigmaMwv/mass, evweight);
                                            l.FillHist("sigmaMOverMH_mlow"+sideband_str,i-1,sigmaMrv/mass_hypothesis, evweight);
                                            l.FillHist("sigmaMOverMH_wrongVtx_mlow"+sideband_str,i-1,sigmaMwv/mass_hypothesis, evweight);
                                            l.FillHist("deltaMOverSigmaM_mlow"+sideband_str,i-1,_deltaMOverSigmaM, evweight);
                                            l.FillHist("deltaMOverSigmaM_wrongVtx_mlow"+sideband_str,i-1,(mass-mass_hypothesis)/sigmaMwv, evweight);
                                            l.FillHist("pho1_phoidMva_mlow"+sideband_str,i-1,phoid_mvaout_lead, evweight);
                                            l.FillHist("pho2_phoidMva_mlow"+sideband_str,i-1,phoid_mvaout_sublead, evweight);
                                            l.FillHist("vtxProb_mlow"+sideband_str,i-1,vtxProb, evweight);
                                        }
                                        l.FillHist("bdtoutput_mlow"+sideband_str,i-1,bdtoutput, evweight);

                                    }

                                    if (cur_type == 0 ){//data
                                        if (VBFevent){
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dlow_VBF_%3.1f",sideband_i,mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                        } else {
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                                        }
                                    }
                                    else if (cur_type > 0 ){// background MC
                                        if (cur_type==6){
                                            //l.rooContainer->InputBinnedDataPoint(Form("zee_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                            //l.rooContainer->InputBinnedDataPoint(Form("zee_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                                        } else {
                                            if (VBFevent){
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dlow_VBF_%3.1f",sideband_i,mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                            } else {
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dlow_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dlow_BDT_grad_%3.1f",sideband_i,mass_hypothesis),category,bdt_grad,evweight);   
                                            }
                                        }
                                    }
                                }
                            }

                            // Loop over N higher sidebands
                            for (int sideband_i = 1 ; sideband_i <= numberOfSidebands ; sideband_i++){

                                histoplace=sideband_i+4;
                                double hypothesisModifier = (1.+sidebandWidth)/(1-sidebandWidth);
                                double mass_hypothesis_high     = (mass_hypothesis*(1.+signalRegionWidth)/(1.-sidebandWidth)+sidebandShift)*(TMath::Power(hypothesisModifier,sideband_i-1));
                                double sideband_boundaries_low = mass_hypothesis_high*(1.-sidebandWidth);
                                double sideband_boundaries_high= mass_hypothesis_high*(1.+sidebandWidth);

                                if ( mass>sideband_boundaries_low && mass<sideband_boundaries_high){
                                    SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis_high,bdtoutput,evweight);
                                    bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                                    bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                                    if (mass_hypothesis == masses[i]) {

                                        if (sideband_i==1) {
                                            sideband_str="";
                                        } else {
                                            sideband_stream << sideband_i;
                                            sideband_str = sideband_stream.str();
                                        }

                                        if (bdtoutput>=0.05) {
                                            l.FillHist("pt_mhigh"+sideband_str,i-1, Higgs.Pt(), evweight);
                                            l.FillHist("logpt_mhigh"+sideband_str,i-1, log10(Higgs.Pt()), evweight);
                                            l.FillHist("ptOverM_mhigh"+sideband_str,i-1, Higgs.Pt()/mass, evweight);
                                            l.FillHist("ptOverMH_mhigh"+sideband_str,i-1, _H_ptOverM, evweight);
                                            l.FillHist("eta_mhigh"+sideband_str,i-1, Higgs.Eta(), evweight);
                                            l.FillHist("deltaPhi_mhigh"+sideband_str,i-1, _d_phi, evweight);
                                            l.FillHist("cosDeltaPhi_mhigh"+sideband_str,i-1, cos(_d_phi), evweight);
                                            l.FillHist("deltaEta_mhigh"+sideband_str,i-1, _d_eta, evweight);
                                            l.FillHist("helicityAngle_mhigh"+sideband_str,i-1, _cos_theta_star, evweight);
                                            l.FillHist("pho1_pt_mhigh"+sideband_str,i-1,lead_p4.Pt(), evweight);
                                            l.FillHist("pho2_pt_mhigh"+sideband_str,i-1,sublead_p4.Pt(), evweight);
                                            l.FillHist("pho1_ptOverM_mhigh"+sideband_str,i-1,lead_p4.Pt()/mass, evweight);
                                            l.FillHist("pho2_ptOverM_mhigh"+sideband_str,i-1,sublead_p4.Pt()/mass, evweight);
                                            l.FillHist("pho1_ptOverMH_mhigh"+sideband_str,i-1,_pho1_ptOverM, evweight);
                                            l.FillHist("pho2_ptOverMH_mhigh"+sideband_str,i-1,_pho2_ptOverM, evweight);
                                            l.FillHist("pho1_eta_mhigh"+sideband_str,i-1,lead_p4.Eta(), evweight);
                                            l.FillHist("pho2_eta_mhigh"+sideband_str,i-1,sublead_p4.Eta(), evweight);
                                            l.FillHist("pho_minr9_mhigh"+sideband_str,i-1,_min_r9, evweight);
                                            l.FillHist("pho1_r9_mhigh"+sideband_str,i-1,l.pho_r9[diphoton_index.first], evweight);
                                            l.FillHist("pho2_r9_mhigh"+sideband_str,i-1,l.pho_r9[diphoton_index.second], evweight);
                                            l.FillHist("maxeta_mhigh"+sideband_str,i-1,_max_eta, evweight);
                                            l.FillHist("deltaMOverMH_mhigh"+sideband_str,i-1, _deltaMOverM, evweight);
                                            l.FillHist("sigmaM_mhigh"+sideband_str,i-1,sigmaMrv, evweight);
                                            l.FillHist("sigmaM_wrongVtx_mhigh"+sideband_str,i-1,sigmaMwv, evweight);
                                            l.FillHist("sigmaMOverM_mhigh"+sideband_str,i-1,sigmaMrv/mass, evweight);
                                            l.FillHist("sigmaMOverM_wrongVtx_mhigh"+sideband_str,i-1,sigmaMwv/mass, evweight);
                                            l.FillHist("sigmaMOverMH_mhigh"+sideband_str,i-1,sigmaMrv/mass_hypothesis, evweight);
                                            l.FillHist("sigmaMOverMH_wrongVtx_mhigh"+sideband_str,i-1,sigmaMwv/mass_hypothesis, evweight);
                                            l.FillHist("deltaMOverSigmaM_mhigh"+sideband_str,i-1,_deltaMOverSigmaM, evweight);
                                            l.FillHist("deltaMOverSigmaM_wrongVtx_mhigh"+sideband_str,i-1,(mass-mass_hypothesis)/sigmaMwv, evweight);
                                            l.FillHist("pho1_phoidMva_mhigh"+sideband_str,i-1,phoid_mvaout_lead, evweight);
                                            l.FillHist("pho2_phoidMva_mhigh"+sideband_str,i-1,phoid_mvaout_sublead, evweight);
                                            l.FillHist("vtxProb_mhigh"+sideband_str,i-1,vtxProb, evweight);
                                        }
                                        l.FillHist("bdtoutput_mhigh"+sideband_str,i-1,bdtoutput, evweight);

                                    }

                                    if (cur_type == 0 ){//data
                                        if (VBFevent){
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dhigh_VBF_%3.1f",sideband_i,mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                        } else {
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                            l.rooContainer->InputBinnedDataPoint(Form("data_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                                        }
                                    }
                                    else if (cur_type > 0 ){// background MC
                                        if (cur_type==6){
                                            //l.rooContainer->InputBinnedDataPoint(Form("zee_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                            //l.rooContainer->InputBinnedDataPoint(Form("zee_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis) ,category,bdt_grad,evweight);
                                        } else {
                                            if (VBFevent){
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dhigh_VBF_%3.1f",sideband_i,mass_hypothesis),category,1.+sidebandWidth+_deltaMOverM,evweight);
                                            } else {
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dhigh_BDT_ada_%3.1f",sideband_i,mass_hypothesis),category,bdt_ada,evweight);
                                                l.rooContainer->InputBinnedDataPoint(Form("bkg_%dhigh_BDT_grad_%3.1f",sideband_i,mass_hypothesis),category,bdt_grad,evweight);
                                            }
                                        }
                                    }
                                }
                            }
                        } 
                    } // end loop over mass hypotheses
                } // end loop over BDT
            } // Background MC or Data

            // histos, 0 -> default, 1->signal region, 2,3,4 -> lower SB , 5,6,7 -> higher SB,8 all
            if (mass>massMin && mass<massMax){

                // Special histogram fill BDT_mH, BDT_mH + .5GeV
                // ------------------------------------------------------------------------------------------------------------------------------------//
                float q_mass_hypothesis = 120.0;
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis,bdtoutput,evweight);
                float bdt_grad_1 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis+1.2,bdtoutput,evweight);
                float bdt_grad_2 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis-1.2,bdtoutput,evweight);
                float bdt_grad_3 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                if (bdtoutput>=0.05){
                    l.FillHist2D("shiftingMH_up_bdt1_bdt2_mH120_1p"  ,0,bdt_grad_1,bdt_grad_2,evweight);
                    l.FillHist2D("shiftingMH_dn_bdt1_bdt2_mH120_1p"  ,0,bdt_grad_1,bdt_grad_3,evweight);
                }

                q_mass_hypothesis = 140.0;
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis,bdtoutput,evweight);
                bdt_grad_1 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis+1.4,bdtoutput,evweight);
                bdt_grad_2 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis-1.4,bdtoutput,evweight);
                bdt_grad_3 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                if (bdtoutput>=0.05){
                    l.FillHist2D("shiftingMH_up_bdt1_bdt2_mH140_1p"  ,0,bdt_grad_1,bdt_grad_2,evweight);
                    l.FillHist2D("shiftingMH_dn_bdt1_bdt2_mH140_1p"  ,0,bdt_grad_1,bdt_grad_3,evweight);
                }

                q_mass_hypothesis = 120.0;
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis,bdtoutput,evweight);
                bdt_grad_1 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis+2.4,bdtoutput,evweight);
                bdt_grad_2 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis-2.4,bdtoutput,evweight);
                bdt_grad_3 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                if (bdtoutput>=0.05){
                    l.FillHist2D("shiftingMH_up_bdt1_bdt2_mH120_2p"  ,0,bdt_grad_1,bdt_grad_2,evweight);
                    l.FillHist2D("shiftingMH_dn_bdt1_bdt2_mH120_2p"  ,0,bdt_grad_1,bdt_grad_3,evweight);
                }

                q_mass_hypothesis = 140.0;
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis,bdtoutput,evweight);
                bdt_grad_1 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis+2.8,bdtoutput,evweight);
                bdt_grad_2 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,q_mass_hypothesis-2.8,bdtoutput,evweight);
                bdt_grad_3 = tmvaReader_->EvaluateMVA( "BDT_grad_123" );
                if (bdtoutput>=0.05){
                    l.FillHist2D("shiftingMH_up_bdt1_bdt2_mH140_2p"  ,0,bdt_grad_1,bdt_grad_2,evweight);
                    l.FillHist2D("shiftingMH_dn_bdt1_bdt2_mH140_2p"  ,0,bdt_grad_1,bdt_grad_3,evweight);
                }
                // ------------------------------------------------------------------------------------------------------------------------------------//

                /*
                  l.FillHist2D("bdtgrad_vs_mass"    ,histoplace,bdt_grad,mass,evweight);
                  l.FillHist2D("bdtgrad_vs_hpt"    ,histoplace,bdt_grad,ptHiggs,evweight);
                  l.FillHist2D("bdtgrad_vs_leadpt"  ,histoplace,bdt_grad,pt_lead,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadpt"  ,histoplace,bdt_grad,pt_sublead,evweight);
                  l.FillHist2D("bdtgrad_vs_leadeta"  ,histoplace,bdt_grad,eta_lead,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadeta"  ,histoplace,bdt_grad,eta_sublead,evweight);
                  l.FillHist2D("bdtgrad_vs_dphi"    ,histoplace,bdt_grad,delta_phi,evweight);
                  l.FillHist2D("bdtgrad_vs_heta"    ,histoplace,bdt_grad,etaHiggs,evweight);
                  l.FillHist2D("bdtgrad_vs_sigmamoverm"  ,histoplace,bdt_grad,sigmaMrv/mass,evweight);
                  l.FillHist2D("bdtgrad_vs_leadr9"  ,histoplace,bdt_grad,lead_r9,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadr9"  ,histoplace,bdt_grad,sublead_r9,evweight);

                  l.FillHist2D("bdtada_vs_mass"    ,histoplace,bdt_ada,mass,evweight);
                  l.FillHist2D("bdtada_vs_hpt"    ,histoplace,bdt_ada,ptHiggs,evweight);
                  l.FillHist2D("bdtada_vs_leadpt"  ,histoplace,bdt_ada,pt_lead,evweight);
                  l.FillHist2D("bdtada_vs_subleadpt"  ,histoplace,bdt_ada,pt_sublead,evweight);
                  l.FillHist2D("bdtada_vs_leadeta"  ,histoplace,bdt_ada,eta_lead,evweight);
                  l.FillHist2D("bdtada_vs_subleadeta"  ,histoplace,bdt_ada,eta_sublead,evweight);
                  l.FillHist2D("bdtada_vs_dphi"    ,histoplace,bdt_ada,delta_phi,evweight);
                  l.FillHist2D("bdtada_vs_heta"    ,histoplace,bdt_ada,etaHiggs,evweight);
                  l.FillHist2D("bdtada_vs_sigmamoverm"  ,histoplace,bdt_ada,sigmaMrv/mass,evweight);
                  l.FillHist2D("bdtada_vs_leadr9"  ,histoplace,bdt_ada,lead_r9,evweight);
                  l.FillHist2D("bdtada_vs_subleadr9"  ,histoplace,bdt_ada,sublead_r9,evweight);

                  l.FillHist2D("bdtgrad_vs_mass"    ,8,bdt_grad,mass,evweight);
                  l.FillHist2D("bdtgrad_vs_hpt"    ,8,bdt_grad,ptHiggs,evweight);
                  l.FillHist2D("bdtgrad_vs_leadpt"  ,8,bdt_grad,pt_lead,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadpt"  ,8,bdt_grad,pt_sublead,evweight);
                  l.FillHist2D("bdtgrad_vs_leadeta"  ,8,bdt_grad,eta_lead,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadeta"  ,8,bdt_grad,eta_sublead,evweight);
                  l.FillHist2D("bdtgrad_vs_dphi"    ,8,bdt_grad,delta_phi,evweight);
                  l.FillHist2D("bdtgrad_vs_heta"    ,8,bdt_grad,etaHiggs,evweight);
                  l.FillHist2D("bdtgrad_vs_sigmamoverm"  ,8,bdt_grad,sigmaMrv/mass,evweight);
                  l.FillHist2D("bdtgrad_vs_leadr9"  ,8,bdt_grad,lead_r9,evweight);
                  l.FillHist2D("bdtgrad_vs_subleadr9"  ,8,bdt_grad,sublead_r9,evweight);

                  l.FillHist2D("bdtada_vs_mass"    ,8,bdt_ada,mass,evweight);
                  l.FillHist2D("bdtada_vs_hpt"    ,8,bdt_ada,ptHiggs,evweight);
                  l.FillHist2D("bdtada_vs_leadpt"  ,8,bdt_ada,pt_lead,evweight);
                  l.FillHist2D("bdtada_vs_subleadpt"  ,8,bdt_ada,pt_sublead,evweight);
                  l.FillHist2D("bdtada_vs_leadeta"  ,8,bdt_ada,eta_lead,evweight);
                  l.FillHist2D("bdtada_vs_subleadeta"  ,8,bdt_ada,eta_sublead,evweight);
                  l.FillHist2D("bdtada_vs_dphi"    ,8,bdt_ada,delta_phi,evweight);
                  l.FillHist2D("bdtada_vs_heta"    ,8,bdt_ada,etaHiggs,evweight);
                  l.FillHist2D("bdtada_vs_sigmamoverm"  ,8,bdt_ada,sigmaMrv/mass,evweight);
                  l.FillHist2D("bdtada_vs_leadr9"  ,8,bdt_ada,lead_r9,evweight);
                  l.FillHist2D("bdtada_vs_subleadr9"  ,8,bdt_ada,sublead_r9,evweight);
                */
            }
            l.FillCounter( "Accepted", weight );
            l.FillCounter( "Smeared", evweight );
            sumaccept += weight;
            sumsmear += evweight;

        }
    }
    if(PADEBUG) 
        cout<<"myFillHistRed END"<<endl;


    // Now do some MC Systematic studies - For MVA Analysis, assume that we only care about Signal Systematics Here, ie make the check cur_type < 0 to increase speed
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
    if( cur_type < 0 && doMCSmearing && !doTraining){
        // fill steps for syst uncertainty study
        float systStep = systRange / (float)nSystSteps;
        // di-photon smearers systematics (These frst smeaers will not effect the selection so can reuse the diphoton_id)
        if (diphoton_id > -1 ) {

            float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
            float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
            TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
            TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
            TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);

            for(std::vector<BaseGenLevelSmearer*>::iterator si=systGenLevelSmearers_.begin(); si!=systGenLevelSmearers_.end(); si++){
                std::vector<double> bdt_grad_errors;
                std::vector<double> bdt_ada_errors;
                std::vector<double> vbf_values;
                std::vector<double> weights;
                std::vector<int>    categories;

                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
                    TLorentzVector Higgs = lead_p4 + sublead_p4;   

                    int category = 0;  // Category always 0 fro MVA analysis
                    double genLevWeightSyst=1; 

                    for(std::vector<BaseGenLevelSmearer *>::iterator sj=genLevelSmearers_.begin(); sj!= genLevelSmearers_.end(); ++sj ) {
                        float swei=1.;
                        if( *si == *sj ) { 
                            (*si)->smearEvent(swei, gP4, l.pu_n, cur_type, syst_shift );
                        } else {
                            (*sj)->smearEvent(swei, gP4, l.pu_n, cur_type, 0. );
                        }
                        genLevWeightSyst *= swei;
                    }

                    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeightSyst;
                    float mass = Higgs.M();
                    float ptHiggs = Higgs.Pt();

                    // Mass Resolution of the Event
                    //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
                    massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
                    //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
                    //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

                    double massResolution = massResolutionCalculator->massResolutionCorrVtx();
                    double vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
                    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
                    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
                    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
                    float vtxProb = 1.-0.49*(vtx_mva+1.0);

                    float bdtoutput;
                    if (bdtTrainingPhilosophy=="MIT"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"MIT");
                        if (bdtoutput < 0.05 ) category = -1;
                    } 
                    else if (bdtTrainingPhilosophy=="UCSD"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"UCSD");
                    }

                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[SignalType(cur_type)];
                    // define the sidebands
                    float sideband_boundaries[2];
                    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                    //Signal Window
                    if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight);
                        float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                        float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                        categories.push_back(category);
                        if (VBFevent){
                            bdt_ada_errors.push_back(-100);
                            bdt_grad_errors.push_back(-100);    
                            vbf_values.push_back(1.5);
                        } else {
                            bdt_ada_errors.push_back(bdt_ada);
                            bdt_grad_errors.push_back(bdt_grad);
                            vbf_values.push_back(-100);
                        }
                        weights.push_back(evweight);

                    } else {
                        categories.push_back(-1);
                        bdt_ada_errors.push_back(-100.);
                        bdt_grad_errors.push_back(-100.);
                        vbf_values.push_back(-100);
                        weights.push_back(0.);
                    }

                }// end loop on systematics steps
                // Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                if (includeVBF) l.rooContainer->InputSystematicSet("sig_VBF_"+currentTypeSignalLabel,(*si)->name(),categories,vbf_values,weights);
                l.rooContainer->InputSystematicSet("sig_BDT_ada_"+currentTypeSignalLabel,(*si)->name(),categories,bdt_ada_errors,weights);
                l.rooContainer->InputSystematicSet("sig_BDT_grad_"+currentTypeSignalLabel,(*si)->name(),categories,bdt_grad_errors,weights);
                // -------------------------------------------------------------------------------------------------------------//

            } // end loop on smearers 


            for(std::vector<BaseDiPhotonSmearer *>::iterator si=systDiPhotonSmearers_.begin(); si!= systDiPhotonSmearers_.end(); ++si ) {
                std::vector<double> bdt_grad_errors;
                std::vector<double> bdt_ada_errors;
                std::vector<double> vbf_values;
                std::vector<double> weights;
                std::vector<int> categories;

                for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                    if( syst_shift == 0. ) { continue; } // skip the central value
                    TLorentzVector Higgs = lead_p4 + sublead_p4;   

                    // restart with 'fresh' wait for this round of systematics
                    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] * genLevWeight;

                    // FIXME pass smeared R9 and di-photon
                    int category = 0; 
                    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);

                    float photon_idMVA1=l.photonIDMVA(diphoton_index.first,l.dipho_vtxind[diphoton_id],lead_p4,"MIT");
                    float photon_idMVA2=l.photonIDMVA(diphoton_index.second,l.dipho_vtxind[diphoton_id],sublead_p4,"MIT");

                    for(std::vector<BaseDiPhotonSmearer *>::iterator sj=diPhotonSmearers_.begin(); sj!= diPhotonSmearers_.end(); ++sj ) {
                        float swei=1.;
                        float pth = Higgs.Pt();
                        if( *si == *sj ) { 
                            (*si)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), photon_idMVA1,photon_idMVA2,syst_shift);
                        } else { 
                            (*sj)->smearDiPhoton( Higgs, *vtx, swei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)) ,photon_idMVA1,photon_idMVA2,0.);
                        }
                        evweight *= swei;
                    }
                    float mass = Higgs.M();
                    float ptHiggs = Higgs.Pt();
                    // Mass Resolution of the Event
                    //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
                    massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
                    //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
                    //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

                    double massResolution = massResolutionCalculator->massResolutionCorrVtx();
                    double vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
                    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
                    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
                    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
                    float vtxProb = 1.-0.49*(vtx_mva+1.0);

                    float bdtoutput;
                    if (bdtTrainingPhilosophy=="MIT"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"MIT",photon_idMVA1,photon_idMVA2);
                        if (bdtoutput < 0.05) category = -1;
                    } 
                    else if (bdtTrainingPhilosophy=="UCSD"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"UCSD");
                    }

                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[SignalType(cur_type)];
                    // define the sidebands
                    float sideband_boundaries[2];
                    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                    //Signal Window
                    if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight);
                        float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                        float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                        categories.push_back(category);
                        if (VBFevent){
                            bdt_ada_errors.push_back(-100);
                            bdt_grad_errors.push_back(-100);    
                            vbf_values.push_back(1.5);
                        } else {
                            bdt_ada_errors.push_back(bdt_ada);
                            bdt_grad_errors.push_back(bdt_grad);
                            vbf_values.push_back(-100);
                        }
                        weights.push_back(evweight);

                    } else {
                        categories.push_back(-1);
                        bdt_ada_errors.push_back(-100.);
                        bdt_grad_errors.push_back(-100.);
                        vbf_values.push_back(-100);
                        weights.push_back(0.);
                    }
                }// end loop on systematics steps

                // Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
                if (includeVBF)  l.rooContainer->InputSystematicSet("sig_VBF_" +currentTypeSignalLabel,(*si)->name(),categories,vbf_values,weights);
                l.rooContainer->InputSystematicSet("sig_BDT_ada_" +currentTypeSignalLabel,(*si)->name(),categories,bdt_ada_errors,weights);
                l.rooContainer->InputSystematicSet("sig_BDT_grad_"+currentTypeSignalLabel,(*si)->name(),categories,bdt_grad_errors,weights);
                // -------------------------------------------------------------------------------------------------------------//

            } // end loop on smearers 

        } // Close If on CiC Selection

        // -------------------------------------------------------------------------------------------------------------------------------------//
        // Now the systematic Steps which can effect the CiC Selection -------------------------------------------------------------------------//

        // loop over the smearers included in the systematics study
        for(std::vector<BaseSmearer *>::iterator  si=systPhotonSmearers_.begin(); si!= systPhotonSmearers_.end(); ++si ) {
            std::vector<double> bdt_grad_errors;
            std::vector<double> bdt_ada_errors;
            std::vector<double> vbf_values;
            std::vector<double> weights;
            std::vector<int> categories;

            // loop over syst shift
            for(float syst_shift=-systRange; syst_shift<=systRange; syst_shift+=systStep ) { 
                if( syst_shift == 0. ) { continue; } // skip the central value
                // smear the photons 
                photonInfoCollection.clear();
                for(int ipho=0; ipho<l.pho_n; ++ipho ) { 
                    std::vector<std::vector<bool> > p;
                    PhotonReducedInfo phoInfo (/// *((TVector3*)l.pho_calopos->At(ipho)), 
                                               *((TVector3*)l.sc_xyz->At(l.pho_scind[ipho])), 
                                               ((TLorentzVector*)l.pho_p4->At(ipho))->Energy(), 
                                               energyCorrected[ipho],
                                               l.pho_isEB[ipho], l.pho_r9[ipho],
                                               l.PhotonCiCSelectionLevel(ipho,l.vtx_std_sel,p,nPhotonCategories_),
                                               (energyCorrectedError!=0?energyCorrectedError[ipho]:0)
                                               );
                    if (l.CheckSphericalPhoton(ipho)) phoInfo.setSphericalPhoton(true);

                    float pweight = 1.;
                    for(std::vector<BaseSmearer *>::iterator  sj=photonSmearers_.begin(); sj!= photonSmearers_.end(); ++sj ) {
                        float sweight = 1.;
                        if( *si == *sj ) {
                            // move the smearer under study by syst_shift
                            (*si)->smearPhoton(phoInfo,sweight,l.run,syst_shift);
                        } else {
                            // for the other use the nominal points
                            (*sj)->smearPhoton(phoInfo,sweight,l.run,0.);
                        }
                        pweight *= sweight;
                    }
                    smeared_pho_energy[ipho] = phoInfo.energy();
                    smeared_pho_r9[ipho] = phoInfo.r9();
                    smeared_pho_weight[ipho] = pweight;
                    photonInfoCollection.push_back(phoInfo);
                }

                // analyze the event
                // FIXME pass smeared R9
                //    int diphoton_id = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
                // VBF-TAGGING -------------------------------------------------------------------------- //
                VBFevent = false;

                int diphoton_id=-1;

                if(includeVBF) {
                    int diphoton_id_vbf = l.DiphotonMITPreSelection(leadEtCutVBF,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
                    if (diphoton_id_vbf >-1){
                        // JET MET Corrections // No need to reset the pointers again, so no need to RescaleJetEnergy a second time
                        diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id_vbf],  l.dipho_subleadind[diphoton_id_vbf] );

                        TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
                        TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id_vbf], l.dipho_vtxind[diphoton_id_vbf], &smeared_pho_energy[0]);
                        float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id_vbf]];
                        float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id_vbf]];
                        TLorentzVector Higgs = lead_p4 + sublead_p4;   

                        float jet1ptcut =0.0;
                        float jet2ptcut =0.0;
                        bool crosscheck = false;
                        std::pair<int,int> highestPtJets(-1,-1);

                        highestPtJets = l.Select2HighestPtJets(lead_p4, sublead_p4, jet1ptcut, jet2ptcut );
                        bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.second!=-1);

                        if(VBFpresel){

                            TLorentzVector* jet1 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.first);
                            TLorentzVector* jet2 = (TLorentzVector*)l.jet_algoPF1_p4->At(highestPtJets.second);
                            TLorentzVector dijet = (*jet1) + (*jet2);

                            myVBFLeadJPt = jet1->Pt();
                            myVBFSubJPt = jet2->Pt();
                            myVBF_Mjj = dijet.M();
                            myVBFdEta = fabs(jet1->Eta() - jet2->Eta());
                            myVBFZep  = fabs(Higgs.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
                            myVBFdPhi = fabs(Higgs.DeltaPhi(dijet));
                            myVBF_Mgg =Higgs.M();


                            // Cannot Get Apply cuts to work -> Need to discuss with C.Palmer, for now, jyst apply the cuts
                            //l.ApplyCutsFill(0,3,evweight, myweight);
                            //VBFevent = l.ApplyCutsFill(0,5,evweight, myweight);
                            //VBFevent = l.ApplyCutsFill(0,1,evweight, myweight);
                            VBFevent = l.ApplyCuts(0,1);
                            if(VBFevent) diphoton_id = diphoton_id_vbf;

                        }
                    }
                }
                // CP // NW VBF Tagging
                // --------------------- END VBF-TAGGING --------------------------------------------------------//
                if (diphoton_id < 0){ // failed to find a VBF 
                    if (bdtTrainingPhilosophy=="MIT"){
                        diphoton_id = l.DiphotonMITPreSelection(leadEtCut,subleadEtCut,applyPtoverM, &smeared_pho_energy[0] ); 
                    } else if (bdtTrainingPhilosophy=="UCSD"){
                        diphoton_id = l.DiphotonCiCSelection(l.phoLOOSE, l.phoLOOSE, leadEtCut, subleadEtCut, 4,applyPtoverM, &smeared_pho_energy[0] ); 
                    }
                }

                if (diphoton_id > -1 ) {

                    diphoton_index = std::make_pair( l.dipho_leadind[diphoton_id],  l.dipho_subleadind[diphoton_id] );
                    float evweight = weight * smeared_pho_weight[diphoton_index.first] * smeared_pho_weight[diphoton_index.second] *genLevWeight;

                    float lead_r9    = l.pho_r9[l.dipho_leadind[diphoton_id]];
                    float sublead_r9 = l.pho_r9[l.dipho_subleadind[diphoton_id]];
                    TLorentzVector lead_p4 = l.get_pho_p4( l.dipho_leadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
                    TLorentzVector sublead_p4 = l.get_pho_p4( l.dipho_subleadind[diphoton_id], l.dipho_vtxind[diphoton_id], &smeared_pho_energy[0]);
                    TVector3 * vtx = (TVector3*)l.vtx_std_xyz->At(l.dipho_vtxind[diphoton_id]);
                    TLorentzVector Higgs = lead_p4 + sublead_p4;   

                    int category = 0; 
                    int selectioncategory = l.DiphotonCategory(diphoton_index.first,diphoton_index.second,Higgs.Pt(),nEtaCategories,nR9Categories,0);
                    if( cur_type != 0 && doMCSmearing ) {
                        for(std::vector<BaseDiPhotonSmearer *>::iterator si=diPhotonSmearers_.begin(); si!= diPhotonSmearers_.end(); ++si ) {
                            float rewei=1.;
                            float pth = Higgs.Pt();
                            (*si)->smearDiPhoton( Higgs, *vtx, rewei, selectioncategory, cur_type, *((TVector3*)l.gv_pos->At(0)), zero_ ,zero_,0.);
                            evweight *= rewei;
                        }
                    }
                    float mass = Higgs.M();

                    float ptHiggs = Higgs.Pt();
                    // Mass Resolution of the Event
                    //massResolutionCalculator->Setup(l,&lead_p4,&sublead_p4,diphoton_index.first,diphoton_index.second,diphoton_id,ptHiggs,mass,eSmearPars,nR9Categories,nEtaCategories);
                    massResolutionCalculator->Setup(l,&photonInfoCollection[diphoton_index.first],&photonInfoCollection[diphoton_index.second],diphoton_id,eSmearPars,nR9Categories,nEtaCategories);
                    //massResolutionCalculator->setSphericalLeadPhoton(l.CheckSphericalPhoton(diphoton_index.first));
                    //massResolutionCalculator->setSphericalSubleadPhoton(l.CheckSphericalPhoton(diphoton_index.second));

                    double massResolution = massResolutionCalculator->massResolutionCorrVtx();
                    double vtx_mva = l.vtx_std_evt_mva->at(diphoton_id);
                    float sigmaMrv = massResolutionCalculator->massResolutionEonly();
                    float sigmaMwv = massResolutionCalculator->massResolutionWrongVtx();
                    float sigmaMeonly = massResolutionCalculator->massResolutionEonly();
                    float vtxProb = 1.-0.49*(vtx_mva+1.0);

                    float bdtoutput;
                    if (bdtTrainingPhilosophy=="MIT"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"MIT");
                        if (bdtoutput < 0.05) category = -1;  // Remove the VBF tagged events
                    } 
                    else if (bdtTrainingPhilosophy=="UCSD"){
                        bdtoutput = l.diphotonMVA(diphoton_index.first,diphoton_index.second,l.dipho_vtxind[diphoton_id],vtxProb,lead_p4,sublead_p4,sigmaMrv,sigmaMwv,sigmaMeonly,"UCSD");
                    }

                    // define hypothesis masses for the sidebands
                    float mass_hypothesis = masses[SignalType(cur_type)];
                    // define the sidebands
                    float sideband_boundaries[2];
                    sideband_boundaries[0] = mass_hypothesis*(1-sidebandWidth);
                    sideband_boundaries[1] = mass_hypothesis*(1+sidebandWidth);

                    //Signal Window
                    if( mass>sideband_boundaries[0] && mass<sideband_boundaries[1]){//Signal mass window cut
                        SetBDTInputVariables(&lead_p4,&sublead_p4,lead_r9,sublead_r9,massResolutionCalculator,vtx_mva,mass_hypothesis,bdtoutput,evweight);
                        float bdt_ada  = tmvaReader_->EvaluateMVA( "BDT_ada_123" );
                        float bdt_grad = tmvaReader_->EvaluateMVA( "BDT_grad_123" );

                        categories.push_back(category);
                        if (VBFevent){
                            bdt_ada_errors.push_back(-100);
                            bdt_grad_errors.push_back(-100);
                            vbf_values.push_back(1.5);
                        } else {
                            bdt_ada_errors.push_back(bdt_ada);
                            bdt_grad_errors.push_back(bdt_grad);
                            vbf_values.push_back(-100.);
                        }
                        weights.push_back(evweight);

                    } else {

                        categories.push_back(-1);
                        bdt_ada_errors.push_back(-100.);
                        bdt_grad_errors.push_back(-100.);
                        vbf_values.push_back(-100.);
                        weights.push_back(0.);
                    }

                } else { // In case CiC selection fails now
                    categories.push_back(-1);
                    bdt_ada_errors.push_back(-100.);
                    bdt_grad_errors.push_back(-100.);
                    vbf_values.push_back(-100.);
                    weights.push_back(0.);
                }

            }// end loop on systematics steps

            // Fill In the Corect Systematic Set ---------------------------------------------------------------------------//
            if (includeVBF) l.rooContainer->InputSystematicSet("sig_VBF_" +currentTypeSignalLabel,(*si)->name(),categories,vbf_values,weights);
            l.rooContainer->InputSystematicSet("sig_BDT_ada_" +currentTypeSignalLabel,(*si)->name(),categories,bdt_ada_errors,weights);
            l.rooContainer->InputSystematicSet("sig_BDT_grad_"+currentTypeSignalLabel,(*si)->name(),categories,bdt_grad_errors,weights);
            // -------------------------------------------------------------------------------------------------------------//

        } // Close on Smearers
    } // End Signal Systematic Study 
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------
}

// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
void MvaAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{
    vtxAna_.setBranchAdresses(t,"vtx_std_");
    vtxAna_.getBranches(t,"vtx_std_",s);
}

// ----------------------------------------------------------------------------------------------------
bool MvaAnalysis::SelectEvents(LoopAll& l, int jentry) 
{
    return true;
}
// ----------------------------------------------------------------------------------------------------

int MvaAnalysis::SignalType(int cur_type){
    int i0 = -1;
    if (doTraining) {
        if (cur_type == -53 || cur_type == -54 ||  cur_type == -55 || cur_type == -56){//121 
            i0 = 0;}
        else if (cur_type == -57 || cur_type == -58 ||  cur_type == -59 || cur_type == -60){//123 
            i0 = 1;}
    } else{
        if (cur_type == -13 || cur_type == -14 ||  cur_type == -15 || cur_type == -16){//105 
            i0 = 0;}
        else if (cur_type == -17 || cur_type == -18 ||  cur_type == -19 || cur_type == -20){//110 
            i0 = 1;}
        else if (cur_type == -21 || cur_type == -22 ||  cur_type == -23 || cur_type == -24){//115 
            i0 = 2;}
        else if (cur_type == -25 || cur_type == -26 ||  cur_type == -27 || cur_type == -28){//120 
            i0 = 3;}
        else if (cur_type == -37 || cur_type == -38 ||  cur_type == -39 || cur_type == -40){//125 
            i0 = 4;}
        else if (cur_type == -29 || cur_type == -30 ||  cur_type == -31 || cur_type == -32){//130 
            i0 = 5;}
        else if (cur_type == -41 || cur_type == -42 ||  cur_type == -43 || cur_type == -44){//135 
            i0 = 6;}
        else if (cur_type == -33 || cur_type == -34 ||  cur_type == -35 || cur_type == -36){//140 
            i0 = 7;}
        //    else if (cur_type == -45 || cur_type == -46 ||  cur_type == -47 || cur_type == -48){//145 
        //        i0 = 8;}
        else if (cur_type == -49 || cur_type == -50 ||  cur_type == -51 || cur_type == -52){//150 
            i0 = 8;}
    }
    return i0;
}

void MvaAnalysis::FillSignalLabelMap(){

    // Basically A Map of the ID (type) to the signal's name which can be filled Now:
    signalLabels[-57]="ggh_123.0";
    signalLabels[-58]="vbf_123.0";
    signalLabels[-60]="wzh_123.0";
    signalLabels[-59]="tth_123.0";
    signalLabels[-53]="ggh_121.0";
    signalLabels[-54]="vbf_121.0";
    signalLabels[-56]="wzh_121.0";
    signalLabels[-55]="tth_121.0";
    signalLabels[-65]="ggh_160.0";
    signalLabels[-66]="vbf_160.0";
    signalLabels[-68]="wzh_160.0";
    signalLabels[-67]="tth_160.0";
    signalLabels[-61]="ggh_155.0";
    signalLabels[-62]="vbf_155.0";
    signalLabels[-64]="wzh_155.0";
    signalLabels[-63]="tth_155.0";
    signalLabels[-49]="ggh_150.0";
    signalLabels[-50]="vbf_150.0";
    signalLabels[-52]="wzh_150.0";
    signalLabels[-51]="tth_150.0";
    signalLabels[-45]="ggh_145.0";
    signalLabels[-46]="vbf_145.0";
    signalLabels[-48]="wzh_145.0";
    signalLabels[-47]="tth_145.0";
    signalLabels[-33]="ggh_140.0";
    signalLabels[-34]="vbf_140.0";
    signalLabels[-36]="wzh_140.0";
    signalLabels[-35]="tth_140.0";
    signalLabels[-41]="ggh_135.0";
    signalLabels[-42]="vbf_135.0";
    signalLabels[-44]="wzh_135.0";
    signalLabels[-43]="tth_135.0";
    signalLabels[-29]="ggh_130.0";
    signalLabels[-30]="vbf_130.0";
    signalLabels[-32]="wzh_130.0";
    signalLabels[-31]="tth_130.0";
    signalLabels[-37]="ggh_125.0";
    signalLabels[-38]="vbf_125.0";
    signalLabels[-40]="wzh_125.0";
    signalLabels[-39]="tth_125.0";
    signalLabels[-25]="ggh_120.0";
    signalLabels[-26]="vbf_120.0";
    signalLabels[-28]="wzh_120.0";
    signalLabels[-27]="tth_120.0";
    signalLabels[-21]="ggh_115.0";
    signalLabels[-22]="vbf_115.0";
    signalLabels[-24]="wzh_115.0";
    signalLabels[-23]="tth_115.0";
    signalLabels[-17]="ggh_110.0";
    signalLabels[-18]="vbf_110.0";
    signalLabels[-19]="wzh_110.0";
    signalLabels[-20]="tth_110.0";
    signalLabels[-13]="ggh_105.0";
    signalLabels[-14]="vbf_105.0";
    signalLabels[-16]="wzh_105.0";
    signalLabels[-15]="tth_105.0";
    signalLabels[-69]="ggh_100.0";
    signalLabels[-70]="vbf_100.0";
    signalLabels[-72]="wzh_100.0";
    signalLabels[-71]="tth_100.0";
}

std::string MvaAnalysis::GetSignalLabel(int id){

    // For the lazy man, can return a memeber of the map rather than doing it yourself
    std::map<int,std::string>::iterator it = signalLabels.find(id);

    if (it!=signalLabels.end()){
        return it->second;

    } else { 

        std::cerr << "No Signal Type defined in map with id - " << id << std::endl;
        return "NULL";
    }

}

void MvaAnalysis::SetBDTInputVariables(TLorentzVector *lead_p4, TLorentzVector *sublead_p4, double lead_r9, double sublead_r9, MassResolution *massResolutionCalculator, double vtx_mva, double mass_hypothesis, double bdtoutput, double evweight, int cat){

    TLorentzVector Higgs = *lead_p4 + *sublead_p4;   
    float mass    = Higgs.M();

    _log_H_pt =  log10( Higgs.Pt());
    _H_eta = fabs(Higgs.Eta());
    _d_phi = fabs(lead_p4->DeltaPhi(*sublead_p4));
    _cos_d_phi = TMath::Cos(lead_p4->Phi()-sublead_p4->Phi());        
    _max_eta = max(fabs(lead_p4->Eta()),fabs(sublead_p4->Eta()));
    _min_r9  = min(lead_r9,sublead_r9);
    _pho1_eta = lead_p4->Eta();
    _pho2_eta = sublead_p4->Eta();

    _mgg = mass;
    _pho1_phi = lead_p4->Phi();
    _pho1_pt = lead_p4->Pt();
    _pho1_r9 = lead_r9;

    _pho2_phi = sublead_p4->Phi();
    _pho2_pt = sublead_p4->Pt();
    _pho2_r9 = sublead_r9;

    _H_pt = Higgs.Pt();
    _Ht = lead_p4->Pt()+sublead_p4->Pt();

    _d_eta = lead_p4->Eta()-sublead_p4->Eta();
    _mod_d_eta = fabs(lead_p4->Eta()-sublead_p4->Eta());
    _cos_theta_star = fabs(lead_p4->E()-sublead_p4->E())/Higgs.P();


    _vtx_prob = 1.-0.49*(vtx_mva+1.0);

    _wt= evweight;

    _pho1_ptOverM = lead_p4->Pt()/mass_hypothesis;
    _pho2_ptOverM = sublead_p4->Pt()/mass_hypothesis;
    _deltaMOverM = (mass-mass_hypothesis)/mass_hypothesis;

    double massResolution = massResolutionCalculator->massResolutionCorrVtx();
    _deltaMOverSigmaM = (mass-mass_hypothesis)/massResolution;
    _sigmaMOverM = massResolution/mass;
    _sigmaMOverM_wrongVtx = massResolutionCalculator->massResolutionWrongVtx()/mass;
    _H_ptOverM    = Higgs.Pt()/mass_hypothesis;
    _cat        = cat;
    _bdtoutput = bdtoutput;
}

void MvaAnalysis::SetBDTInputTree(TTree *tree){
    mvaFile_->cd();
    TBranch *b_log_H_pt             = tree->Branch("log_H_pt", &_log_H_pt , "log_H_pt/F");;
    TBranch *b_H_ptOverM            = tree->Branch("H_ptOverM", &_H_ptOverM , "_H_ptOverM/F");;
    TBranch *b_H_eta                = tree->Branch("H_eta", &_H_eta , "H_eta/F");;
    TBranch *b_d_phi                = tree->Branch("d_phi", &_d_phi,"d_phi/F");
    TBranch *b_max_eta              = tree->Branch("max_eta", &_max_eta , "max_eta/F");;
    TBranch *b_min_r9               = tree->Branch("min_r9", &_min_r9 , "min_r9/F");;
    TBranch *b_pho1_eta             = tree->Branch("pho1_eta", &_pho1_eta , "pho1_eta/F");
    TBranch *b_pho2_eta             = tree->Branch("pho2_eta", &_pho2_eta , "pho2_eta/F");
    TBranch *b_pho1_ptOverM         = tree->Branch("pho1_ptOverM", &_pho1_ptOverM , "pho1_ptOverM/F");;
    TBranch *b_pho2_ptOverM         = tree->Branch("pho2_ptOverM", &_pho2_ptOverM , "pho2_ptOverM/F");;
    TBranch *b_deltaMOverM          = tree->Branch("deltaMOverM", &_deltaMOverM,"deltaMOverM/F");
    TBranch *b_deltaMOverSigmaM     = tree->Branch("deltaMOverSigmaM", &_deltaMOverSigmaM,"deltaMOverSigmaM/F");
    TBranch *b_sigmaMOverM          = tree->Branch("sigmaMOverM", &_sigmaMOverM,"sigmaMOverM/F");
    TBranch *b_sigmaMOverM_wrongVtx = tree->Branch("sigmaMOverM_wrongVtx", &_sigmaMOverM_wrongVtx,"sigmaMOverM_wrongVtx/F");
    TBranch *b_mgg                  = tree->Branch("mgg", &_mgg, "mgg/F");
    TBranch *b_pho1_phi             = tree->Branch("pho1_phi", &_pho1_phi , "pho1_phi/F");
    TBranch *b_pho1_pt              = tree->Branch("pho1_pt", &_pho1_pt , "pho1_pt/F");;
    TBranch *b_pho1_r9              = tree->Branch("pho1_r9", &_pho1_r9 , "_pho1_r9/F");;
    TBranch *b_pho2_phi             = tree->Branch("pho2_phi", &_pho2_phi , "pho2_phi/F");
    TBranch *b_pho2_pt              = tree->Branch("pho2_pt", &_pho2_pt , "pho2_pt/F");;
    TBranch *b_pho2_r9              = tree->Branch("pho2_r9", &_pho2_r9 , "_pho2_r9/F");;
    TBranch *b_H_pt                 = tree->Branch("H_pt", &_H_pt , "H_pt/F");;
    TBranch *b_Ht                   = tree->Branch("Ht", &_Ht , "Ht/F");;
    TBranch *b_d_eta                = tree->Branch("d_eta", &_d_eta,"d_eta/F");
    TBranch *b_mod_d_eta            = tree->Branch("mod_d_eta", &_mod_d_eta,"mod_d_eta/F");
    TBranch *b_cos_theta_star       = tree->Branch("cos_theta_star", &_cos_theta_star , "cos_theta_star/F");;
    TBranch *b_vtx_prob             = tree->Branch("vtx_prob", &_vtx_prob, "vtx_prob/F");
    TBranch *b_wt                   = tree->Branch("wt", &_wt, "wt/F");
    TBranch *b_cat                  = tree->Branch("cat", &_cat, "cat/I");
    TBranch *b_sideband             = tree->Branch("sideband", &_sideband, "sideband/I");
    TBranch *b_bdtoutput            = tree->Branch("bdtoutput", &_bdtoutput, "bdtoutput/F");
}

void MvaAnalysis::ResetAnalysis(){
    // Reset Random Variable on the EnergyResolution Smearer
    eResolSmearer->resetRandom();
}


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
