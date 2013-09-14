#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TArrow.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TColor.h"
#include "TROOT.h"
#include "Math/DistFunc.h"

void setTDRStyle();

int extractSignificanceStats(
			     bool unblind=true,
			     TString legALT="2^{+}_{m}(gg)",
			     TString nameALT="2pmgg",
			     TString input="qmu_*.root",
           int rebin=500,
           float xlow=-10.,
           float xhigh=10.
			     ){

	float labelCLs=0;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  setTDRStyle();

  TChain* t = new TChain("q");
  t->Add(input);

  float q,m,w;
  int type;
  t->SetBranchAddress("q",&q);
  t->SetBranchAddress("mh",&m);
  t->SetBranchAddress("weight",&w);
  t->SetBranchAddress("type",&type);

  unsigned long ntoysALT = t->Draw("","type>0","goff");
  unsigned long ntoysSM  = t->Draw("","type<0","goff");

  TH1F *hSM=new TH1F("hSM;S = -2 #times ln(L_{ALT}/L_{SM});Number of Toys","",32000,xlow,xhigh);
  TH1F *hALT=new TH1F("hALT;S = -2 #times ln(L_{ALT}/L_{SM});Number of Toys","",32000,xlow,xhigh);
  TH1F *hObs=new TH1F("hObserved","",32000,xlow,xhigh);

  std::vector<float> vSM, vALT,vObs;
  int counterData=0;
  for(int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);

    if(i==0) cout<<"MASS in the TREE = "<<m<<endl<<endl;

    q *= -2.0;
    if(type<0){ //SM hypothesis 
      hSM->Fill(q);
      vSM.push_back(q);
    }
    else if(type>0){//ALT hypothesis
      hALT->Fill(q);
      vALT.push_back(q);
    }
    else{
      counterData++;
      hObs->Fill(q);
      //cout<<"DATA -q ="<<q<<endl;
      vObs.push_back(q);
    }
  }//end loop on tree entries
  //cout<<"Finished to loop, sorting vectors "<<vSM.size()<<" "<<vALT.size()<<" "<<vObs.size()<<endl;
  sort(vSM.begin(),vSM.end());//sort in ascending order
  sort(vALT.begin(),vALT.end()); 
  sort(vObs.begin(),vObs.end());

  if(vALT.at(0)>vSM.at(ntoysSM-1)){
    cout<<"Swapped distributions !!! The alternative model should stay on the negative side of the significance."<<endl;
    cout<<"Please edit the code and change the sign of q when filling histos and vectors in the loop on tree entries"<<endl;
    return 1;
  }

  if((vSM.size()!= ntoysSM)||(vALT.size()!= ntoysALT)){
    cout<<"Mismatch in size of vectors and #entries of histograms ! vSM.size()="<< vSM.size() <<"  ntoysSM="<<ntoysSM<<endl;
    cout<<"Mismatch in size of vectors and #entries of histograms ! vALT.size()="<< vALT.size() <<"  ntoysALT="<<ntoysALT<<endl;
    return 1;
  }

  float medianSM=vSM.at(int(ntoysSM/2));
  float medianALT=vALT.at(int(ntoysALT/2));
  cout<<"Toys generated SM: "<<ntoysSM<<"\t"<<ntoysALT<<endl;
  cout<<"Mean of SM/ALT hypothesis: "<<hSM->GetMean()<<"\t"<<hALT->GetMean()<<endl;
  cout<<"RMS  of SM/ALT hypothesis: "<<hSM->GetRMS()<<"\t"<<hALT->GetRMS()<<endl;
  cout<<"Median of SM/ALT hypothesis: "<<medianSM<<"\t"<<medianALT<<endl;

  const float step=0.05;
  float coverage=0.0;
  float diff=10.0;
  float cut=vALT.at(0)-step;
  float crosspoint=-99.0;
  int startSM=ntoysSM-1, startALT=0;
  //cout<<"Starting to loop with cut at "<<cut<<endl;

  while(cut<=vSM.at(ntoysSM-1)+step){
    float cutSM=-1.0,cutALT=-1.0;

    for(long iSM=startSM;iSM>=0;iSM--){      
      if(vSM.at(iSM)<cut){//gotcha
	cutSM=ntoysSM-iSM;
	break;
      }
    }

    for(unsigned long iALT=startALT;iALT<ntoysALT;iALT++){
      if(vALT.at(iALT)>cut){//gotcha
	cutALT=iALT;
	break;
      }
    }

    if(cutSM>=0 && cutALT>=0){
      float fracSM=(ntoysSM-cutSM)/ntoysSM;
      float fracALT=(ntoysALT-cutALT)/ntoysALT;
      if(fabs(fracSM-fracALT)<diff){
	diff=fabs(fracSM-fracALT);
	coverage=fabs(fracSM+fracALT)/2.0;
	crosspoint=cut;
      }
    }
    cut+=step;
  }//end while loop
 
  //cout<<"Finished loop on vector elements, min diff is "<<diff<<", looped until cut_fin="<<cut<<endl;
  cout<<"q value where SM and ALT distributions have same area on opposite sides: "<<crosspoint<<"  Coverage="<<coverage<< "  Separation="<<1.0-coverage <<endl;
  float separation=2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from tail prob: "<<separation<<endl<<endl<<endl;
  

  float integralSM=hSM->Integral();
  float integralALT=hALT->Integral();
 
  float tailSM=hSM->Integral(1,hSM->FindBin(medianALT))/integralSM;
  float tailALT=hALT->Integral(hALT->FindBin(medianSM),hALT->GetNbinsX())/integralALT;
  cout<<"Median point prob SM : "<<tailSM <<"  ("<<ROOT::Math::normal_quantile_c(tailSM, 1.0) <<" sigma)"<<endl;
  cout<<"Median point prob ALT: "<<tailALT<<"  ("<<ROOT::Math::normal_quantile_c(tailALT,1.0) <<" sigma)"<<"\t Exp 1-CLS: "<< 1.0-tailALT/0.5 << endl;
  labelCLs=tailALT/0.5;

  diff=10.0;
  coverage=0.0;
  for(int i=1;i<hSM->GetNbinsX();i++){
    
    float fracSM=hSM->Integral(1,i) / integralSM;
    float fracALT=hALT->Integral(i,hALT->GetNbinsX()) / integralALT;
    if(fabs(fracSM-fracALT)<diff){
      diff=fabs(fracSM-fracALT);
      coverage=(fracSM+fracALT)/2.0;
    }

  }

  float sepH= 2*ROOT::Math::normal_quantile_c(1.0 - coverage, 1.0);
  cout<<"Separation from histograms = "<<sepH<<" with coverage "<<coverage<<endl;
  double observedQ=vObs.at(0);

  if(unblind){
    if ( (*(vObs.end()-1) / *vObs.begin() - 1.0) > 0.1/100 ){
	cout<<"Ooops ! The size of the vector with the observed separation is "<<vObs.size()<<" !"<<endl;
	cout << "Found more than 0.1% difference between the first and the last: " << *vObs.begin() << " -- " << *(vObs.end()-1) << endl;
	return -1;
    }

    else{
      float obsTailSM = hSM->Integral(1,hSM->FindBin(observedQ))/integralSM;
      float obsTailALT=hALT->Integral(hALT->FindBin(observedQ),hALT->GetNbinsX())/integralALT;
      cout<<"P(SM < Obs) : "<<obsTailSM <<"  ("<<ROOT::Math::normal_quantile_c(obsTailSM,1.0) <<" sigma)"<<endl;
      cout<<"P(ALT > Obs): "<<obsTailALT<<"  ("<<ROOT::Math::normal_quantile_c(obsTailALT,1.0) <<" sigma)"<<endl;

      float obsCLsRatio = obsTailALT / (1.0 - obsTailSM);
      cout<<"CLs criterion P(ALT > Obs) / P(SM > Obs) : "<<obsCLsRatio<<"  ("<<ROOT::Math::normal_quantile_c(obsCLsRatio,1.0) <<" sigma)"<< "\t Obs 1-CLS="<<1.0-obsCLsRatio <<endl;
      labelCLs=obsCLsRatio;
    }

    cout << "\n\nOBSERVED SIGNIFICANCE" << endl;
    cout << "observation: " << vObs[0] << endl;
    cout << "bin: " << hObs->GetMaximumBin() << endl;
    cout << " --------------- " << endl;
    double obsPvalSM = 1-hSM->Integral(0,hObs->GetMaximumBin())/integralSM;
    cout << "pvalue SM : " << obsPvalSM << endl;
    cout << "signif SM : " << ROOT::Math::normal_quantile_c(obsPvalSM,1.0) << endl;
    double obsPvalALT =  hALT->Integral(0,hObs->GetMaximumBin())/integralALT;
    cout << "pvalue ALT: " << obsPvalALT << endl;
    cout << "signif ALT: " << ROOT::Math::normal_quantile_c(obsPvalALT,1.0) << endl<<endl<<endl;

   }//end if unblinding

  //Plotting
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c","c",500,500);
  c1->cd();
  //hSM->Rebin(20e7/ntoysSM);
  //hALT->Rebin(20e7/ntoysALT);
  hSM->Rebin(rebin);
  hALT->Rebin(rebin);
  TString xtitle=" -2 #times ln(L_{"+legALT+"} / L_{0^{+}})";
  hSM->SetXTitle(xtitle);
  hSM->SetYTitle("Probability density");
  hALT->SetXTitle(xtitle);
  hALT->SetYTitle("Probability density");
  hSM->SetLineColor(kRed+2);
  hSM->SetLineStyle(2);
  hSM->SetFillColor(798);
  hSM->SetLineWidth(2);
  hALT->SetFillColor(kAzure+7);
  hALT->SetLineColor(kBlue);
  hALT->SetLineWidth(1);
  hALT->SetFillStyle(3001);

  hObs->SetLineColor(kRed);
  hObs->SetLineWidth(2);
  
  hSM->GetXaxis()->SetRangeUser(-30.0,30.0);
  hSM->GetXaxis()->SetLabelFont(42);
  hSM->GetXaxis()->SetLabelOffset(0.007);
  hSM->GetXaxis()->SetLabelSize(0.045);
  hSM->GetXaxis()->SetTitleSize(0.05);
  hSM->GetXaxis()->SetTitleOffset(1.15);
  hSM->GetXaxis()->SetTitleFont(42);
  hSM->GetYaxis()->SetLabelFont(42);
  hSM->GetYaxis()->SetLabelOffset(0.007);
  hSM->GetYaxis()->SetLabelSize(0.045);
  hSM->GetYaxis()->SetTitleSize(0.05);
  hSM->GetYaxis()->SetTitleOffset(1.8);
  hSM->GetYaxis()->SetTitleFont(42); 
  //TGaxis::SetMaxDigits(2); 
  hSM->Scale(1./hSM->Integral("width"));
  hALT->Scale(1./hALT->Integral("width"));
  float maxhSM=hSM->GetBinContent(hSM->GetMaximumBin());
  float maxhALT=hALT->GetBinContent(hALT->GetMaximumBin());
  if(maxhALT>maxhSM){
    hSM->SetMaximum(maxhALT*1.3);
    hALT->SetMaximum(maxhALT*1.3);
  }
  else{
    hSM->SetMaximum(maxhSM*1.3);
    hALT->SetMaximum(maxhSM*1.3);
  }

  TH1F* hSM2= (TH1F*)(hSM->Clone("hSM2")); 
  hSM2->SetFillStyle(0);
  hSM2->SetFillColor(0);
  hSM2->SetLineColor(kRed+2);
  hSM2->SetLineStyle(2);
  hSM2->SetLineWidth(2);

  hSM->Draw("HIST");
  hALT->Draw("HISTsames");
  hSM2->Draw("HISTsames");

  TArrow *obsArrow=0;
  if(unblind){
    obsArrow=new TArrow(observedQ,hSM->GetMaximum()/5.0,observedQ,0.0,0.05,"|->");
    obsArrow->SetLineColor(kRed);
    obsArrow->SetLineWidth(2.0);
    obsArrow->Draw();
  }

  TLegend *leg = new TLegend(0.63,0.73,0.88,0.93);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(hSM, "0^{+}","f");
  leg->AddEntry(hALT, legALT,"f");
  if(unblind){
    leg->AddEntry(hObs,"CMS data","L");
    leg->AddEntry((TObject*)0,Form("(CL_{s}^{obs.} = %2.1f%%)",labelCLs*100),"");
  }else{
	  leg->AddEntry((TObject*)0,"","");
	  leg->AddEntry((TObject*)0,Form("(CL_{s}^{exp.} = %2.1f%%)",labelCLs*100),"");
  }

  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.18);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);

  TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  pt->AddText(0.01,0.5,"CMS preliminary");
  //text = pt->AddText(0.3,0.6,"#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}  #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  pt->AddText(0.6,0.6,"#sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  pt->Draw();   
  
  
  leg->Draw();

  c1->SaveAs(nameALT+".png");
  c1->SaveAs(nameALT+".pdf");
  c1->SaveAs(nameALT+".eps");
  c1->SaveAs(nameALT+".root");
 
 return 0;
}//end main

void setTDRStyle() {
  //TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(750); //Height of canvas
  gStyle->SetCanvasDefW(1050); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

// For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

// For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

// For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  //  gStyle->SetErrorMarker(20);
  gStyle->SetErrorX(0.);
  
  gStyle->SetMarkerStyle(20);

//For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

//For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.010);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.25);
  gStyle->SetStatW(0.15);
  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

// Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.04);

// For the Global title:

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.005);
  // gStyle->SetTitleH(0); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

// For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  //gStyle->SetNdivisions(505, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

// Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

// Postscript options:
  gStyle->SetPaperSize(20.,20.);
  // gStyle->SetLineScalePS(Float_t scale = 3);
  // gStyle->SetLineStyleString(Int_t i, const char* text);
  // gStyle->SetHeaderPS(const char* header);
  // gStyle->SetTitlePS(const char* pstitle);

  // gStyle->SetBarOffset(Float_t baroff = 0.5);
  // gStyle->SetBarWidth(Float_t barwidth = 0.5);
  // gStyle->SetPaintTextFormat(const char* format = "g");
  // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // gStyle->SetTimeOffset(Double_t toffset);
  // gStyle->SetHistMinimumZero(kTRUE);

  Int_t n=35;
  Int_t *colors = new Int_t[n];
  for (Int_t i =0; i < n; i++) {
    //colors[i] = i+61;
    //colors[i] = 70-i;
    colors[i] = i+63;
    }
  /*colors[0]=61
  colors[1]=62;
  colors[2]=63;
  colors[3]=;
  colors[4]=;
  colors[5]=;
  colors[6]=;
  colors[7]=;
  colors[8]=;
  colors[9]=;
  colors[10]=;*/
  gStyle->SetPalette(n, colors);
  
  //gStyle->SetPalette(1, 0);
  //gStyle->cd();
}
/*
void damnedPlot(){
TCanvas *c1 =  new TCanvas("paddy","paddy",700,700)
  c1->Divide(1,2,0,0);
 c1->cd(1).SetBottomMargin(0.001);
 c1->cd(1).SetTopMargin(0.01);
 c1->cd(1).SetRightMargin(0.01);
 c1.cd(2);
 c1.cd(2).SetTopMargin(0.001);
 c1.cd(2).SetRightMargin(0.01);

}
*/
