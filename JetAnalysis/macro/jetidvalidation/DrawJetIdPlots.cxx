#include <algorithm>

void DrawJetIdPlots(string etaRange  = "TK", 
		    string ptRange   = "pt20to30",
		    string inputdir  = "rootfiles/V00-00-08",
		    string outputdir = "./")
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  gStyle->SetOptStat(1110);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleBorderSize(0);
 
  string infile1 = inputdir+"/testjetid_mc_"+etaRange+"_"+ptRange+".root";
  string infile2 = inputdir+"/testjetid_data_"+etaRange+"_"+ptRange+".root";
  string dirName = outputdir+"/"+etaRange+"_"+ptRange;
  string fileType = "png";

  
  gSystem->mkdir(dirName.c_str(),true);
  
  cout << "Producing validation plots for: " << infile1 << " and " << infile2 << endl;

  
  TFile* f[2]; 
  f[0] = new TFile(infile1.c_str()); //MC
  f[1] = new TFile(infile2.c_str()); //DATA

   

  TH1F *hjetPt[2];
  for (int i=0;i<2;i++)
    hjetPt[i] = (TH1F*)f[i]->Get("hjetPt") ; 
  
  float nEvents0 = hjetPt[0]->GetSumOfWeights();
  float nEvents1 = hjetPt[1]->GetSumOfWeights();
  float s        = nEvents1/nEvents0;


  // Define list of object names
  const int nObj=16;
  string objName[nObj]={"hNvtx",
			"hjetPt",
			"hjetEta",
// 			"hnCharged",
// 			"hnNeutrals",
			"hdR2Mean",
			"hdR2Mean",
			"hdR2Mean",
			"hdR2Mean",
			"hfrac01",
			"hfrac02",
			"hfrac03",
			"hfrac04",
			"hfrac05",
			"hbeta",
			"hbetaStar",
			"hsimpleDiscriminant",
			"hfullDiscriminant"
  };

  string objNamePU[nObj];
  string objNameNoPU[nObj];
  for (int i = 0; i < nObj; i++){
    objNamePU[i]   = objName[i]+"_PU" ;
    objNameNoPU[i] = objName[i]+"_NoPU" ;
//     cout << objName[i]<< endl;
//     cout << objNamePU[i]<< endl;
//     cout << objNameNoPU[i]<< endl;
  }
 
  
  string objTitle[nObj]={"Nvtx       ",
			 "jetPt      ",
			 "jetEta     ",
			 "nCharged   ",
			 "nNeutrals  ",
			 "dR2Mean    ",
			 "dRMean     ",
			 "hfrac01    ",
			 "hfrac02    ",
			 "hfrac03    ",
			 "hfrac04    ",
			 "hfrac05    ",
			 "beta       ",
			 "betaStar   ",
			 "simpleDiscriminant",
			 "fullDiscriminant"
  };

  char *labelX[nObj]={"number of vertices",
		      "jet p^{T} (GeV)",
		      "jet #eta",
		      "nCharged   ",
		      "nNeutrals  ",
		      "dR2Mean    ",
		      "dRMean     ",
		      "frac01     ",
		      "frac02     ",
		      "frac03     ",
		      "frac04     ",
		      "frac05     ",
		      "beta       ",
		      "betaStar   ",
 		      "mva        ",
 		      "mva        "
   };
  
  char *labelY[1]={"Number of entries"};

  double xMin[nObj]={0,   // nvtx
		     0.,  // jet pt
		     -5., // jet eta 
		     0, 0,// particles multiplicity
		     0, 0,// dR
		     0,   // frac01
		     0,   // frac02
		     0,   // frac03
		     0,   // frac04
		     0,   // frac05
		     0, 0,// beta 
		     -1.,-1. // mva
  };
  
  double xMax[nObj]={40,   // nvtx
		     100., // jet pt
		     5.,   // jet eta 
		     50,   // particles multiplicity
		     0.5, 0.5, // dR
		     1,    // frac01
		     1,    // frac02
		     1,    // frac03
		     1,    // frac04
		     1,    // frac05
		     1, 1, // beta 
		     1.,1. //mva
  };
  
  int reBin[nObj]  = {1,   // nvtx
		      1,   // jet pt
		      1,   // jet eta 
		      1,   // particles multiplicity
		      2, 2,// dR
		      2,   // frac01
		      2,   // frac02
		      2,   // frac03
		      2,   // frac04
		      2,   // frac05
		      5, 5,// beta 
		      2, 2 //mva
  };

  int optLogY[nObj] = {0,   // nvtx
		       0,   // jet pt
		       0,   // jet eta 
		       0,   // particles multiplicity
		       0, 0,// dR
		       0,   // frac01
		       0,   // frac02
		       0,   // frac03
		       0,   // frac04
		       0,   // frac05
		       1, 1,// beta 
		       0,0  //mva
  };
  
  TH1D* h[2][100]; 
  TH1D* hPU[2][100]; 
  TH1D* hNoPU[2][100]; 
  TCanvas *c[100];
  TPaveStats *st[100];
  TLegend *legend[100];
 
  int iHisto = 0;
  while(iHisto<nObj){
    for (int ifile=0;ifile<2;ifile++){ 

      //cout << ifile<< "  " << iHisto << "  " << objName[iHisto] << endl;
      h[ifile][iHisto] = (TH1D*)f[ifile]->Get(objName[iHisto].c_str());
      h[ifile][iHisto]->Rebin(reBin[iHisto]);

      if (iHisto > 0){
	hPU[ifile][iHisto] = (TH1D*)f[ifile]->Get(objNamePU[iHisto].c_str());
	hNoPU[ifile][iHisto] = (TH1D*)f[ifile]->Get(objNameNoPU[iHisto].c_str());
	hPU[ifile][iHisto]->Rebin(reBin[iHisto]);
	hNoPU[ifile][iHisto]->Rebin(reBin[iHisto]);
      }

     
      if (ifile == 0) {
	// open a new canvas
	c[iHisto] = new TCanvas(objName[iHisto].c_str(),objName[iHisto].c_str(),50+iHisto*20,50+iHisto*5,500,400);
	c[iHisto]->cd();
	// customize and plot
	h[ifile][iHisto]->GetYaxis()->SetTitleOffset(1.4);
	h[ifile][iHisto]->GetYaxis()->SetTitle(labelY[0]);
	h[ifile][iHisto]->GetXaxis()->SetTitle(labelX[iHisto]);
	h[ifile][iHisto]->SetFillColor(kGray);
	h[ifile][iHisto]->SetFillStyle(1001);
	h[ifile][iHisto]->SetTitle(objTitle[iHisto].c_str());
	h[ifile][iHisto]->Scale(s);
	h[ifile][iHisto]->GetXaxis()->SetRangeUser(xMin[iHisto],xMax[iHisto]);
	h[ifile][iHisto]->Draw();

	if (iHisto>0){
	  hPU[ifile][iHisto]->SetLineColor(kRed);
	  hPU[ifile][iHisto]->SetLineStyle(2);
	  hPU[ifile][iHisto]->SetLineWidth(2);
	  hPU[ifile][iHisto]->Scale(s);
	  hPU[ifile][iHisto]->Draw("same");
	  
	  hNoPU[ifile][iHisto]->SetLineColor(kBlue);
	  hNoPU[ifile][iHisto]->SetLineStyle(2);
	  hNoPU[ifile][iHisto]->SetLineWidth(2);
	  hNoPU[ifile][iHisto]->Scale(s);
	  hNoPU[ifile][iHisto]->Draw("same");
	}
	
      }
      if (ifile == 1) {
	
	h[ifile][iHisto]->SetMarkerStyle(20);
	h[ifile][iHisto]->SetMarkerSize(0.7);
	h[ifile][iHisto]->SetMarkerColor(kBlack);
	h[ifile][iHisto]->Draw("epsames");
	
	// set the log-scale and range 
	float maxy = max (h[ifile][iHisto]->GetMaximum(),h[ifile-1][iHisto]->GetMaximum() );
	float miny = h[ifile][iHisto]->GetMinimum();
	if (optLogY[iHisto]) {
	  c[iHisto]->SetLogy();
	  h[0][iHisto]->SetMaximum(maxy*100);
	  h[0][iHisto]->SetMinimum(0.1);
	}
	else  h[0][iHisto]->SetMaximum(maxy*1.4);
	
	c[iHisto]->Update();
	
	// stat. box
	st[iHisto]= (TPaveStats*)(h[ifile][iHisto]->GetListOfFunctions()->FindObject("stats"));
	st[iHisto]->SetY1NDC(0.72); //new x start position
	st[iHisto]->SetY2NDC(0.85); //new x end position
	st[iHisto]->SetTextColor(kGray+1);
	st[iHisto]->Draw();
	
	//legend 
	legend[iHisto] = new TLegend(0.5, 0.65, 0.85, 0.85);
	legend[iHisto]->SetFillStyle(0);
	legend[iHisto]->SetBorderSize(0);
	legend[iHisto] ->AddEntry(h[1][0], "DATA","LP");
	legend[iHisto] ->AddEntry(h[0][0], "MC","FL");
	if (iHisto>0){
	  legend[iHisto] ->AddEntry(hPU[0][1], "MC - PU","LP");
	  legend[iHisto] ->AddEntry(hNoPU[0][1], "MC - no PU","LP");
	}
	legend[iHisto] ->Draw("same");	
      }
    }
    

    string myname =  (dirName+"/"+objName[iHisto]+"."+fileType).c_str();
    c[iHisto]->Print(myname.c_str(),fileType.c_str());
    
    
    iHisto++;
    
    
  }
  

  
 
}

