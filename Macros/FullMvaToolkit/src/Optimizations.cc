#include "../interface/Optimizations.h"
#include "../interface/LocalRegression.h"

using namespace std;

Optimizations::Optimizations(TH2F *signalHisto, TH2F *backgroundHisto){

	// set some defaults 
	maxNumberOfBins = 8;
        nNewBins = 100;
	predefMin = 20;	
	
	nFinalBins = 0;
	delta = 0.00001;
	threshold = 0.001; // % imporovement on significance for extra bin
	
	targetS2d = (TH2F*)signalHisto->Clone();
	targetB2d = (TH2F*)backgroundHisto->Clone();

	n2dbinsX = signalHisto->GetNbinsX();
  	n2dbinsY = signalHisto->GetNbinsY();

}

void Optimizations::smoothHistograms(double bandwidths,double bandwidthb, int mode){

	if (mode == 2) std::cout << "Running Polynomial Regression with TMVA! (slows down the result)" <<std::endl;
	std::cout << "Smoothing Signal Histogram with Local Regression" << std::endl;
	smoothedHistogram(targetS2d,bandwidths,mode);
	std::cout << "Smoothing Background Histogram with Local Regression" << std::endl;
	smoothedHistogram(targetB2d,bandwidthb,mode);
} 

void Optimizations::runOptimization(){
	
	// Make by vectors before running optimization
	signalVector1=new double[maxNumberOfBins+2];
	backgroundVector1=new double[maxNumberOfBins+2];

	TH2F *hsoverb = (TH2F*)targetS2d->Clone();
	hsoverb->Divide(targetB2d);

	double maximumSoverB = findMaximum(hsoverb);

 	// Now we make it a Log(S/B) plot 
  	for (int k=1;k<=n2dbinsX;k++){
    		for (int l=1;l<=n2dbinsY;l++){
		  hsoverb->SetBinContent(k,l,defx(hsoverb->GetBinContent(k,l)));
    		}
  	}

	SBscale =  1./(defx(maximumSoverB));
	hsoverb->Scale(SBscale);

	FinalHist = (TH2F*) hsoverb->Clone();
	FinalHist->SetName("Category_Map");

	// Create the Integral histograms

	targetS = new TH1F("targetS","",nNewBins,0,1);
 	targetB = new TH1F("targetB","",nNewBins,0,1);

 	for (int bin = 1;bin<=nNewBins;bin++){

 	  double retVals = 0;
  	  double retValb = 0;

  	  double lval = (1./nNewBins)*(bin-1);
  	  double hval = (1./nNewBins)*bin;

  	  for (int k=1;k<=n2dbinsX;k++){
    	   for (int l=1;l<=n2dbinsY;l++){

      		  double ns = targetS2d->GetBinContent(k,l);
      		  double nb = targetB2d->GetBinContent(k,l);
		  double binsb;
		  if (nb < 0.0001){
      		  	binsb = 1 - delta;
		  } else {

      		  	binsb = SBscale*defx(ns/nb) - delta;
		  }

      		  if ( binsb < hval && binsb > lval ){
        		retVals+=ns;
        		retValb+=nb;
      		  }
    	   }
  	  }
	
  	  targetS->SetBinContent(bin,retVals);
          targetB->SetBinContent(bin,retValb);

 	}

	// now find the optimial ranges 
	std::vector<double> myBins = significanceOptimizedBinning();
	nFinalBins = myBins.size()-1;

	int rangecounter = 1;
	for (std::vector<double>::iterator it = myBins.begin();it!=myBins.end();it++){

        // Loop Through the bins and set the values of the new histogram (remember, h2 already scaled) 
        double lval = *it;

        for (int k=1;k<=n2dbinsX;k++){
          for (int l=1;l<=n2dbinsY;l++){
		
      		  double ns = targetS2d->GetBinContent(k,l);
      		  double nb = targetB2d->GetBinContent(k,l);
		  double binsb;
		  if (nb < 0.0001){
      		  	binsb = 1 - delta; // All the S/B 
		  } else {

      		  	binsb = SBscale*defx(ns/nb) - delta;
		  }

                if ( binsb > lval ){ // Note will never be > 1
                        FinalHist->SetBinContent(k,l,((double)rangecounter/nFinalBins)-(0.5/nFinalBins));  
                }
          }
        }
        rangecounter++;
 }

}

// These are the work horses of the optimization.

// find the maximum entry in a 2D hist
double Optimizations::findMaximum(TH2F *hist){

  double max = 0.;

  for (int k=1;k<=n2dbinsX;k++){
    for (int l=1;l<=n2dbinsY;l++){

      if ( hist->GetBinContent(k,l)>max ){
        max = hist->GetBinContent(k,l);
      }
    }
  }
 return max;
}

void Optimizations::getIntegralBetweenRanges(double *s,double *b,int l,     int h){
 
   *s = targetS->Integral(l,h);
   *b = targetB->Integral(l,h);
 
}

bool Optimizations::checkMinBkg(double bkg){

        return ( bkg > predefMin );

}

double Optimizations::calculateSigMulti(double *s1, double *b1, int nchannel){

  //    int nchannel=s1.size();
  double sterm=0;
  double logterms=0;
  for (int i=0;i<nchannel;i++){
    if (!checkMinBkg(b1[i])) return 0.;
    logterms+=(s1[i]+b1[i])*TMath::Log((s1[i]+b1[i])/b1[i]);
    sterm+=s1[i];
  }
  double sig =  1.4142*TMath::Sqrt(logterms - sterm);
  return sig;
}

// RECURSIVE FUNCTION ------------------------------------------------- //
void Optimizations::maxSigScan(double *maximumSignificance,int *frozen_counters,int *chosen_counters, int N,int *counters, int movingCounterIndex){


  if (movingCounterIndex >=N ) std::cout << "Should never have got here!!!" <<std::endl;
  if (movingCounterIndex < 0) return;
  if (counters[movingCounterIndex] < 2) std::cout << "WHAT IS GOING ON?? " <<  movingCounterIndex << " " << counters[movingCounterIndex]<<std::endl;
  int nBins = nNewBins;

  double significance_now;
  if (not sweepmode){
    // N is number of boundaries
    int m=1;
    for (int k=movingCounterIndex+1;k<N;k++) {
      counters[k]=counters[movingCounterIndex]+m;
      m++;
    }
  } else {
    // N is number of boundaries
    int m=1;
    for (int k=movingCounterIndex+1;k<N;k++) {
      int newpoint = (frozen_counters[k] - g_step > 1) ? frozen_counters[k] - g_step:2;
      counters[k]= (newpoint > counters[movingCounterIndex]+m )? newpoint: counters[movingCounterIndex]+m ;
      m++;
    }
  }

  if (counters[movingCounterIndex] < 2) std::cout << "WHAT IS GOING ON?? " <<  movingCounterIndex << " " << counters[movingCounterIndex]<<std::endl;

  if ( movingCounterIndex==N-1) {	
    if (not sweepmode){
      for (;counters[N-1]<=nBins;counters[N-1]+=g_step){
        for (int j=0;j<=N-1;j++){
          if (j==0){

	    
          getIntegralBetweenRanges(&signalVector1[j],&backgroundVector1[j],1,counters[j]-1);

          } else {
          getIntegralBetweenRanges(&signalVector1[j],&backgroundVector1[j],counters[j-1],counters[j]-1);

          }
        }

        getIntegralBetweenRanges(&signalVector1[N],&backgroundVector1[N],counters[N-1],nBins);
        //backgroundVector1[N]=(hb->Integral(counters[N-1],nBins));
        significance_now = calculateSigMulti(signalVector1,backgroundVector1,N+1);

        if (significance_now>*maximumSignificance){
          *maximumSignificance=significance_now;
          for (int j=0;j<N;j++){
            chosen_counters[j]=counters[j];	
            if (chosen_counters[j] < 0) std::cout << "Freak OUT !!! - " << j << "  " << chosen_counters[j] <<std::endl;
          }
        }
      }
      maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex-1);
    } else { // fine scanning
      int currmax = (nBins < frozen_counters[N-1] + g_step )? nBins : frozen_counters[N-1] + g_step;
      for (;counters[N-1]<=currmax;counters[N-1]++){
        for (int j=0;j<=N-1;j++){
          if (j==0){
            getIntegralBetweenRanges(&signalVector1[j],&backgroundVector1[j],1,counters[j]-1);

          } else {
            getIntegralBetweenRanges(&signalVector1[j],&backgroundVector1[j],counters[j-1],counters[j]-1);

          }
        }	
        getIntegralBetweenRanges(&signalVector1[N],&backgroundVector1[N],counters[N-1],nBins);
        significance_now = calculateSigMulti(signalVector1,backgroundVector1,N+1);
        if (significance_now>*maximumSignificance){
          *maximumSignificance=significance_now;
          for (int j=0;j<N;j++){
            chosen_counters[j]=counters[j];
          }
        }
      }
      maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex-1);
    }
  }


  else if (counters[movingCounterIndex]+1 <= nBins-(N-movingCounterIndex)){

    if (not sweepmode){
      if (counters[movingCounterIndex]+g_step <= nBins-(N-movingCounterIndex)){
        counters[movingCounterIndex]+=g_step;
      } else {
        counters[movingCounterIndex]=nBins-(N-movingCounterIndex);
      }
      int m=1;
      for (int k=movingCounterIndex+1;k<N;k++) {
        counters[k]=counters[movingCounterIndex]+m;
        m++;
      }
      maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex+1);
    } else {
      // N is number of boundaries
      if (counters[movingCounterIndex]+1 <= frozen_counters[movingCounterIndex]+g_step){
        counters[movingCounterIndex]++;
        int m=1;
        for (int k=movingCounterIndex+1;k<N;k++) {
          int newpoint = (frozen_counters[k] - g_step > 1) ? frozen_counters[k] - g_step:2;
          if (newpoint < 0) {std::cout << "Whaaaaaa? "<< newpoint <<std::endl;}
          counters[k]= (newpoint > counters[movingCounterIndex]+m )? newpoint: counters[movingCounterIndex]+m ;
          m++;
        }
        maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex+1);
      } else {
        if (movingCounterIndex>0){
          maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex-1);
        } else {
          return;
        }
      }
    }

  }

  else { // got to the end,

    if (movingCounterIndex>0){
      maxSigScan(maximumSignificance,frozen_counters,chosen_counters,N,counters,movingCounterIndex-1);
    } else {
      return;
    }
  }
}
// RECURSIVE FUNCTION ------------------------------------------------- //

std::vector<double> Optimizations::significanceOptimizedBinning(){


//  TH2F *hsnew = (TH2F*) hs->Clone();
//  TH2F *hbnew = (TH2F*) hb->Clone();

  // Here is the Maximisation Algorithm
  int 	*counters, *chosen_counters,*frozen_counters;
  double 	highestMaxSignificance=0;
  int 	chosenN=1;
  int 	*finalCounters=NULL ;
	
  // Optimizes the number of big steps to take 
  g_step = (int)TMath::Exp(TMath::Log(nNewBins/2)/2);
  if (g_step < 1) g_step=1;

  int Retry=0;

  for (int N=1;N<maxNumberOfBins;N++){  // Refuse to go beyond 8 Bins, will take forever
    sweepmode=0;	// First perform Broad Scan with optimized step size (g_step)
    bool skipBroad = false;
    if ( nNewBins < (N-1+2+Retry) ) {std::cout << "Forced to perform Fine scan since all the Retries failed to find a nice minimum :("<<std::endl; skipBroad=true;}

    double maximumSignificance=0;
    counters = new int[N];
    chosen_counters = new int[N];
    frozen_counters = new int[N];
    for (int c=0;c<N;c++) counters[c]=c+2+Retry;	// init to starting values
    for (int c=0;c<N;c++) frozen_counters[c]=c+2+Retry; // init to starting values
    for (int c=0;c<N;c++) chosen_counters[c]=c+2+Retry; // init to starting values

    double diff;
    clock_t start;

    std::cout << "Performing Fully optimized Scan"	<<std::endl;
    start=clock();
    if (!skipBroad) maxSigScan(&maximumSignificance,frozen_counters,chosen_counters,N,counters,N-1);

    sweepmode=1;	// Now do Fine scan after having found rough maximum
    for (int c=0;c<N;c++) counters[c]=chosen_counters[c]; // init to rough guess
    for (int c=0;c<N;c++) frozen_counters[c]=chosen_counters[c];  // init to rough guess

    // For full scanning, need to move lowest boundary to lowest point now
    int resetpoint = (2>frozen_counters[0]-g_step) ? 2 : frozen_counters[0]-g_step;
    counters[0]=resetpoint;

    maxSigScan(&maximumSignificance,frozen_counters,chosen_counters,N,counters,N-1);

    diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
    std::cout << Form("Finished, time taken = %3.5f",diff)<<std::endl;
    std::cout << "N Bins, Max Significance -> " << N+1 << " "<<maximumSignificance << std::endl;
    /*std::cout << "Boundaries at: [ ";
    for (int cc=0;cc<N; cc++) std::cout << hsnew->GetBinLowEdge(chosen_counters[cc]) << " , ";
    std::cout << " ]" << std::endl;
    */

    if (maximumSignificance < highestMaxSignificance){
         
         std::cout << "Looks like the Broad Scan found a local maxmimum and got stuck there (can happen if the initial number of bins is small), Try again " <<std::endl;
         cout << "Chosen bins:   "; 
         for (int c=0;c<N;c++) cout << chosen_counters[c] << "    ";
         cout << endl;
         N--; Retry++;
         continue;
 
    }
    else Retry =0;


    if ((maximumSignificance-highestMaxSignificance)/highestMaxSignificance > threshold){
      highestMaxSignificance = maximumSignificance ;
      finalCounters= new int[N];
      chosenN = N;
      for (int cc=0;cc<N;cc++) finalCounters[cc]=chosen_counters[cc];
    } else {
      
      break;
    }

  }

  std::vector<double> newbinEdges;
  
  newbinEdges.push_back(-1); // Pretend he's a BDT
  std::cout << "-1,";
  for (int newguy=0;newguy<chosenN;newguy++){
    double val = (1./nNewBins)*finalCounters[newguy];
    std::cout << val << ", ";
    newbinEdges.push_back(val);
  }
  std::cout << " 1" << std::endl;

  newbinEdges.push_back(1); // Upper Edge

  return newbinEdges;
}


