morph_code = """#include <iostream>
#include <math.h>
#include "TH1F.h"
#include "TH1D.h"
#include <TROOT.h>
using namespace std;

TH1D *th1dmorph(Char_t *chname="TH1D-interpolated", 
		Char_t *chtitle="Interpolated histogram",
		TH1D *hist1=0,TH1D *hist2=0,
		Double_t par1=0,Double_t par2=1,Double_t parinterp=0,
		Double_t morphedhistnorm=-1,
		Int_t idebug=0)
{
  //--------------------------------------------------------------------------
  // Author           : Alex Read 
  // Version 0.3 of ROOT implementation, 08.05.2011
  //
  // Changes 0.3->0.31 17.07.2011:
  //    o Squashed bug that gave errors for histograms with holes.
  //
  // Downloaded from https://wiki.uio.no/mn/fys/epf/index.php/Linear_interpolation_of_histograms


  // Return right away if one of the input histograms doesn't exist.
  if(!hist1) {
    cout << "ERROR! th1morph says first input histogram doesn't exist." << endl;
    return(0);
  }
  if(!hist2) {
    cout << "ERROR! th1morph says second input histogram doesn't exist." << endl;
    return(0);
  }
  
  // Extract bin parameters of input histograms 1 and 2. We haven't implemented
  // nonuniform binning (yet?).

  Int_t nb1 = hist1->GetNbinsX();
  Int_t nb2 = hist2->GetNbinsX();
  Double_t xmin1 = hist1->GetXaxis()->GetXmin();
  Double_t xmin2 = hist2->GetXaxis()->GetXmin();
  Double_t xmax1 = hist1->GetXaxis()->GetXmax();
  Double_t xmax2 = hist2->GetXaxis()->GetXmax();
  if (idebug > 0) {
    cout << nb1 << " " << xmin1 << " " << xmax1 << endl;
    cout << nb2 << " " << xmin2 << " " << xmax2 << endl;
  }

// ......The weights (wt1,wt2) are the complements of the "distances" between 
//       the values of the parameters at the histograms and the desired 
//       interpolation point. For example, wt1=0, wt2=1 means that the 
//       interpolated histogram should be identical to input histogram 2.
//       Check that they make sense. If par1=par2 then we can choose any
//       valid set of wt1,wt2 so why not take the average?

  Double_t wt1,wt2;
  if (par2 != par1) {
    wt1 = 1. - (parinterp-par1)/(par2-par1);
    wt2 = 1. + (parinterp-par2)/(par2-par1);
  }
  else { 
    wt1 = 0.5;
    wt2 = 0.5;
  }

  //......Give a warning if this is an extrapolation.

  if (wt1 < 0 || wt1 > 1. || wt2 < 0. || wt2 > 1. || fabs(1-(wt1+wt2)) 
      > 1.0e-4) {
    cout << "Warning! th1morph: This is an extrapolation!! Weights are "
	 << wt1 << " and " << wt2 << " (sum=" << wt1+wt2 << ")" << endl;
  }
  if (idebug >= 1) cout << "th1morph - Weights: " << wt1 << " " << wt2 << endl;

  //
  //......Perform the interpolation of histogram bin parameters. Use
  //      assignments instead of computation when input binnings
  //      are identical to assure best possible precision.

  Double_t xminn=-1,xmaxn=-1;
  Int_t nbn=0;
  Double_t wtmin;

  wtmin = wt1; if (wt2 < wt1) wtmin = wt2;

  if (wtmin >= 0) {
    if (xmin1 == xmin2) {
      xminn = xmin1;
    } else {
      xminn = wt1*xmin1 + wt2*xmin2;
    }
    if (xmax1 == xmax2) {
      xmaxn = xmax1;
    } else {
      xmaxn = wt1*xmax1 + wt2*xmax2;
    }
    if (nb1 == nb2) {
      nbn = nb1;
    } else {
      nbn   = wt1*nb1   + wt2*nb2;
    }
  }
  //......If one of the weights is zero, then use the binnings of the
  //      histogram with nonzero weight.
  //      but reasonable with the histogram bin parameters.
  else {
    if (wt1 == 0) {
      xminn = xmin2; xmaxn = xmax2; nbn = nb2;
    } else if (wt2 == 0) {
      xminn = xmin1; xmaxn = xmax1; nbn = nb1;
    }
  }
  if (idebug >= 1) cout << "New hist: " << nbn << " " << xminn << " " 
			<< xmaxn << endl;

  // Treatment for empty histograms: Return an empty histogram
  // with interpolated bins.

  if (hist1->GetSum() <= 0 || hist2->GetSum() <=0 ) {
    cout << "Warning! th1morph detects an empty input histogram. Empty interpolated histogram returned: " 
	 <<endl << "         " << chname << " - " << chtitle << endl;
    TH1D *morphedhist = (TH1D *)gROOT->FindObject(chname);
    if (morphedhist) delete morphedhist;
    morphedhist = new TH1D(chname,chtitle,nbn,xminn,xmaxn);
    return(morphedhist);
  }
  if (idebug >= 1) cout << "Input histogram content sums: " 
			<< hist1->GetSum() << " " << hist2->GetSum() << endl;
// *         
// *......Extract the single precision histograms into double precision arrays
// *      for the interpolation computation. The offset is because sigdis(i)
// *      describes edge i (there are nbins+1 of them) while dist1/2
// *      describe bin i. Be careful, ROOT does not use C++ convention to
// *      number bins: dist1[ibin] is content of bin ibin where ibin runs from
// *      1 to nbins. We allocate some extra space for the derived distributions
// *      because there may be as many as nb1+nb2+2 edges in the intermediate 
// *      interpolated cdf described by xdisn[i] (position of edge i) and 
// *      sigdisn[i] (cummulative probability up this edge) before we project 
// *      into the final binning.

  Double_t *dist1=hist1->GetArray(); 
  Double_t *dist2=hist2->GetArray();
  Double_t *sigdis1 = new Double_t[1+nb1];
  Double_t *sigdis2 = new Double_t[1+nb2];
  Double_t *sigdisn = new Double_t[2+nb1+nb2];
  Double_t *xdisn = new Double_t[2+nb1+nb2];
  Double_t *sigdisf = new Double_t[nbn+1];

  for(Int_t i=0;i<2+nb1+nb2;i++) xdisn[i] = 0; // Start with empty edges
  sigdis1[0] = 0; sigdis2[0] = 0; // Start with cdf=0 at left edge

  for(Int_t i=1;i<nb1+1;i++) {   // Remember, bin i has edges at i-1 and 
    sigdis1[i] = dist1[i];       // i and i runs from 1 to nb.
  }
  for(Int_t i=1;i<nb2+1;i++) {
    sigdis2[i] = dist2[i];
  }

  if (idebug >= 3) {
    for(Int_t i=0;i<nb1+1;i++) {
      cout << i << " dist1" << dist1[i] << endl;
    }
    for(Int_t i=0;i<nb1+1;i++) {
      cout << i << " dist2" << dist1[i] << endl;
    }
  }
  
//......Normalize the distributions to 1 to obtain pdf's and integrate 
//      (sum) to obtain cdf's.

  Double_t total = 0, norm1, norm2;
  for(Int_t i=0;i<nb1+1;i++) {
    total += sigdis1[i];
  }
  if (idebug >=1) cout << "Total histogram 1: " <<  total << endl;
  for(Int_t i=1;i<nb1+1;i++) {
    sigdis1[i] = sigdis1[i]/total + sigdis1[i-1];
  }
  norm1 = total;
  
  total = 0.;
  for(Int_t i=0;i<nb2+1;i++) {
    total += sigdis2[i];
  }
  if (idebug >=1) cout << "Total histogram 22: " <<  total << endl;
  for(Int_t i=1;i<nb2+1;i++) {
    sigdis2[i] = sigdis2[i]/total + sigdis2[i-1];
  }
  norm2 = total;

// *
// *......We are going to step through all the edges of both input
// *      cdf's ordered by increasing value of y. We start at the
// *      lower edge, but first we should identify the upper ends of the
// *      curves. These (ixl1, ixl2) are the first point in each cdf from 
// *      above that has the same integral as the last edge.
// *

  Int_t ix1l = nb1;
  Int_t ix2l = nb2;
  while(sigdis1[ix1l-1] >= sigdis1[ix1l]) {
    ix1l = ix1l - 1;
  }
  while(sigdis2[ix2l-1] >= sigdis2[ix2l]) {
    ix2l = ix2l - 1;
  }

// *
// *......Step up to the beginnings of the curves. These (ix1, ix2) are the
// *      first non-zero points from below.

  Int_t ix1 = -1;
  do {
    ix1 = ix1 + 1;
  } while(sigdis1[ix1+1] <= sigdis1[0]);

  Int_t ix2 = -1;
  do {
    ix2 = ix2 + 1;
  } while(sigdis2[ix2+1] <= sigdis2[0]);

  if (idebug >= 1) {
    cout << "First and last edge of hist1: " << ix1 << " " << ix1l << endl;
    cout << "   " << sigdis1[ix1] << " " << sigdis1[ix1+1] << endl;
    cout << "First and last edge of hist2: " << ix2 << " " << ix2l << endl;
    cout << "   " << sigdis2[ix2] << " " << sigdis2[ix2+1] << endl;
  }

  //.......Need bin widths

  Double_t dx1=(xmax1-xmin1)/double(nb1);
  Double_t dx2=(xmax2-xmin2)/double(nb2);
  Double_t dx=(xmaxn-xminn)/double(nbn);

  //......The first interpolated point should be computed now.

  Int_t nx3 = 0;
  Double_t x1,x2,x;
  x1 = xmin1 + double(ix1)*dx1;
  x2 = xmin2 + double(ix2)*dx2;
  x = wt1*x1 + wt2*x2;
  xdisn[nx3] = x;
  sigdisn[nx3] = 0;
  if(idebug >= 1) {
    cout << "First interpolated point: " << xdisn[nx3] << " " 
	 << sigdisn[nx3] << endl;
    cout << "                          " << x1 << " <= " << x << " <= " 
	 << x2 << endl;
  }

  //......Loop over the remaining point in both curves. Getting the last
  //      points may be a bit tricky due to limited floating point 
  //      precision.

  if (idebug >= 1) {
      cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l 
	   << ", ix2=" << ix2 << ", ix2l=" << ix2l << endl;
      for(Int_t i=ix1;i<=ix1l;i++) {cout << "   1: " << i << " " << sigdis1[i] << endl;}
      for(Int_t i=ix2;i<=ix2l;i++) {cout << "   2: " << i << " " << sigdis2[i] << endl;}
  }

  Double_t yprev = -1; // The probability y of the previous point, it will 
                       //get updated and used in the loop.
  Double_t y,x20,x21,y20,y21; // Interpolation points along cdfs 0,1,2
  Double_t x10,x11,y10,y11;

  while(ix1 < ix1l | ix2 < ix2l) {
    if (idebug >= 1 ) cout << "----Top of while with ix1=" << ix1 
			   << ", ix1l=" << ix1l << ", ix2=" << ix2 
			   << ", ix2l=" << ix2l << endl;

    //......Increment to the next lowest point. Step up to the next
    //      kink in case there are several empty (flat in the integral)
    //      bins.

    Int_t i12type = -1; // Tells which input distribution we need to 
                        // see next point of.

    if ((sigdis1[ix1+1] <= sigdis2[ix2+1] || ix2 == ix2l) && ix1 < ix1l) {
      ix1 = ix1 + 1;
      // try to fix empty bin holes!!!!
      //      while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
      // 	ix1 = ix1 + 1;
      //      }
      //empty bin fix??? while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
      while(sigdis1[ix1+1] < sigdis1[ix1] && ix1 < ix1l) {
 	ix1 = ix1 + 1;
      }
      i12type = 1;
    } else if (ix2 < ix2l) {
      ix2 = ix2 + 1;
      //empty bin fix ?? while(sigdis2[ix2+1] <= sigdis2[ix2] && ix2 < ix2l) {
      while(sigdis2[ix2+1] < sigdis2[ix2] && ix2 < ix2l) {
 	ix2 = ix2 + 1;
      }
      i12type = 2;
    }
    if (i12type == 1) {
      if (idebug >= 3) {
	cout << "Pair for i12type=1: " << ix1 << " " << ix2 << " " << sigdis2[ix2] << " " 
	     << sigdis1[ix1] << " " << sigdis2[ix2+1] << endl;
      }
      x1 = xmin1 + double(ix1)*dx1 ;
      y = sigdis1[ix1];
      x20 = double(ix2)*dx2 + xmin2;
      x21 = x20 + dx2;
      y20 = sigdis2[ix2];
      y21 = sigdis2[ix2+1];

      //......Calculate where the cummulative probability y in distribution 1
      //      intersects between the 2 points from distribution 2 which 
      //      bracket it.

      if (y21 > y20) {
	x2 = x20 + (x21-x20)*(y-y20)/(y21-y20);
      } 
      else {
	x2 = x20;
      }
    } else {
      if (idebug >= 3) {
	cout << "Pair for i12type=2: " << sigdis1[ix1] << " " << sigdis2[ix2] 
	     << " " << sigdis1[ix1+1] << endl;
      }
      x2 = xmin2 + double(ix2)*dx2 ;
      y = sigdis2[ix2];
      x10 = double(ix1)*dx1 + xmin1;
      x11 = x10 + dx1;
      y10 = sigdis1[ix1];
      y11 = sigdis1[ix1+1];

      //......Calculate where the cummulative probability y in distribution 2
      //      intersects between the 2 points from distribution 1 which 
      //      brackets it.

      if (y11 > y10) {
	x1 = x10 + (x11-x10)*(y-y10)/(y11-y10);
      } else {
	x1 = x10;
      }
    }

    //......Interpolate between the x's in the 2 distributions at the 
    //      cummulative probability y. Store the (x,y) for provisional 
    //      edge nx3 in (xdisn[nx3],sigdisn[nx3]). nx3 grows for each point
    //      we add the the arrays. Note: Should probably turn the pair into 
    //      a structure to make the code more object-oriented and readable.

    x = wt1*x1 + wt2*x2;
    if (y >= yprev) { // bugfix for empty bins?!?!?!?!
      nx3 = nx3+1;
      if (idebug >= 1) {
	cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" 
	     << nx3 << ", x= " << x << ", y=" << y << ", yprev=" << yprev 
	     << endl;
      }
      yprev = y;
      xdisn[nx3] = x;
      sigdisn[nx3] = y;
      if(idebug >= 1) {
	cout << "    ix1=" << ix1 << ", ix2= " << ix2 << ", i12type= " 
	     << i12type << ", sigdis1[ix1]=" << sigdis1[ix1] << endl;
	cout << "        " << ", nx3=" << nx3 << ", x=" << x << ", y= " 
	     << sigdisn[nx3] << endl;
      }
    }
  }
  if (idebug >=3) for (Int_t i=0;i<=nx3;i++) {
    cout << " nx " << i << " " << xdisn[i] << " " << sigdisn[i] << endl;
  }

  // *......Now we loop over the edges of the bins of the interpolated
  // *      histogram and find out where the interpolated cdf 3
  // *      crosses them. This projection defines the result and will
  // *      be stored (after differention and renormalization) in the
  // *      output histogram.
  // *
  // *......We set all the bins following the final edge to the value
  // *      of the final edge.

  x = xminn + double(nbn)*dx;
  Int_t ix = nbn;

  if (idebug >= 1) cout << "------> Any final bins to set? " << x << " " 
			<< xdisn[nx3] << endl;
  while(x >= xdisn[nx3]) {
    sigdisf[ix] = sigdisn[nx3];
    if (idebug >= 2) cout << "   Setting final bins" << ix << " " << x 
			  << " " << sigdisf[ix] << endl;
    ix = ix-1;
    x = xminn + double(ix)*dx;
  }
  Int_t ixl = ix + 1;
  if (idebug >= 1) cout << " Now ixl=" << ixl << " ix=" << ix << endl;

  // *
  // *......The beginning may be empty, so we have to step up to the first
  // *      edge where the result is nonzero. We zero the bins which have
  // *      and upper (!) edge which is below the first point of the
  // *      cummulative distribution we are going to project to this
  // *      output histogram binning.
  // *

  ix = 0;
  x = xminn + double(ix+1)*dx;
  if (idebug >= 1) cout << "Start setting initial bins at x=" << x << endl;
  while(x <= xdisn[0]) {
    sigdisf[ix] = sigdisn[0];
    if (idebug >= 1) cout << "   Setting initial bins " << ix << " " << x 
			  << " " << xdisn[1] << " " << sigdisf[ix] << endl;
    ix = ix+1;
    x = xminn + double(ix+1)*dx;
  }
  Int_t ixf = ix;

  if (idebug >= 1)
    cout << "Bins left to loop over:" << ixf << "-" << ixl << endl;

  // *......Also the end (from y to 1.0) often comes before the last edge
  // *      so we have to set the following to 1.0 as well.

  Int_t ix3 = 0; // Problems with initial edge!!!
  for(ix=ixf;ix<ixl;ix++) {
    x = xminn + double(ix)*dx;
    if (x < xdisn[0]) {
      y = 0;
    } else if (x > xdisn[nx3]) {
      y = 1.;
    } else {
      while(xdisn[ix3+1] <= x && ix3 < 2*nbn) {
	ix3 = ix3 + 1;
      }
      if (xdisn[ix3+1]-x > 1.1*dx2) { // Empty bin treatment
	y = sigdisn[ix3+1]; //y = sigdisn[ix3+1];
	if(idebug>=1) cout << "Empty bin treatment " << ix3+1 << " " << y << " " << sigdisn[ix3+2] << " " << sigdisn[ix3-1] << endl;  
      }
      else if (xdisn[ix3+1] > xdisn[ix3]) { // Normal bins
	y = sigdisn[ix3] + (sigdisn[ix3+1]-sigdisn[ix3])
	  *(x-xdisn[ix3])/(xdisn[ix3+1]-xdisn[ix3]);
      } else {  // Is this ever used?
	y = 0;
	cout << "Warning - th1morph: This probably shoudn't happen! " 
	     << endl;
	cout << "Warning - th1morph: Zero slope solving x(y)" << endl;
      }
    }
    sigdisf[ix] = y;
    if (idebug >= 3) {
      cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3] << ", x=" 
	   << x << ", next xdisn=" << xdisn[ix3+1] << endl;
      cout << "   cdf n=" << sigdisn[ix3] << ", y=" << y << ", next point=" 
	   << sigdisn[ix3+1] << endl;
    }
  }

  //......Differentiate interpolated cdf and return renormalized result in 
  //      new histogram. 

  TH1D *morphedhist = (TH1D *)gROOT->FindObject(chname);
  if (morphedhist) delete morphedhist;
  morphedhist = new TH1D(chname,chtitle,nbn,xminn,xmaxn);
 
  Double_t norm = morphedhistnorm;
  // norm1, norm2, wt1, wt2 are computed before the interpolation
  if (norm <= 0) {
    if (norm1 == norm2) {
      norm = norm1;
    } else {
      norm   = wt1*norm1   + wt2*norm2;
    }
  }

  for(Int_t ixx=nbn-1;ixx>-1;ixx--) {
    x = xminn + double(ixx)*dx;
    y =  sigdisf[ixx+1]-sigdisf[ixx];
    if (y<0) cout << "huh??? " << ixx << " " << sigdisf[ixx] << " " << sigdisf[ixx+1] << endl;
    morphedhist->SetBinContent(ixx+1,y*norm);
  }
  
  //......Clean up the temporary arrays we allocated.

  delete sigdis1; delete sigdis2; 
  delete sigdisn; delete xdisn; delete sigdisf;

  //......All done, return the result.
  //morphedhist->Draw("same");
  return(morphedhist);
}

TH1F *th1fmorph(Char_t *chname="TH1F-interpolated", 
		Char_t *chtitle="Interpolated histogram",
		TH1F *hist1=0,TH1F *hist2=0,
		Double_t par1=0,Double_t par2=1,Double_t parinterp=0,
		Double_t morphedhistnorm=-1,
		Int_t idebug=0)
{
  //--------------------------------------------------------------------------
  // Author           : Alex Read 
  // Version 0.3 of ROOT implementation, 08.05.2011
  //
  // Changes 0.3->0.31 17.07.2011:
  //    o Squashed bug that gave errors for histograms with holes.
  //
  // Downloaded from https://wiki.uio.no/mn/fys/epf/index.php/Linear_interpolation_of_histograms


  // Return right away if one of the input histograms doesn't exist.
  if(!hist1) {
    cout << "ERROR! th1morph says first input histogram doesn't exist." << endl;
    return(0);
  }
  if(!hist2) {
    cout << "ERROR! th1morph says second input histogram doesn't exist." << endl;
    return(0);
  }
  
  // Extract bin parameters of input histograms 1 and 2. We haven't implemented
  // nonuniform binning (yet?).

  Int_t nb1 = hist1->GetNbinsX();
  Int_t nb2 = hist2->GetNbinsX();
  Double_t xmin1 = hist1->GetXaxis()->GetXmin();
  Double_t xmin2 = hist2->GetXaxis()->GetXmin();
  Double_t xmax1 = hist1->GetXaxis()->GetXmax();
  Double_t xmax2 = hist2->GetXaxis()->GetXmax();
  if (idebug > 0) {
    cout << nb1 << " " << xmin1 << " " << xmax1 << endl;
    cout << nb2 << " " << xmin2 << " " << xmax2 << endl;
  }

// ......The weights (wt1,wt2) are the complements of the "distances" between 
//       the values of the parameters at the histograms and the desired 
//       interpolation point. For example, wt1=0, wt2=1 means that the 
//       interpolated histogram should be identical to input histogram 2.
//       Check that they make sense. If par1=par2 then we can choose any
//       valid set of wt1,wt2 so why not take the average?

  Double_t wt1,wt2;
  if (par2 != par1) {
    wt1 = 1. - (parinterp-par1)/(par2-par1);
    wt2 = 1. + (parinterp-par2)/(par2-par1);
  }
  else { 
    wt1 = 0.5;
    wt2 = 0.5;
  }

  //......Give a warning if this is an extrapolation.

  if (wt1 < 0 || wt1 > 1. || wt2 < 0. || wt2 > 1. || fabs(1-(wt1+wt2)) 
      > 1.0e-4) {
    cout << "Warning! th1morph: This is an extrapolation!! Weights are "
	 << wt1 << " and " << wt2 << " (sum=" << wt1+wt2 << ")" << endl;
  }
  if (idebug >= 1) cout << "th1morph - Weights: " << wt1 << " " << wt2 << endl;

  //
  //......Perform the interpolation of histogram bin parameters. Use
  //      assignments instead of computation when input binnings
  //      are identical to assure best possible precision.

  Double_t xminn=-1,xmaxn=-1;
  Int_t nbn=0;
  Double_t wtmin;

  wtmin = wt1; if (wt2 < wt1) wtmin = wt2;

  if (wtmin >= 0) {
    if (xmin1 == xmin2) {
      xminn = xmin1;
    } else {
      xminn = wt1*xmin1 + wt2*xmin2;
    }
    if (xmax1 == xmax2) {
      xmaxn = xmax1;
    } else {
      xmaxn = wt1*xmax1 + wt2*xmax2;
    }
    if (nb1 == nb2) {
      nbn = nb1;
    } else {
      nbn   = wt1*nb1   + wt2*nb2;
    }
  }
  //......If one of the weights is zero, then use the binnings of the
  //      histogram with nonzero weight.
  //      but reasonable with the histogram bin parameters.
  else {
    if (wt1 == 0) {
      xminn = xmin2; xmaxn = xmax2; nbn = nb2;
    } else if (wt2 == 0) {
      xminn = xmin1; xmaxn = xmax1; nbn = nb1;
    }
  }
  if (idebug >= 1) cout << "New hist: " << nbn << " " << xminn << " " 
			<< xmaxn << endl;

  // Treatment for empty histograms: Return an empty histogram
  // with interpolated bins.

  if (hist1->GetSum() <= 0 || hist2->GetSum() <=0 ) {
    cout << "Warning! th1morph detects an empty input histogram. Empty interpolated histogram returned: " 
	 <<endl << "         " << chname << " - " << chtitle << endl;
    TH1F *morphedhist = (TH1F *)gROOT->FindObject(chname);
    if (morphedhist) delete morphedhist;
    morphedhist = new TH1F(chname,chtitle,nbn,xminn,xmaxn);
    return(morphedhist);
  }
  if (idebug >= 1) cout << "Input histogram content sums: " 
			<< hist1->GetSum() << " " << hist2->GetSum() << endl;
// *         
// *......Extract the single precision histograms into double precision arrays
// *      for the interpolation computation. The offset is because sigdis(i)
// *      describes edge i (there are nbins+1 of them) while dist1/2
// *      describe bin i. Be careful, ROOT does not use C++ convention to
// *      number bins: dist1[ibin] is content of bin ibin where ibin runs from
// *      1 to nbins. We allocate some extra space for the derived distributions
// *      because there may be as many as nb1+nb2+2 edges in the intermediate 
// *      interpolated cdf described by xdisn[i] (position of edge i) and 
// *      sigdisn[i] (cummulative probability up this edge) before we project 
// *      into the final binning.

  Float_t *dist1=hist1->GetArray(); 
  Float_t *dist2=hist2->GetArray();
  Double_t *sigdis1 = new Double_t[1+nb1];
  Double_t *sigdis2 = new Double_t[1+nb2];
  Double_t *sigdisn = new Double_t[2+nb1+nb2];
  Double_t *xdisn = new Double_t[2+nb1+nb2];
  Double_t *sigdisf = new Double_t[nbn+1];

  for(Int_t i=0;i<2+nb1+nb2;i++) xdisn[i] = 0; // Start with empty edges
  sigdis1[0] = 0; sigdis2[0] = 0; // Start with cdf=0 at left edge

  for(Int_t i=1;i<nb1+1;i++) {   // Remember, bin i has edges at i-1 and 
    sigdis1[i] = dist1[i];       // i and i runs from 1 to nb.
  }
  for(Int_t i=1;i<nb2+1;i++) {
    sigdis2[i] = dist2[i];
  }

  if (idebug >= 3) {
    for(Int_t i=0;i<nb1+1;i++) {
      cout << i << " dist1" << dist1[i] << endl;
    }
    for(Int_t i=0;i<nb1+1;i++) {
      cout << i << " dist2" << dist1[i] << endl;
    }
  }
  
//......Normalize the distributions to 1 to obtain pdf's and integrate 
//      (sum) to obtain cdf's.

  Double_t total = 0, norm1, norm2;
  for(Int_t i=0;i<nb1+1;i++) {
    total += sigdis1[i];
  }
  if (idebug >=1) cout << "Total histogram 1: " <<  total << endl;
  for(Int_t i=1;i<nb1+1;i++) {
    sigdis1[i] = sigdis1[i]/total + sigdis1[i-1];
  }
  norm1 = total;
  
  total = 0.;
  for(Int_t i=0;i<nb2+1;i++) {
    total += sigdis2[i];
  }
  if (idebug >=1) cout << "Total histogram 22: " <<  total << endl;
  for(Int_t i=1;i<nb2+1;i++) {
    sigdis2[i] = sigdis2[i]/total + sigdis2[i-1];
  }
  norm2 = total;

// *
// *......We are going to step through all the edges of both input
// *      cdf's ordered by increasing value of y. We start at the
// *      lower edge, but first we should identify the upper ends of the
// *      curves. These (ixl1, ixl2) are the first point in each cdf from 
// *      above that has the same integral as the last edge.
// *

  Int_t ix1l = nb1;
  Int_t ix2l = nb2;
  while(sigdis1[ix1l-1] >= sigdis1[ix1l]) {
    ix1l = ix1l - 1;
  }
  while(sigdis2[ix2l-1] >= sigdis2[ix2l]) {
    ix2l = ix2l - 1;
  }

// *
// *......Step up to the beginnings of the curves. These (ix1, ix2) are the
// *      first non-zero points from below.

  Int_t ix1 = -1;
  do {
    ix1 = ix1 + 1;
  } while(sigdis1[ix1+1] <= sigdis1[0]);

  Int_t ix2 = -1;
  do {
    ix2 = ix2 + 1;
  } while(sigdis2[ix2+1] <= sigdis2[0]);

  if (idebug >= 1) {
    cout << "First and last edge of hist1: " << ix1 << " " << ix1l << endl;
    cout << "   " << sigdis1[ix1] << " " << sigdis1[ix1+1] << endl;
    cout << "First and last edge of hist2: " << ix2 << " " << ix2l << endl;
    cout << "   " << sigdis2[ix2] << " " << sigdis2[ix2+1] << endl;
  }

  //.......Need bin widths

  Double_t dx1=(xmax1-xmin1)/double(nb1);
  Double_t dx2=(xmax2-xmin2)/double(nb2);
  Double_t dx=(xmaxn-xminn)/double(nbn);

  //......The first interpolated point should be computed now.

  Int_t nx3 = 0;
  Double_t x1,x2,x;
  x1 = xmin1 + double(ix1)*dx1;
  x2 = xmin2 + double(ix2)*dx2;
  x = wt1*x1 + wt2*x2;
  xdisn[nx3] = x;
  sigdisn[nx3] = 0;
  if(idebug >= 1) {
    cout << "First interpolated point: " << xdisn[nx3] << " " 
	 << sigdisn[nx3] << endl;
    cout << "                          " << x1 << " <= " << x << " <= " 
	 << x2 << endl;
  }

  //......Loop over the remaining point in both curves. Getting the last
  //      points may be a bit tricky due to limited floating point 
  //      precision.

  if (idebug >= 1) {
      cout << "----BEFORE while with ix1=" << ix1 << ", ix1l=" << ix1l 
	   << ", ix2=" << ix2 << ", ix2l=" << ix2l << endl;
      for(Int_t i=ix1;i<=ix1l;i++) {cout << "   1: " << i << " " << sigdis1[i] << endl;}
      for(Int_t i=ix2;i<=ix2l;i++) {cout << "   2: " << i << " " << sigdis2[i] << endl;}
  }

  Double_t yprev = -1; // The probability y of the previous point, it will 
                       //get updated and used in the loop.
  Double_t y,x20,x21,y20,y21; // Interpolation points along cdfs 0,1,2
  Double_t x10,x11,y10,y11;

  while(ix1 < ix1l | ix2 < ix2l) {
    if (idebug >= 1 ) cout << "----Top of while with ix1=" << ix1 
			   << ", ix1l=" << ix1l << ", ix2=" << ix2 
			   << ", ix2l=" << ix2l << endl;

    //......Increment to the next lowest point. Step up to the next
    //      kink in case there are several empty (flat in the integral)
    //      bins.

    Int_t i12type = -1; // Tells which input distribution we need to 
                        // see next point of.

    if ((sigdis1[ix1+1] <= sigdis2[ix2+1] || ix2 == ix2l) && ix1 < ix1l) {
      ix1 = ix1 + 1;
      // try to fix empty bin holes!!!!
      //      while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
      // 	ix1 = ix1 + 1;
      //      }
      //empty bin fix??? while(sigdis1[ix1+1] <= sigdis1[ix1] && ix1 < ix1l) {
      while(sigdis1[ix1+1] < sigdis1[ix1] && ix1 < ix1l) {
 	ix1 = ix1 + 1;
      }
      i12type = 1;
    } else if (ix2 < ix2l) {
      ix2 = ix2 + 1;
      //empty bin fix ?? while(sigdis2[ix2+1] <= sigdis2[ix2] && ix2 < ix2l) {
      while(sigdis2[ix2+1] < sigdis2[ix2] && ix2 < ix2l) {
 	ix2 = ix2 + 1;
      }
      i12type = 2;
    }
    if (i12type == 1) {
      if (idebug >= 3) {
	cout << "Pair for i12type=1: " << ix1 << " " << ix2 << " " << sigdis2[ix2] << " " 
	     << sigdis1[ix1] << " " << sigdis2[ix2+1] << endl;
      }
      x1 = xmin1 + double(ix1)*dx1 ;
      y = sigdis1[ix1];
      x20 = double(ix2)*dx2 + xmin2;
      x21 = x20 + dx2;
      y20 = sigdis2[ix2];
      y21 = sigdis2[ix2+1];

      //......Calculate where the cummulative probability y in distribution 1
      //      intersects between the 2 points from distribution 2 which 
      //      bracket it.

      if (y21 > y20) {
	x2 = x20 + (x21-x20)*(y-y20)/(y21-y20);
      } 
      else {
	x2 = x20;
      }
    } else {
      if (idebug >= 3) {
	cout << "Pair for i12type=2: " << sigdis1[ix1] << " " << sigdis2[ix2] 
	     << " " << sigdis1[ix1+1] << endl;
      }
      x2 = xmin2 + double(ix2)*dx2 ;
      y = sigdis2[ix2];
      x10 = double(ix1)*dx1 + xmin1;
      x11 = x10 + dx1;
      y10 = sigdis1[ix1];
      y11 = sigdis1[ix1+1];

      //......Calculate where the cummulative probability y in distribution 2
      //      intersects between the 2 points from distribution 1 which 
      //      brackets it.

      if (y11 > y10) {
	x1 = x10 + (x11-x10)*(y-y10)/(y11-y10);
      } else {
	x1 = x10;
      }
    }

    //......Interpolate between the x's in the 2 distributions at the 
    //      cummulative probability y. Store the (x,y) for provisional 
    //      edge nx3 in (xdisn[nx3],sigdisn[nx3]). nx3 grows for each point
    //      we add the the arrays. Note: Should probably turn the pair into 
    //      a structure to make the code more object-oriented and readable.

    x = wt1*x1 + wt2*x2;
    if (y >= yprev) { // bugfix for empty bins?!?!?!?!
      nx3 = nx3+1;
      if (idebug >= 1) {
	cout << " ---> y > yprev: i12type=" << i12type << ", nx3=" 
	     << nx3 << ", x= " << x << ", y=" << y << ", yprev=" << yprev 
	     << endl;
      }
      yprev = y;
      xdisn[nx3] = x;
      sigdisn[nx3] = y;
      if(idebug >= 1) {
	cout << "    ix1=" << ix1 << ", ix2= " << ix2 << ", i12type= " 
	     << i12type << ", sigdis1[ix1]=" << sigdis1[ix1] << endl;
	cout << "        " << ", nx3=" << nx3 << ", x=" << x << ", y= " 
	     << sigdisn[nx3] << endl;
      }
    }
  }
  if (idebug >=3) for (Int_t i=0;i<=nx3;i++) {
    cout << " nx " << i << " " << xdisn[i] << " " << sigdisn[i] << endl;
  }

  // *......Now we loop over the edges of the bins of the interpolated
  // *      histogram and find out where the interpolated cdf 3
  // *      crosses them. This projection defines the result and will
  // *      be stored (after differention and renormalization) in the
  // *      output histogram.
  // *
  // *......We set all the bins following the final edge to the value
  // *      of the final edge.

  x = xminn + double(nbn)*dx;
  Int_t ix = nbn;

  if (idebug >= 1) cout << "------> Any final bins to set? " << x << " " 
			<< xdisn[nx3] << endl;
  while(x >= xdisn[nx3]) {
    sigdisf[ix] = sigdisn[nx3];
    if (idebug >= 2) cout << "   Setting final bins" << ix << " " << x 
			  << " " << sigdisf[ix] << endl;
    ix = ix-1;
    x = xminn + double(ix)*dx;
  }
  Int_t ixl = ix + 1;
  if (idebug >= 1) cout << " Now ixl=" << ixl << " ix=" << ix << endl;

  // *
  // *......The beginning may be empty, so we have to step up to the first
  // *      edge where the result is nonzero. We zero the bins which have
  // *      and upper (!) edge which is below the first point of the
  // *      cummulative distribution we are going to project to this
  // *      output histogram binning.
  // *

  ix = 0;
  x = xminn + double(ix+1)*dx;
  if (idebug >= 1) cout << "Start setting initial bins at x=" << x << endl;
  while(x <= xdisn[0]) {
    sigdisf[ix] = sigdisn[0];
    if (idebug >= 1) cout << "   Setting initial bins " << ix << " " << x 
			  << " " << xdisn[1] << " " << sigdisf[ix] << endl;
    ix = ix+1;
    x = xminn + double(ix+1)*dx;
  }
  Int_t ixf = ix;

  if (idebug >= 1)
    cout << "Bins left to loop over:" << ixf << "-" << ixl << endl;

  // *......Also the end (from y to 1.0) often comes before the last edge
  // *      so we have to set the following to 1.0 as well.

  Int_t ix3 = 0; // Problems with initial edge!!!
  for(ix=ixf;ix<ixl;ix++) {
    x = xminn + double(ix)*dx;
    if (x < xdisn[0]) {
      y = 0;
    } else if (x > xdisn[nx3]) {
      y = 1.;
    } else {
      while(xdisn[ix3+1] <= x && ix3 < 2*nbn) {
	ix3 = ix3 + 1;
      }
      if (xdisn[ix3+1]-x > 1.1*dx2) { // Empty bin treatment
	y = sigdisn[ix3+1]; //y = sigdisn[ix3+1];
	if(idebug>=1) cout << "Empty bin treatment " << ix3+1 << " " << y << " " << sigdisn[ix3+2] << " " << sigdisn[ix3-1] << endl;  
      }
      else if (xdisn[ix3+1] > xdisn[ix3]) { // Normal bins
	y = sigdisn[ix3] + (sigdisn[ix3+1]-sigdisn[ix3])
	  *(x-xdisn[ix3])/(xdisn[ix3+1]-xdisn[ix3]);
      } else {  // Is this ever used?
	y = 0;
	cout << "Warning - th1morph: This probably shoudn't happen! " 
	     << endl;
	cout << "Warning - th1morph: Zero slope solving x(y)" << endl;
      }
    }
    sigdisf[ix] = y;
    if (idebug >= 3) {
      cout << ix << ", ix3=" << ix3 << ", xdisn=" << xdisn[ix3] << ", x=" 
	   << x << ", next xdisn=" << xdisn[ix3+1] << endl;
      cout << "   cdf n=" << sigdisn[ix3] << ", y=" << y << ", next point=" 
	   << sigdisn[ix3+1] << endl;
    }
  }

  //......Differentiate interpolated cdf and return renormalized result in 
  //      new histogram. 

  TH1F *morphedhist = (TH1F *)gROOT->FindObject(chname);
  if (morphedhist) delete morphedhist;
  morphedhist = new TH1F(chname,chtitle,nbn,xminn,xmaxn);
 
  Double_t norm = morphedhistnorm;
  // norm1, norm2, wt1, wt2 are computed before the interpolation
  if (norm <= 0) {
    if (norm1 == norm2) {
      norm = norm1;
    } else {
      norm   = wt1*norm1   + wt2*norm2;
    }
  }

  for(Int_t ixx=nbn-1;ixx>-1;ixx--) {
    x = xminn + double(ixx)*dx;
    y =  sigdisf[ixx+1]-sigdisf[ixx];
    if (y<0) cout << "huh??? " << ixx << " " << sigdisf[ixx] << " " << sigdisf[ixx+1] << endl;
    morphedhist->SetBinContent(ixx+1,y*norm);
  }
  
  //......Clean up the temporary arrays we allocated.

  delete sigdis1; delete sigdis2; 
  delete sigdisn; delete xdisn; delete sigdisf;

  //......All done, return the result.
  //morphedhist->Draw("same");
  return(morphedhist);
}
"""

import tempfile
file = tempfile.NamedTemporaryFile(mode='w+',prefix='th1morph',suffix='.C',delete=False)
file.write(morph_code)
file.flush()

import ROOT as r
if not r.gSystem.CompileMacro(file.name):
    raise Exception('Could not compile embedded morphing macro')
file.close()


import numpy as np
def makeBand(n, p, m):
    top = n.Clone(n.GetName()+'_top')
    bot = n.Clone(n.GetName()+'_bottom')
    for f in np.linspace(0,1,1+ 20 ):
        if f<0.5:
            start, end = (m,n)
            tmin, tmax = (0.0, 0.5)
        else:
            start, end = (n,p)
            tmin, tmax = (0.5, 1.0)

        name = n.GetName()+'interpolated'+str(f)
        intermediate = r.th1fmorph(name,name, start, end, tmin, tmax, f)

        for i in xrange(0, n.GetNbinsX() + 1):
            vals = [ getattr(h, 'GetBinContent')(i) for h in (start,intermediate,end) ]
            if max(vals) > top.GetBinContent(i):
                top.SetBinContent(i, max(vals))
            if min(vals) < bot.GetBinContent(i):
                bot.SetBinContent(i, min(vals))

    top.SetMarkerStyle(r.kFullTriangleDown)
    top.SetMarkerColor(r.kBlue)
    bot.SetMarkerStyle(r.kFullTriangleUp)
    bot.SetMarkerColor(r.kRed)
    return (top,bot)

def makeToys():
    print 'Generating toys'
    nominal = r.TH1F('nominal','toy data',100,-5,5)
    nominal.Sumw2()
    plus = nominal.Clone('plus')
    minus = nominal.Clone('minus')

    for i in xrange(100000):
        minus.Fill( r.gRandom.Gaus() )
        nominal.Fill( r.gRandom.Gaus() + 0.75, 0.975 )
        plus.Fill( r.gRandom.Gaus() + 1.50, 0.950 )

    return (nominal,plus,minus)

# import operator
# def checkInterpolation(n,p,m):
#     nInt = r.th1fmorph('nominal_interpolated','Interpolation check',m,p,0.0,1.0,0.5)
#     nInt.SetLineColor(r.kGreen)
#     nInt.SetLineStyle(r.kDashed)
    
#     testspecs =  (('Chi2Test','WWUFOF'),('KolmogorovTest','UO'))
#     tests = [ (test, getattr(n, test)(nInt, args) ) for (test, args) in testspecs ]
#     checks = [ test[1] < 0.05 for test in tests ]
#     if reduce(operator.or_, checks, False):
#         print tests
#         raise Exception("Does not make sense to interpolate " + n.GetName())

#     return nInt

def dryRun():
    histos = makeToys()
    for h in histos:
        h.Print()
        print h.GetName(),h.GetTitle()

    r.gROOT.SetBatch(True)
    c = r.TCanvas('c','bands',800,400)
    c.Divide(2,1)
    c.cd(1)
    nominal, plus, minus = histos
    minus.Draw()
    nominal.Draw('same')
    plus.Draw('same')

    c.cd(2)    
    nominal.Draw()
    envelope = makeBand(*histos)
    for h in envelope:
        h.Draw('same')

    c.SaveAs("dryRun.png")

class FormatDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

def getConfig(name):
    cfg = __import__(name, globals(), locals(), [], -1)
    return (cfg.templates, cfg.variations)

import os
def openFiles(path,tag):
    fin = r.TFile.Open(path,'READ')
    if not fin:
        raise Exception('Could not open '+path+' for reading.')
    fname, fext = os.path.splitext(path)
    outpath = fname+tag+'_envelopes'+fext
    fout = r.TFile.Open(outpath,'RECREATE')
    if not fout:
        raise Exception('Could not open '+outpath+' for writing.')
    return(fin, fout)

def closeFiles(files):
    for file in files:
        file.Close()

def getFromFile(file,name):
    obj = file.Get(name)
    if not obj:
        raise Exception('Could not find '+name+' in input file.')
    return obj


if __name__ == '__main__':

    r.gROOT.Reset()

    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [options] template_spec.py input.root",
                          version="%prog 1.0")
    parser.add_option("-n", "--dry-run",
                      action="store_true",
                      dest="dry_run",
                      default=False,
                      help="Run simple tests on internally generated toy histograms")
    (opts, args) = parser.parse_args()

    if opts.dry_run:
        dryRun()
        import sys
        sys.exit()

    if len(args) != 2:
        parser.print_help()
        parser.error("Please specify two files.")
        
    for path, ext in zip(args, ('.py','.root') ):
        fname, fext = os.path.splitext(path)
        if fext != ext:
            print fname, fext
            parser.print_help()
            parser.error(path+' is not of type '+ext+ '\nPlease specify the files in the right order and type.')


    cfgpath, rootpath = args
    cfgfname, cfgfext = os.path.splitext(cfgpath)
    
    config = (templates,variations) = getConfig(cfgfname)

    files = (fin,fout) = openFiles(rootpath,'_'+os.path.basename(cfgfname))

    import itertools as it
    import string
    fmt = string.Formatter()

    for template in templates:
        if '{syst}' not in template:
                raise Exception('"%s" does not reference {syst} so there is nothing to morph.' % template)
        loci = filter(lambda item: item not in (None, 'syst'), [token[1] for token in fmt.parse(template)])
        variants = it.product(*[ variations[locus] for locus in loci ])
        for variant in variants:
            fmtdict = dict( (k,v) for (k,v) in zip(loci, variant) )
            nosystTemplate = fmt.vformat(template, (), FormatDict(fmtdict))
            histnames = [nosystTemplate.format(syst=syst) for syst in variations['syst'] ]
            print 'Processing '+histnames[0]
            nominal, up, down = histos = [ getFromFile(fin,histname) for histname in histnames ]
            try:
                envelope = makeBand(*histos)
                nominal.Write() 
                for boundary in envelope:
                    boundary.Write()
            except Exception as e:
                print e.message

    closeFiles(files)
