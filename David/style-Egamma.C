/***********************
  To use this:
  root[0] .L style-Egamma.C 
  root[1] setEgammaStyle()
***********************/
#include "TStyle.h"

void egammaGrid(bool gridOn) {
  egammaStyle->SetPadGridX(gridOn);
  egammaStyle->SetPadGridY(gridOn);
}

void fixOverlay() {
  gPad->RedrawAxis();
}

void egammaLegend() {

  // leg->SetLegendBorderSize(0);
  // leg->SetLegendFillColor(0);
  // leg->SetLegendFont(142);
  // leg->SetTextSize(0.05);
  // leg->SetLineColor(1);
  // leg->SetLineStyle(1);
  // leg->SetLineWidth(1);
  // leg->SetFillColor(0);
  // leg->SetFillStyle(1001);
}

void egammaPrelim() {

  // TPaveText *pt = new TPaveText(0.15,0.9,0.77,0.98,"brNDC");
  // pt->SetBorderSize(0);
  // pt->SetFillStyle(0);
  // pt->SetTextAlign(12);
  // pt->SetTextFont(132);
  // pt->SetTextSize(0.04);
  // TText *text = pt->AddText("CMS Preliminary #sqrt{s} = 8 TeV,  L = 19.6 fb^{-1}");
}

void setEgammaStyle() {
  TStyle *egammaStyle = new TStyle("egammaStyle","Style for P-TDR");

  // For the palette:
  egammaStyle->SetPalette(1);

  // For the grid
  egammaStyle->SetPadGridX(1);
  egammaStyle->SetPadGridY(1);

  // For the canvas: 
  egammaStyle->SetCanvasBorderMode(0);
  egammaStyle->SetCanvasBorderSize(2);     
  egammaStyle->SetCanvasColor(0);
  //egammaStyle->SetCanvasDefH(600);  
  //egammaStyle->SetCanvasDefW(600);  
  egammaStyle->SetCanvasDefX(0);    
  egammaStyle->SetCanvasDefY(0);

  // For the pad:
  egammaStyle->SetPadBorderMode(0);
  egammaStyle->SetPadColor(kWhite);
  egammaStyle->SetPadGridX(false);
  egammaStyle->SetPadGridY(false);
  egammaStyle->SetGridColor(0);
  egammaStyle->SetGridStyle(3);
  egammaStyle->SetGridWidth(1);

  // For the frame:
  egammaStyle->SetFrameBorderMode(0);
  egammaStyle->SetFrameBorderSize(1);
  egammaStyle->SetFrameFillColor(0);
  egammaStyle->SetFrameFillStyle(0);
  egammaStyle->SetFrameLineColor(1);
  egammaStyle->SetFrameLineStyle(1);
  egammaStyle->SetFrameLineWidth(1);

  // For the histos:
  egammaStyle->SetHistLineColor(1);
  egammaStyle->SetHistLineStyle(0);
  egammaStyle->SetHistLineWidth(2);      // chiara: era 1
  egammaStyle->SetEndErrorSize(2);
  egammaStyle->SetErrorX(0.);

  // for the graphs
  egammaStyle->SetMarkerColor(1);
  egammaStyle->SetMarkerStyle(20);
  egammaStyle->SetMarkerSize(1.1);        // chiara: era 1.2

  // For the fit/function:
  egammaStyle->SetOptFit(1);
  egammaStyle->SetFitFormat("5.4g");
  egammaStyle->SetFuncColor(2);
  egammaStyle->SetFuncStyle(1);
  egammaStyle->SetFuncWidth(2);          // chiara: era 1

  //For the date:
  egammaStyle->SetOptDate(0);

  // For the statistics box:
  egammaStyle->SetOptFile(0);
  egammaStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  egammaStyle->SetStatColor(kWhite);
  egammaStyle->SetStatFont(42);
  egammaStyle->SetStatFontSize(0.025);
  egammaStyle->SetStatTextColor(1);
  egammaStyle->SetStatFormat("6.4g");
  egammaStyle->SetStatBorderSize(1);
  egammaStyle->SetStatH(0.1);
  egammaStyle->SetStatW(0.15);

  // For the Global title:
  egammaStyle->SetOptTitle(1);    // 0=No Title
  egammaStyle->SetTitleFont(42);
  egammaStyle->SetTitleColor(1);
  egammaStyle->SetTitleTextColor(1);
  egammaStyle->SetTitleFillColor(10);
  egammaStyle->SetTitleFontSize(0.05);      

  // For the axis titles:
  egammaStyle->SetTitleColor(1, "XYZ");
  egammaStyle->SetTitleFont(42, "XYZ");
  egammaStyle->SetTitleSize(0.05, "XYZ");        // chiara: era 0.06
  // egammaStyle->SetTitleXOffset(0.9);          // chiara
  // egammaStyle->SetTitleYOffset(1.25);         // chiara

  // For the axis labels:
  egammaStyle->SetLabelColor(1, "XYZ");
  egammaStyle->SetLabelFont(42, "XYZ");
  egammaStyle->SetLabelOffset(0.007, "XYZ");    // chiara 
  egammaStyle->SetLabelSize(0.04, "XYZ");       // chiara: era 0.05

  // For the axis:
  egammaStyle->SetAxisColor(1, "XYZ");
  egammaStyle->SetStripDecimals(kTRUE);
  egammaStyle->SetTickLength(0.03, "XYZ");
  egammaStyle->SetNdivisions(510, "XYZ");
  egammaStyle->SetPadTickX(0);  // 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom
  egammaStyle->SetPadTickY(1);

  // Change for log plots:
  egammaStyle->SetOptLogx(0);
  egammaStyle->SetOptLogy(0);
  egammaStyle->SetOptLogz(0);

  // Postscript options:
  egammaStyle->SetPaperSize(20.,20.);

  egammaStyle->cd();
}
