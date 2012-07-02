#include "Plotting.h"

#include <vector>
#include <boost/foreach.hpp>
#include <cassert>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <iostream>


using namespace std;

//----------------------------------------------------------------------
Plotting::Plotting(HistoContainer *_histoContainer) : 
  histoContainer(_histoContainer)
{
  gStyle->SetOptStat(0);
}

//----------------------------------------------------------------------

void
Plotting::plotAll()
{
  vector<string> allNames = histoContainer->getAllNames();

  BOOST_FOREACH(string name, allNames)
  {
    plot(name);
  }
}

//----------------------------------------------------------------------

void 
Plotting::plot(const std::string &plotName)
{
  // find all categories for this plot
  int numCategories = this->histoContainer->getNumCategories(plotName);
  assert(numCategories >= 0);

  for (int cat = 0; cat < numCategories; ++cat)
  {
    // check whether this category actually exists
    if (! histoContainer->hasCategory(plotName, cat))
      continue;

    plot(plotName, cat);
  }
}

//----------------------------------------------------------------------

bool
Plotting::isSignal(int itype)
{
  // TODO: should eventually use a set to store the selected signal itypes
  return std::find(signalItypes.begin(), signalItypes.end(), itype) != signalItypes.end();
}

//----------------------------------------------------------------------
void
Plotting::plot(const std::string &plotName, int cat)
{
  // loop over all itypes
  vector<int> itypes = this->histoContainer->getItypes(plotName, cat);

  // split itypes into data, signal and background
  vector<int> sigItypes, bgItypes, dataItypes;

  BOOST_FOREACH(int itype, itypes)
  {
    if (isSignal(itype))
      sigItypes.push_back(itype);
    else if (isBackground(itype))
      bgItypes.push_back(itype);
    else if (isData(itype))
      dataItypes.push_back(itype);
  } // loop over all itypes

  cout << dataItypes.size() << " itypes for data" << endl;

  // calculate the sum histograms for all three
  TH1 *sig = histoContainer->makeHistoSum(plotName,"sig",sigItypes, cat);
  TH1 *bg = histoContainer->makeHistoSum(plotName,"bg",bgItypes, cat);
  TH1 *data = histoContainer->makeHistoSum(plotName,"data",dataItypes, cat);

  // TODO: what happens for non-1D histograms here ?
  new TCanvas;

  double maxY = max(max(sig->GetMaximum(), bg->GetMaximum()), data->GetMaximum());

  maxY *= 1.1;

  sig->SetMaximum(maxY);
  bg->SetMaximum(maxY);
  data->SetMaximum(maxY);

  sig->SetLineColor(kRed);
  bg->SetLineColor(kBlue);

  sig->SetLineWidth(2);
  bg->SetLineWidth(2);

  data->SetMarkerStyle(20);

  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);

  //----------
  data->Draw("E");
  legend->AddEntry(data,("data " + boost::str(boost::format("%.1f") % data->GetSumOfWeights())).c_str());

  bg->Draw("same,hist");
  legend->AddEntry(bg,("BG " + boost::str(boost::format("%.1f") % bg->GetSumOfWeights())).c_str());


  sig->Draw("same,hist");
  legend->AddEntry(sig,("sig " + boost::str(boost::format("%.2f") % sig->GetSumOfWeights())).c_str());


  //----------

  legend->SetFillColor(kWhite);
  legend->Draw();

  //----------

  string outputName = plotName + "_cat" + boost::lexical_cast<string>(cat) + ".png";
  gPad->SaveAs(outputName.c_str());

}

//----------------------------------------------------------------------

void
Plotting::addSignalItype(int itype)
{
  this->signalItypes.push_back(itype);
}


//----------------------------------------------------------------------
