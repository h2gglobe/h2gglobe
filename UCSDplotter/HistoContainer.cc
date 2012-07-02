#include "HistoContainer.h"
#include <utility>
#include <iostream>
#include <cassert>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>


#include <stdexcept>

#define HCDEBUG 1

using namespace std;

//----------------------------------------------------------------------

HistoContainer::HistoContainer() :
  total_scale(1)
{}

//----------------------------------------------------------------------
HistoContainer::HistoContainer(int setTo, std::string setName) {
  setHistVal(setTo);
  setHistNam(setName);
}

//----------------------------------------------------------------------

HistoContainer::~HistoContainer() 
{}

//----------------------------------------------------------------------
void HistoContainer::setHistVal(int setTo) {
  histVal = setTo;
}

//----------------------------------------------------------------------

void HistoContainer::setHistNam(std::string setName) {
  histNam = setName;
}

//----------------------------------------------------------------------
int HistoContainer::getHistVal() {
  return histVal;
}

//----------------------------------------------------------------------

void HistoContainer::setScale(float scale){
  total_scale = scale;
}

//----------------------------------------------------------------------

std::string HistoContainer::fullName(const string &name, int itype, int cat)
{
  string suffix = "itype" + boost::lexical_cast<string>(itype);
  return makeSumName(name, suffix, cat);
}

//----------------------------------------------------------------------

std::string
HistoContainer::makeSumName(const std::string &name, const std::string &suffix, int cat)
{
  return name + "_" + "cat" + boost::lexical_cast<string>(cat) +
      "_" + suffix + histNam;
}

//----------------------------------------------------------------------
map<int, TH1 *> *
HistoContainer::getMapForNameAndCategory(const std::string &name, int category)
{
  map<string, vector< map<int, TH1 *> > >::iterator it = histos.find(name);

  if (it == histos.end())
    {    // a plot with this name was not booked
       cout << "NAME NOT FOUND" << endl;
    return NULL;
    }
  // get category
  if (category < 0 || category >= it->second.size())
    {    // category out of range
      cout << "CATEGORY OUT OF RANGE" << endl;
    return NULL;
    }

  return &it->second[category];
}

//----------------------------------------------------------------------
TH1 *
HistoContainer::createHistogram(const std::string &name, const HistoDef &histoDef)
{
  switch (histoDef.type)
  {
  case HISTO_1D:
    {
      TH1F *histo_temp = new TH1F(name.c_str(), name.c_str(), histoDef.xbins, histoDef.xmin, histoDef.xmax);
      histo_temp->Sumw2();
      histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
      histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
      histo_temp->SetDirectory(0);

      return histo_temp;
    }

  case HISTO_2D:
  {
    TH2F *histo_temp = new TH2F(name.c_str(), name.c_str(),
                                histoDef.xbins,
                                histoDef.xmin,
                                histoDef.xmax,

                                histoDef.ybins,
                                histoDef.ymin,
                                histoDef.ymax);

    histo_temp->Sumw2();
    histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
    histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
    histo_temp->SetDirectory(0);

    return histo_temp;
  }

  case HISTO_PROFILE:
  {
    TProfile *histo_temp = new TProfile(name.c_str(), name.c_str(),
      histoDef.xbins,
      histoDef.xmin,
      histoDef.xmax,

      // histoDef.ybins,
      histoDef.ymin,
      histoDef.ymax);

      histo_temp->Sumw2();
      histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
      histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
      histo_temp->SetDirectory(0);

      return histo_temp;
  }

  default:
    cerr << "internal error on line " << __LINE__ << " in file " << __FILE__ << endl;
    exit(1);
  }
}

//----------------------------------------------------------------------

TH1F*
HistoContainer::Get1Dhisto(const std::string &name, int itype, int category)
{
  map<int, TH1 *> *histoMap = getMapForNameAndCategory(name, category);
  if (histoMap == NULL)
    // something was wrong
    return NULL;

  // try to find itype
  map<int, TH1 *>::iterator it2 = histoMap->find(itype);

  if (it2 != histoMap->end())
  {
    assert(it2->second != NULL);
    TH1F *retval = dynamic_cast<TH1F*>(it2->second);
    assert(retval != NULL);
    return retval;
  }

  // book a new histogram for this itype
  std::string modName = fullName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TH1F *histo_temp = dynamic_cast<TH1F*>(createHistogram(modName, histoDef));
  assert(histo_temp != NULL);
  (*histoMap)[itype] = histo_temp;

  return histo_temp;
}

//----------------------------------------------------------------------

TH2F *
HistoContainer::Get2Dhisto(const std::string &name, int itype, int category)
{
  map<int, TH1 *> *histoMap = getMapForNameAndCategory(name, category);
  if (histoMap == NULL)
    // something was wrong
    return NULL;

  // try to find itype
  map<int, TH1 *>::iterator it2 = histoMap->find(itype);

  if (it2 != histoMap->end())
  {
    assert(it2->second != NULL);
    TH2F *retval = dynamic_cast<TH2F*>(it2->second);
    assert(retval != NULL);
    return retval;
  }

  // book a new histogram for this itype
  std::string modName = fullName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TH2F *histo_temp = dynamic_cast<TH2F*>(createHistogram(modName, histoDef));
  assert(histo_temp != NULL);
  (*histoMap)[itype] = histo_temp;

  return histo_temp;
}

//----------------------------------------------------------------------
TProfile *
HistoContainer::GetProfileHisto(const std::string &name, int itype, int category)
{
  map<int, TH1 *> *histoMap = getMapForNameAndCategory(name, category);
  if (histoMap == NULL)
    // something was wrong
    return NULL;

  // try to find itype
  map<int, TH1 *>::iterator it2 = histoMap->find(itype);

  if (it2 != histoMap->end())
  {
    assert(it2->second != NULL);
    TProfile *retval = dynamic_cast<TProfile*>(it2->second);
    assert(retval != NULL);
    return retval;
  }

  // book a new histogram for this itype
  std::string modName = fullName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TProfile *histo_temp = dynamic_cast<TProfile*>(createHistogram(modName, histoDef));
  assert(histo_temp != NULL);
  (*histoMap)[itype] = histo_temp;

  return histo_temp;
}

//----------------------------------------------------------------------
void
HistoContainer::Add(const std::string &name, const HistoDef &histoDef, int categories)
{
  // ensure there is not already a histogram with this name
  assert(histos[name].size() == 0);
  names.push_back(name);

  for (int i=0; i<categories; i++) 
  {
    // add empty containers to remember the number of categories
    histos[name].push_back(map<int, TH1*>());
  }

  histoDefs[name] = histoDef;

}

//----------------------------------------------------------------------

void HistoContainer::Add(const std::string &name, const std::string &xaxis, const std::string &yaxis, int categories,int bins, float xmin, float xmax)
{
  // store parameters
  HistoDef histoDef;

  histoDef.type = HISTO_1D;

  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = bins;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  // leave the y axis parameters unset

  Add(name, histoDef, categories);
}

//----------------------------------------------------------------------

void HistoContainer::Add(const std::string &name, const std::string &xaxis, const std::string &yaxis, int categories, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) 
{
  // store parameters
  HistoDef histoDef;

  histoDef.type = HISTO_2D;

  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = binsx;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  histoDef.ybins = binsy;
  histoDef.ymin = ymin;
  histoDef.ymax = ymax;

  Add(name, histoDef, categories);

}

//----------------------------------------------------------------------

void HistoContainer::Add(char* name, char* xaxis, char* yaxis, int categories, int binsx, 
			 float xmin, float xmax,
			 float ymin, float ymax)
{
  // store parameters
  HistoDef histoDef;

  histoDef.type = HISTO_PROFILE;

  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = binsx;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  // histoDef.ybins = binsy;
  histoDef.ymin = ymin;
  histoDef.ymax = ymax;
  
  Add(name, histoDef, categories);
} 

//----------------------------------------------------------------------

void HistoContainer::Fill(const std::string &name, int itype, int category, float value, float weight) {

  // std::string modName = ModifiedName(name);

  TH1F* histo = Get1Dhisto(name, itype, category);
  if (histo == NULL)
  {
    // could be that the plot was not booked or the category number is out of the
    // specified range
    cerr << "ERROR: trying to fill non-existing plot '" << name << "'" << " for category " << category << endl;
    return;
  }

  histo->Fill(value, total_scale * weight);
}

//----------------------------------------------------------------------

void HistoContainer::Fill2D(const std::string &name, int itype, int category, float valuex, float valuey, float weight) 
{ 
  // first search the 2D histograms
  TH2F* histo = Get2Dhisto(name, itype, category);
  if (histo != NULL)
    {
      histo->Fill(valuex, valuey, total_scale * weight);
      return;
    }

  // then the profile histograms
  TProfile *prof = GetProfileHisto(name, itype, category);

  if (prof == NULL)
  {
    // could be that the plot was not booked or the category number is out of the
    // specified range
    cerr << "ERROR: trying to fill non-existing plot '" << name << "'" << " for category " << category << endl;
    return;
  }

  prof->Fill(valuex, valuey, total_scale * weight);

//  //std::string modName = ModifiedName(name);
//  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(name);
//  if (itp != hp.end()) {
//    (itp->second)[category].Fill(valuex, valuey, total_scale*weight);
//    return;
//  }
//
//  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(std::string(name));
//  if (it2 != h2.end()) {
//    (it2->second)[category].Fill(valuex, valuey, total_scale*weight);
//    return;
//  }
//
 //  cerr << "ERROR: trying to fill non-existing plot '" << name << "'" << " for category " << category << endl;
}

//----------------------------------------------------------------------

void HistoContainer::Save() {
  map<string, vector<map<int, TH1 *> > >::iterator it;
  for (it = histos.begin(); it != histos.end(); ++it)
  {
    for (unsigned int i = 0; i < (it->second).size(); i++)
    {
      map<int, TH1*> &tmp = (it->second)[i];
      for (map<int, TH1 *>::iterator it2 = tmp.begin(); it2 != tmp.end(); ++it2)
      {
        it2->second->Write();
      } // loop over itypes
    } // loop over categories
  } // loop over plot names


}

//----------------------------------------------------------------------

// int HistoContainer::getDimension(int n) {
// 
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end())
//     return 1;
// 
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end())
//     return 2;
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end())
//     return 3;
// }

//----------------------------------------------------------------------

int HistoContainer::getNumCategories(const std::string &name)
{
   map<string, vector<map<int, TH1 *> > >::iterator it = histos.find(name);
   if (it != histos.end())
     return it->second.size();

   return -1;
 }

//----------------------------------------------------------------------

// int HistoContainer::nbins(int n, bool isX) {
//   
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end()) {
//     if (isX)
//       return (it->second)[0].GetNbinsX();
//     else
//       return -1;
//   }
//   
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end()) {
//     if (isX)
//       return (it2->second)[0].GetNbinsX();
//     else
//       return (it2->second)[0].GetNbinsY();
//   }
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end()) {
//     if (isX)
//       return (itp->second)[0].GetNbinsX();
//     else
//       return -1;
//   }
// 
//   return -1;
// }

//----------------------------------------------------------------------

// float HistoContainer::max(int n, bool isX) {
//   
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end()) {
//     if (isX)
//       return (it->second)[0].GetXaxis()->GetXmax();
//     else
//       return (it->second)[0].GetYaxis()->GetXmax();
//   }
//   
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end()) {
//     if (isX)
//       return (it2->second)[0].GetXaxis()->GetXmax();
//     else
//       return (it2->second)[0].GetYaxis()->GetXmax();
//   }
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end()) {
//     if (isX)
//       return (itp->second)[0].GetXaxis()->GetXmax();
//     else
//       return (itp->second)[0].GetYaxis()->GetXmax();
//   }
// 
//   return -1;
// }
// 
// //----------------------------------------------------------------------
// 
// float HistoContainer::min(int n, bool isX) {
//   
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end()) {
//     if (isX)
//       return (it->second)[0].GetXaxis()->GetXmin();
//     else
//       return (it->second)[0].GetYaxis()->GetXmin();
//   }
//   
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end()) {
//     if (isX)
//       return (it2->second)[0].GetXaxis()->GetXmin();
//     else
//       return (it2->second)[0].GetYaxis()->GetXmin();
//   }
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end()) {
//     if (isX)
//       return (itp->second)[0].GetXaxis()->GetXmin();
//     else
//       return (itp->second)[0].GetYaxis()->GetXmin();
//   }
// 
//   return -1;
// }

//----------------------------------------------------------------------

// std::string HistoContainer::axisName(int n, bool isX) {
//   
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end()) {
//     if (isX)
//       return (it->second)[0].GetXaxis()->GetTitle();
//     else
//       return (it->second)[0].GetYaxis()->GetTitle();
//   }
//   
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end()) {
//     if (isX)
//       return (it2->second)[0].GetXaxis()->GetTitle();
//     else
//       return (it2->second)[0].GetYaxis()->GetTitle();
//   }
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end()) {
//     if (isX)
//       return (itp->second)[0].GetXaxis()->GetTitle();
//     else
//       return (itp->second)[0].GetYaxis()->GetTitle();
//   }
// 
//   return "";
// }

//----------------------------------------------------------------------

std::vector<int>
HistoContainer::getItypes(const std::string &plotName, int cat)
{
  assert(cat >= 0);
  {
    map<string, vector<map<int, TH1 *> > >::iterator it1 = histos.find(plotName);
    if (it1 != histos.end())
    {
      vector<map<int, TH1 *> > &tmp = it1->second;
      assert(cat < tmp.size());

      vector<int> retval;
      for (map<int, TH1*>::const_iterator it = tmp[cat].begin(); it != tmp[cat].end(); ++it)
        retval.push_back(it->first);

      return retval;
    }
  }

  throw std::runtime_error("could not find plot '" + plotName + "' cat " + boost::lexical_cast<string>(cat));

}
//----------------------------------------------------------------------

TH1*
HistoContainer::makeHistoSum(const std::string &name, const std::string &sumName, const std::vector<int> &itypes, int cat)
{
  string fullName = this->makeSumName(name, sumName, cat);

  // check what type of histogram the histograms to be summed are
  // (take the type from the first histogram)
  //
  // we also have use cases where we e.g. only have data but no MC
  //assert(itypes.size() > 0);

  assert(this->histoDefs.count(name) > 0);

  // book the sum histogram
  TH1 *retval = this->createHistogram(fullName, histoDefs[name]);

  BOOST_FOREACH(int itype, itypes)
  {
    map<int, TH1*>::iterator it = histos[name][cat].find(itype);
    if (it == histos[name][cat].end())
      // itype never appeared for this category / plot
      continue;

    retval->Add(it->second);

  } // loop over all itypes

  // should have found at least one histogram
  return retval;
}

//----------------------------------------------------------------------

bool
HistoContainer::hasCategory(const std::string &plotName, int cat)
{
  if (cat < 0 || cat >= histos[plotName].size())
    // out of range anyway
    return false;

  // check whether something was created/filled for this category
  return histos[plotName][cat].size() != 0;
}

//----------------------------------------------------------------------
