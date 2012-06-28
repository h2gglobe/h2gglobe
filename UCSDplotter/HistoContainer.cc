#include "HistoContainer.h"
#include <utility>
#include <iostream>
#include <cassert>

#include <boost/foreach.hpp>

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

std::string HistoContainer::ModifiedName(const char* name, int itype, int cat) 
{
  char* modName= new char[500];
  sprintf(modName,"%s_cat%d_itype%d%s",name, cat, itype, histNam.c_str());
  return modName;
}

//----------------------------------------------------------------------
TH1F *
HistoContainer::Get1Dhisto(const std::string &name, int itype, int category)
{
  map<string, vector< map<int, TH1F *> > >::iterator it = h1.find(name);

  if (it == h1.end()) 
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
  map<int, TH1F *> &histoMap = it->second[category];

  // try to find itype
  map<int, TH1F *>::iterator it2 = histoMap.find(itype);

  if (it2 != histoMap.end())
    return it2->second;

  // book a new histogram for this itype
  std::string modName = ModifiedName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TH1F *histo_temp = new TH1F(modName.c_str(), modName.c_str(), 
                              histoDef.xbins, 
                              histoDef.xmin, 
                              histoDef.xmax);

  histoMap[itype] = histo_temp;

  histo_temp->Sumw2();
  histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
  histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
  histo_temp->SetDirectory(0);

  return histo_temp;
}

//----------------------------------------------------------------------

TH2F *
HistoContainer::Get2Dhisto(const std::string &name, int itype, int category)
{
  map<string, vector< map<int, TH2F *> > >::iterator it = h2.find(name);

  if (it == h2.end()) 
    // a plot with this name was not booked 
    return NULL;

  // get category
  if (category < 0 || category >= it->second.size())
    // category out of range
    return NULL;

  map<int, TH2F *> &histoMap = it->second[category];

  // try to find itype
  map<int, TH2F *>::iterator it2 = histoMap.find(itype);

  if (it2 != histoMap.end())
    return it2->second;

  // book a new histogram for this itype
  std::string modName = ModifiedName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TH2F *histo_temp = new TH2F(modName.c_str(), modName.c_str(), 
                              histoDef.xbins, 
                              histoDef.xmin, 
                              histoDef.xmax,

                              histoDef.ybins, 
                              histoDef.ymin, 
                              histoDef.ymax);

  histoMap[itype] = histo_temp;

  histo_temp->Sumw2();
  histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
  histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
  histo_temp->SetDirectory(0);

  return histo_temp;
}

//----------------------------------------------------------------------
TProfile *
HistoContainer::GetProfileHisto(const std::string &name, int itype, int category)
{
  map<string, vector< map<int, TProfile *> > >::iterator it = hp.find(name);

  if (it == hp.end()) 
    // a plot with this name was not booked 
    return NULL;

  // get category
  if (category < 0 || category >= it->second.size())
    // category out of range
    return NULL;

  map<int, TProfile *> &histoMap = it->second[category];

  // try to find itype
  map<int, TProfile *>::iterator it2 = histoMap.find(itype);

  if (it2 != histoMap.end())
    return it2->second;

  // book a new histogram for this itype
  std::string modName = ModifiedName(name.c_str(), itype, category);

  // get parameters for histogram booking 
  HistoDef histoDef = histoDefs[name];

  TProfile *histo_temp = new TProfile(modName.c_str(), modName.c_str(), 
                              histoDef.xbins, 
                              histoDef.xmin, 
                              histoDef.xmax,

                              // histoDef.ybins, 
                              histoDef.ymin, 
                              histoDef.ymax);

  histoMap[itype] = histo_temp;

  histo_temp->Sumw2();
  histo_temp->GetXaxis()->SetTitle(histoDef.xtitle.c_str());
  histo_temp->GetYaxis()->SetTitle(histoDef.ytitle.c_str());
  histo_temp->SetDirectory(0);

  return histo_temp;
}

//----------------------------------------------------------------------

void HistoContainer::Add(const std::string &name, const std::string &xaxis, const std::string &yaxis, int categories,int bins, float xmin, float xmax) 
{
  // ensure there is not already a histogram with this name
  assert(h1[name].size() == 0);

  for (int i=0; i<categories; i++) 
  {
    // add empty containers to remember the number of categories
    h1[name].push_back(map<int, TH1F*>());
  }

  names.push_back(name);
  // h1[std::string(name)] = temp;

  // store parameters
  HistoDef histoDef;
  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = bins;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  // leave the y axis parameters unset

  histoDefs[name] = histoDef;
}

//----------------------------------------------------------------------

void HistoContainer::Add(const std::string &name, const std::string &xaxis, const std::string &yaxis, int categories, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) 
{

  // ensure there is not already a histogram with this name
  assert(h2[name].size() == 0);

  for (int i=0; i<categories; i++) 
  {
    // add empty containers to remember the number of categories
    h2[name].push_back(map<int, TH2F*>());
  }

  names.push_back(name);
  // h1[std::string(name)] = temp;

  // store parameters
  HistoDef histoDef;
  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = binsx;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  histoDef.ybins = binsy;
  histoDef.ymin = ymin;
  histoDef.ymax = ymax;
  
  histoDefs[name] = histoDef;

}

//----------------------------------------------------------------------

void HistoContainer::Add(char* name, char* xaxis, char* yaxis, int categories, int binsx, 
			 float xmin, float xmax,
			 float ymin, float ymax) {



  // ensure there is not already a histogram with this name
  assert(hp[name].size() == 0);

  for (int i=0; i<categories; i++) 
  {
    // add empty containers to remember the number of categories
    hp[name].push_back(map<int, TProfile*>());
  }

  names.push_back(name);
  // h1[std::string(name)] = temp;

  // store parameters
  HistoDef histoDef;
  histoDef.xtitle = xaxis;
  histoDef.ytitle = yaxis;
  histoDef.xbins = binsx;
  histoDef.xmin = xmin;
  histoDef.xmax = xmax;

  // histoDef.ybins = binsy;
  histoDef.ymin = ymin;
  histoDef.ymax = ymax;
  
  histoDefs[name] = histoDef;

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
  // TH1
  {
    map<string, vector<map<int, TH1F *> > >::iterator it;
    for (it = h1.begin(); it != h1.end(); ++it) 
      {
        for (unsigned int i=0; i<(it->second).size(); i++) 
          {
            map<int, TH1F* > &tmp = (it->second)[i];
            for (map<int, TH1F *>::iterator it2 = tmp.begin(); it2 != tmp.end(); ++it2)
              {
                it2->second->Write();
              } // loop over itypes
          } // loop over categories
      } // loop over plot names
  }

  // TH2
  {
    map<string, vector<map<int, TH2F *> > >::iterator it;
    for (it = h2.begin(); it != h2.end(); ++it) 
      {
        for (unsigned int i=0; i<(it->second).size(); i++) 
          {
            map<int, TH2F* > &tmp = (it->second)[i];
            for (map<int, TH2F *>::iterator it2 = tmp.begin(); it2 != tmp.end(); ++it2)
              {
                it2->second->Write();
              } // loop over itypes
          } // loop over categories
      } // loop over plot names
  }

  // TProfile
  {
    map<string, vector<map<int, TProfile *> > >::iterator it;
    for (it = hp.begin(); it != hp.end(); ++it) 
      {
        for (unsigned int i=0; i<(it->second).size(); i++) 
          {
            map<int, TProfile* > &tmp = (it->second)[i];
            for (map<int, TProfile *>::iterator it2 = tmp.begin(); it2 != tmp.end(); ++it2)
              {
                it2->second->Write();
              } // loop over itypes
          } // loop over categories
      } // loop over plot names
  }


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

// int HistoContainer::ncat(int n) {
//   
//   std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
//   if (it != h1.end())
//     return it->second.size();
//   
//   std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
//   if (it2 != h2.end())
//     return it2->second.size();
//   
//   std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
//   if (itp != hp.end())
//     return itp->second.size(); 
//   
//   return -1;
// }

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
