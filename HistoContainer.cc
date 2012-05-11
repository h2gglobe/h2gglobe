#include "HistoContainer.h"
#include <utility>
#include <iostream>

#define HCDEBUG 1


HistoContainer::HistoContainer() {}

HistoContainer::HistoContainer(int setTo, std::string setName) {
  setHistVal(setTo);
  setHistNam(setName);
}

HistoContainer::~HistoContainer() 
{}

void HistoContainer::setHistVal(int setTo) {
  histVal = setTo;
}

void HistoContainer::setHistNam(std::string setName) {
  histNam = setName;
}

int HistoContainer::getHistVal() {
  return histVal;
}

void HistoContainer::setScale(float scale){
  total_scale = scale;
}

std::string HistoContainer::ModifiedName(char* name, int i) {
  char* modName= new char[500];
  sprintf(modName,"%s_cat%d_%s",name,i,histNam.c_str());
  std::string output(modName);
  return output;
}

void HistoContainer::Add(char* name, char* xaxis, char* yaxis, int categories,int bins, float xmin, float xmax) {

  std::vector<TH1F> temp;
  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i);
    TH1F histo_temp(modName.c_str(), modName.c_str(), bins, xmin, xmax);
    histo_temp.GetXaxis()->SetTitle(xaxis);
    histo_temp.GetYaxis()->SetTitle(yaxis);
    histo_temp.SetDirectory(0);
    temp.push_back(histo_temp);
    temp.back().SetDirectory(0);
    //temp.push_back(TH1F(modName.c_str(), modName.c_str(), bins, xmin, xmax));
  }
  names.push_back(name);
  h1[std::string(name)] = temp;
}


void HistoContainer::Add(char* name, char* xaxis, char* yaxis, int categories, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
  
  std::vector<TH2F> temp;
  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i); 
    TH2F histo_temp(modName.c_str(), modName.c_str(), binsx, xmin, xmax, binsy, ymin, ymax);
    histo_temp.GetXaxis()->SetTitle(xaxis);
    histo_temp.GetYaxis()->SetTitle(yaxis);
    temp.push_back(histo_temp);
    //temp.push_back(TH2F(modName.c_str(), modName.c_str(), binsx, xmin, xmax, binsy, ymin, ymax));
  }
  names.push_back(name);
  h2[std::string(name)] = temp;
}

void HistoContainer::Add(char* name, char* xaxis, char* yaxis, int categories, int binsx, 
			 float xmin, float xmax,
			 float ymin, float ymax) {

  std::vector<TProfile> temp;
  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i);
    TProfile histo_temp(modName.c_str(), modName.c_str(), binsx, xmin, xmax, ymin, ymax);
    histo_temp.GetXaxis()->SetTitle(xaxis);
    histo_temp.GetYaxis()->SetTitle(yaxis);
    temp.push_back(histo_temp);
    //temp.push_back(TProfile(modName.c_str(), modName.c_str(), binsx, xmin, xmax, ymin, ymax));
  }
  names.push_back(name);
  hp[std::string(name)] = temp;
} 

void HistoContainer::Fill(std::string name, int category, float value) {
  Fill(name, category, value, 1.0);
}

void HistoContainer::Fill(std::string name, int category, float value, float weight) {

  //std::string modName = ModifiedName(name);
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(std::string(name));
  /// if( category == 0 )std::cerr << " HistoContainer::Fill " << total_scale << " " << weight << std::endl;
  if (it != h1.end()) {
    (it->second)[category].Fill(value, total_scale*weight);
    return;
  }
}

void HistoContainer::Fill2D(std::string name, int category, float valuex, float valuey) { 
  Fill2D(name, category, valuex, valuey, 1.0);
}

void HistoContainer::Fill2D(std::string name, int category, float valuex, float valuey, float weight) { 

  //std::string modName = ModifiedName(name);
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(std::string(name));
  if (itp != hp.end()) {
    (itp->second)[category].Fill(valuex, valuey, total_scale*weight);
    return;
  }

  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(std::string(name));
  if (it2 != h2.end()) {
    (it2->second)[category].Fill(valuex, valuey, total_scale*weight);
    return;
  }
}

void HistoContainer::Save() {
  std::map<std::string, std::vector<TH1F> >::iterator it;
  for (it = h1.begin(); it != h1.end(); ++it) {
    for (unsigned int i=0; i<(it->second).size(); i++) {
      (it->second)[i].Write();
    }
  }

  std::map<std::string, std::vector<TH2F> >::iterator it2;
  for (it2 = h2.begin(); it2 != h2.end(); ++it2) {
    for (unsigned int i=0; i<(it2->second).size(); i++) {
      (it2->second)[i].Write(); 
    }
  }

  std::map<std::string, std::vector<TProfile> >::iterator itp;
  for (itp = hp.begin(); itp != hp.end(); ++itp) {
    for (unsigned int i=0; i<(itp->second).size(); i++) {
      (itp->second)[i].Write();
    }
  }
}

int HistoContainer::getDimension(int n) {

  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end())
    return 1;

  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end())
    return 2;
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end())
    return 3;
}

int HistoContainer::ncat(int n) {
  
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end())
    return it->second.size();
  
  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end())
    return it2->second.size();
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end())
    return itp->second.size(); 
  
  return -1;
}

int HistoContainer::nbins(int n, bool isX) {
  
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end()) {
    if (isX)
      return (it->second)[0].GetNbinsX();
    else
      return -1;
  }
  
  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end()) {
    if (isX)
      return (it2->second)[0].GetNbinsX();
    else
      return (it2->second)[0].GetNbinsY();
  }
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end()) {
    if (isX)
      return (itp->second)[0].GetNbinsX();
    else
      return -1;
  }

  return -1;
}

float HistoContainer::max(int n, bool isX) {
  
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end()) {
    if (isX)
      return (it->second)[0].GetXaxis()->GetXmax();
    else
      return (it->second)[0].GetYaxis()->GetXmax();
  }
  
  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end()) {
    if (isX)
      return (it2->second)[0].GetXaxis()->GetXmax();
    else
      return (it2->second)[0].GetYaxis()->GetXmax();
  }
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end()) {
    if (isX)
      return (itp->second)[0].GetXaxis()->GetXmax();
    else
      return (itp->second)[0].GetYaxis()->GetXmax();
  }

  return -1;
}

float HistoContainer::min(int n, bool isX) {
  
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end()) {
    if (isX)
      return (it->second)[0].GetXaxis()->GetXmin();
    else
      return (it->second)[0].GetYaxis()->GetXmin();
  }
  
  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end()) {
    if (isX)
      return (it2->second)[0].GetXaxis()->GetXmin();
    else
      return (it2->second)[0].GetYaxis()->GetXmin();
  }
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end()) {
    if (isX)
      return (itp->second)[0].GetXaxis()->GetXmin();
    else
      return (itp->second)[0].GetYaxis()->GetXmin();
  }

  return -1;
}

std::string HistoContainer::axisName(int n, bool isX) {
  
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(names[n]);
  if (it != h1.end()) {
    if (isX)
      return (it->second)[0].GetXaxis()->GetTitle();
    else
      return (it->second)[0].GetYaxis()->GetTitle();
  }
  
  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(names[n]);
  if (it2 != h2.end()) {
    if (isX)
      return (it2->second)[0].GetXaxis()->GetTitle();
    else
      return (it2->second)[0].GetYaxis()->GetTitle();
  }
  
  std::map<std::string, std::vector<TProfile> >::iterator itp = hp.find(names[n]);
  if (itp != hp.end()) {
    if (isX)
      return (itp->second)[0].GetXaxis()->GetTitle();
    else
      return (itp->second)[0].GetYaxis()->GetTitle();
  }

  return "";
}
