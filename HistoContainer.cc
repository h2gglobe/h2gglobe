#include "HistoContainer.h"
#include <utility>
#include <iostream>

HistoContainer::HistoContainer() 
{}

HistoContainer::~HistoContainer() 
{}

void HistoContainer::Add(const char* name, int bins, float xmin, float xmax) {
  TH1F temp(name, name, bins, xmin, xmax);
  h1[name] = temp;
}


void HistoContainer::Add(const char* name, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
		       
  TH2F temp(name, name, binsx, xmin, xmax, binsy, ymin, ymax);
  h2[name] = temp;
}

void HistoContainer::Add(const char* name, int binsx, float xmin, float xmax,
			 float ymin, float ymax) {
		       
  TProfile temp(name, name, binsx, xmin, xmax, ymin, ymax);
  hp[name] = temp;
} 

void HistoContainer::Fill(const char* name, float value) {

  std::map<const char*, TH1F>::const_iterator it = h1.find(name);
  if (it != h1.end()) {
    h1[name].Fill(value);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is not a TH1F." << std::endl;
}

void HistoContainer::Fill(const char* name, float valuex, float valuey) { 

  std::map<const char*, TProfile>::const_iterator itp = hp.find(name);
  if (itp != hp.end()) {
    hp[name].Fill(valuex, valuey);
    return;
  }

  std::map<const char*, TH2F>::const_iterator it2 = h2.find(name);
  if (it2 != h2.end()) {
    h2[name].Fill(valuex, valuey);
    return;
  }

  std::cerr << "ERROR !: histogram " << name << " is nor a TH2F nor a TProfile." << std::endl;
}

void HistoContainer::Save() {
  std::map<const char*, TH1F>::const_iterator it;
  for (it = h1.begin(); it != h1.end(); ++it)
    it->second.Write();

  std::map<const char*, TH2F>::const_iterator it2;
  for (it2 = h2.begin(); it2 != h2.end(); ++it2)
    it2->second.Write(); 

  std::map<const char*, TProfile>::const_iterator itp;
  for (itp = hp.begin(); itp != hp.end(); ++itp)
    itp->second.Write();
}
