#include "HistoContainer.h"
#include <utility>
#include <iostream>

#define HCDEBUG 1


HistoContainer::HistoContainer() {}

HistoContainer::HistoContainer(int setTo) {
  setHistVal(setTo);
}

HistoContainer::~HistoContainer() 
{}

void HistoContainer::setHistVal(int setTo) {
  histVal = setTo;
}

int HistoContainer::getHistVal() {
  return histVal;
}


std::string HistoContainer::ModifiedName(char* name, int i) {
  char* modName= new char[500];
  sprintf(modName, "%s_cat%d_%d", name, i, histVal);
  std::string output(modName);
  return output;
}

void HistoContainer::Add(char* name, int categories,int bins, float xmin, float xmax) {

  std::vector<TH1F> temp;

  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i);
    temp.push_back(TH1F(modName.c_str(), modName.c_str(), bins, xmin, xmax));
  }
  
  h1[std::string(name)] = temp;
}


void HistoContainer::Add(char* name, int categories, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
  
  std::vector<TH2F> temp;
  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i);
    temp.push_back(TH2F(modName.c_str(), modName.c_str(), binsx, xmin, xmax, binsy, ymin, ymax));
  }
  
  h2[std::string(name)] = temp;
}

void HistoContainer::Add(char* name, int categories, int binsx, 
			 float xmin, float xmax,
			 float ymin, float ymax) {

  std::vector<TProfile> temp;
  for (int i=0; i<categories; i++) {
    std::string modName = ModifiedName(name, i);
    temp.push_back(TProfile(modName.c_str(), modName.c_str(), binsx, xmin, xmax, ymin, ymax));
  }

  hp[std::string(name)] = temp;
} 

void HistoContainer::Fill(std::string name, int category, float value) {
  Fill(name, category, value, 1.0);
}

void HistoContainer::Fill(std::string name, int category, float value, float weight) {

  //std::string modName = ModifiedName(name);
  std::map<std::string, std::vector<TH1F> >::iterator it = h1.find(std::string(name));
 
  if (it != h1.end()) {
    (it->second)[category].Fill(value, weight);
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
    (itp->second)[category].Fill(valuex, valuey, weight);
    return;
  }

  std::map<std::string, std::vector<TH2F> >::iterator it2 = h2.find(std::string(name));
  if (it2 != h2.end()) {
    (it2->second)[category].Fill(valuex, valuey, weight);
    return;
  }
}

void HistoContainer::Save() {
  std::map<std::string, std::vector<TH1F> >::iterator it;
  for (it = h1.begin(); it != h1.end(); ++it) {
    for (unsigned int i=0; i<(it->second).size(); i++) {
      (it->second)[i].Scale(total_scale);
      (it->second)[i].Write();
    }
  }

  std::map<std::string, std::vector<TH2F> >::iterator it2;
  for (it2 = h2.begin(); it2 != h2.end(); ++it2) {
    for (unsigned int i=0; i<(it2->second).size(); i++) {
      (it2->second)[i].Scale(total_scale);
      (it2->second)[i].Write(); 
    }
  }

  std::map<std::string, std::vector<TProfile> >::iterator itp;
  for (itp = hp.begin(); itp != hp.end(); ++itp) {
    for (unsigned int i=0; i<(itp->second).size(); i++) {
      (itp->second)[i].Scale(total_scale);
      (itp->second)[i].Write();
    }
  }
}
