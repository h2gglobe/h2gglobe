#include "HistoContainer.h"
#include <utility>
#include <iostream>

#define HCDEBUG 1


HistoContainer::HistoContainer() {
  h1 = new std::map<std::string, TH1F*>();
}

HistoContainer::HistoContainer(int setTo) 
{
  setHistVal(setTo);
  h1 = new std::map<std::string, TH1F*>();
}

HistoContainer::~HistoContainer() 
{}

void HistoContainer::setHistVal(int setTo) {
  histVal = setTo;
}

int HistoContainer::getHistVal() {
  return histVal;
}


std::string HistoContainer::ModifiedName(char* name) {
  char* modName= new char[500];
  sprintf(modName, "%s_%d", name, histVal);
  std::string output(modName);
  return output;
}

void HistoContainer::Add(char* name, int bins, float xmin, float xmax) {
  std::string modName = ModifiedName(name);
  std::pair<std::string, TH1F*> it(modName, new TH1F(modName.c_str(), modName.c_str(), bins, xmin, xmax));

  (*h1)[it.first] = it.second;
}


void HistoContainer::Add( char* name, int binsx, float xmin, float xmax,
			 int binsy, float ymin, float ymax) {
		       
  std::string modName = ModifiedName(name);
  TH2F temp(modName.c_str(), modName.c_str(), binsx, xmin, xmax, binsy, ymin, ymax);
  h2.insert(std::pair< std::string, TH2F>(modName, temp));
}

void HistoContainer::Add( char* name, int binsx, float xmin, float xmax,
			 float ymin, float ymax) {
		       
  std::string modName = ModifiedName(name);
  TProfile temp(modName.c_str(), modName.c_str(), binsx, xmin, xmax, ymin, ymax);
  hp.insert(std::pair< std::string, TProfile>(modName, temp));
} 

void HistoContainer::Fill(char* name, float value) {

  std::string modName = ModifiedName(name);
  std::map<std::string, TH1F*>::iterator it = h1->find(modName);
 
  if (it != h1->end()) {
    (*h1)[modName]->Fill(value);
    return;
  }

  if(HCDEBUG)std::cerr << "ERROR !: histogram " << modName << " is not a TH1F." << std::endl;
}

void HistoContainer::Fill( char* name, float valuex, float valuey) { 

  std::string modName = ModifiedName(name);
  std::map< std::string, TProfile>::const_iterator itp = hp.find(modName);
  if (itp != hp.end()) {
    hp[modName].Fill(valuex, valuey);
    return;
  }

  std::map< std::string, TH2F>::const_iterator it2 = h2.find(modName);
  if (it2 != h2.end()) {
    h2[modName].Fill(valuex, valuey);
    return;
  }

  if(HCDEBUG)std::cerr << "ERROR !: histogram " << modName << " is nor a TH2F nor a TProfile." << std::endl;
}

void HistoContainer::Save() {
  std::map<std::string, TH1F*>::const_iterator it;
  for (it = h1->begin(); it != h1->end(); ++it) {
    std::cout << it->second->GetName() << std::endl;
    it->second->Write();
  }

  std::map< std::string, TH2F>::const_iterator it2;
  for (it2 = h2.begin(); it2 != h2.end(); ++it2)
    it2->second.Write(); 

  std::map< std::string, TProfile>::const_iterator itp;
  for (itp = hp.begin(); itp != hp.end(); ++itp)
    itp->second.Write();
}
