#include "TreeContainer.h"
#include <utility>
#include <iostream>

#define TCDEBUG 0

TreeContainer::TreeContainer() {}

TreeContainer::TreeContainer(int setTo, std::string setName, std::string setDirNam) {
  setTreeVal(setTo);
  setTreeNam(setName);
  setDirName(setDirNam);
}

TreeContainer::~TreeContainer() 
{ 
}

void TreeContainer::setTreeVal(int setTo) {
  treeVal = setTo;
}

void TreeContainer::setDirName(std::string setName) {
  dirName = setName;
};

void TreeContainer::setTreeNam(std::string setName) {
  treeNam = setName;
  tr_ = new TTree(setName.c_str(),setName.c_str()); tr_->SetDirectory(0);
  if (TCDEBUG) std::cout << "TreeContainer -- Made new Tree Called " << treeNam << std::endl;
}

int TreeContainer::getTreeVal() {
  return treeVal;
}

void TreeContainer::setScale(float scale){
  total_scale=scale;
}

void TreeContainer::AddTreeBranch(std::string name,int type){

  if (type==0){	// Int_t
    int_branches.insert(std::pair<std::string,int> (name,-999) );	
    tr_->Branch(name.c_str(),&(int_branches[name]),Form("%s/Int_t",name.c_str()));
    if (TCDEBUG)	std::cout << "TreeContainer -- Creating Int Branch " << name << std::endl;
  }
  else if (type==1){	// Unsigned Int_t
    uint_branches.insert(std::pair<std::string,unsigned int> (name,999.) );	
    tr_->Branch(name.c_str(),&(uint_branches[name]),Form("%s/UInt_t",name.c_str()));
    if (TCDEBUG)	std::cout << "TreeContainer -- Creating UInt Branch " << name << std::endl;
  }
  else if (type==2){	// Float_t
    float_branches.insert(std::pair<std::string,float> (name,-999.) );
    tr_->Branch(name.c_str(),&(float_branches[name]),Form("%s/Float_t",name.c_str()));
    if (TCDEBUG)	std::cout << "TreeContainer -- Creating Float Branch " << name << std::endl;
  }
  else if (type==3){	// Double_t
    double_branches.insert(std::pair<std::string,double> (name,-999.) );	
    tr_->Branch(name.c_str(),&(double_branches[name]),Form("%s/Double_t",name.c_str()));
    if (TCDEBUG)	std::cout << "TreeContainer -- Creating Double Branch " << name << std::endl;
  }
  else if (type==4) {     // std::string
    string_branches.insert(std::pair<std::string, std::string> (name, ""));	
    //tr_->Branch(name.c_str(), &(string_branches[name]),Form("%s/C",name.c_str()));
    tr_->Branch(name.c_str(), "std::string", &(string_branches[name]));
    if (TCDEBUG)	std::cout << "TreeContainer -- Creating String Branch " << name << std::endl;
  }
  else if (type==5){
    bool_branches.insert(std::pair<std::string,bool> (name,false) );
    tr_->Branch(name.c_str(),&(bool_branches[name]),Form("%s/Bool_t",name.c_str()));
    if (TCDEBUG) std::cout << "TreeContainer -- Creating Bool Branch " << name << std::endl;
    
  } else { 
    std::cerr << "TreeContainer -- No Type " << type << std::endl;
  }
}

void TreeContainer::FillFloat(std::string name, float x){

  std::map<std::string,float>::iterator it = float_branches.find(name);
  if (it!=float_branches.end()){
    (*it).second = x;
  } else {
    std::cerr << "TreeContainer -- No Float Tree trying Double" << name << std::endl;
    FillDouble(name,x);
  }
}

void TreeContainer::FillInt(std::string name, int x){
  
  std::map<std::string,int>::iterator it = int_branches.find(name);
  if (it!=int_branches.end()){
    (*it).second = x;
  } else {
    std::cerr << "TreeContainer -- No Int Tree " << name << std::endl;
  }
}

void TreeContainer::FillUInt(std::string name, unsigned int x){
  
  std::map<std::string,unsigned int>::iterator it = uint_branches.find(name);
  if (it!=uint_branches.end()){
    (*it).second = x;
  } else {
    std::cerr << "TreeContainer -- No UInt Tree " << name << std::endl;
  }
}

void TreeContainer::FillDouble(std::string name, double x){
  
  std::map<std::string,double>::iterator it = double_branches.find(name);
  if (it!=double_branches.end()){
    (*it).second = x;
  } else {	
    std::cerr << "TreeContainer -- No Double Tree trying Float " << name << std::endl;
    FillFloat(name,x);
  }
}

void TreeContainer::FillString(std::string name, std::string x) {
  std::map<std::string,std::string>::iterator it = string_branches.find(name);
  if (it!=string_branches.end()){
    (*it).second = x;
  } else {	
    std::cerr << "TreeContainer -- No Double Tree " << name << std::endl;
  }
}

void TreeContainer::FillBool(std::string name, bool x){
  std::map<std::string,bool>::iterator it = bool_branches.find(name);
  if (it!=bool_branches.end()){
    (*it).second = x;
  } else {
    std::cerr << "TreeeContainer -- No Bool Tree " << name << std::endl;
  }
}

void TreeContainer::FillTree(){
  tr_->Fill();
  resetDefaults();
}

void TreeContainer::Save(TFile* f) {
  if (dirName != "") {
    
    TDirectory* dir = (TDirectory*)f->Get(dirName.c_str());
    if (!dir) dir = f->mkdir(dirName.c_str(), dirName.c_str());
    dir->cd();
    tr_->Write(0,TObject::kWriteDelete);
    delete tr_;
    f->cd();
  } else {
    tr_->Write(0,TObject::kWriteDelete);
    delete tr_;
  }
}

void TreeContainer::resetDefaults(){
  for (std::map<std::string,int>::iterator it = int_branches.begin();it!=int_branches.end() ;it++){
    (*it).second=-999;
  }
  for (std::map<std::string,unsigned int>::iterator it = uint_branches.begin();it!=uint_branches.end() ;it++){
    (*it).second=999;
  }
  for (std::map<std::string,float>::iterator it = float_branches.begin();it!=float_branches.end() ;it++){
    (*it).second=-999.;
  }
  for (std::map<std::string,double>::iterator it = double_branches.begin();it!=double_branches.end() ;it++){
    (*it).second=-999.;
  }
  for (std::map<std::string,bool>::iterator it = bool_branches.begin();it!=bool_branches.end() ;it++){
    (*it).second=-999.;
  }
}
