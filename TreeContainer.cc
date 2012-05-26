#include "TreeContainer.h"
#include <utility>
#include <iostream>

#define TCDEBUG 1

TreeContainer::TreeContainer() {}

TreeContainer::TreeContainer(int setTo, std::string setName) {
  setTreeVal(setTo);
  setTreeNam(setName);
}

TreeContainer::~TreeContainer() 
{ 
}

void TreeContainer::setTreeVal(int setTo) {
  treeVal = setTo;
}

void TreeContainer::setTreeNam(std::string setName) {
  treeNam = setName;
  tr_ = new TTree(setName.c_str(),setName.c_str()); tr_->SetDirectory(0);
  if (TCDEBUG) std::cout << "TreeContainer -- Made new Tree Called " << treeNam << std::endl;
}

int TreeContainer::getTreeVal() {
  return treeVal;
}

void TreeContainer::AddTreeBranch(std::string name,int type){

	if (type==0){	// Int_t
		int_branches.insert(std::pair<std::string,int> (name,-999) );	
		tr_->Branch(name.c_str(),&(int_branches[name]),Form("%s/Int_t",name.c_str()));
		if (TCDEBUG)	std::cout << "TreeContainer -- Creating Int Branch " << name << std::endl;
	}
	else if (type==1){	// Float_t
		float_branches.insert(std::pair<std::string,float> (name,-999.) );
		tr_->Branch(name.c_str(),&(float_branches[name]),Form("%s/Float_t",name.c_str()));
		if (TCDEBUG)	std::cout << "TreeContainer -- Creating Float Branch " << name << std::endl;
	}
	else if (type==2){	// Double_t
		double_branches.insert(std::pair<std::string,double> (name,-999.) );	
		tr_->Branch(name.c_str(),&(double_branches[name]),Form("%s/Double_t",name.c_str()));
		if (TCDEBUG)	std::cout << "TreeContainer -- Creating Double Branch " << name << std::endl;
	} else { 
		std::cerr << "TreeContainer -- No Type " << type << std::endl;
	}
}

void TreeContainer::FillFloat(std::string name, float x){

	std::map<std::string,float>::iterator it = float_branches.find(name);
	if (it!=float_branches.end()){
		(*it).second = x;
	} else {
		std::cerr << "TreeContainer -- No Float Tree " << name << std::endl;
	}
}
void TreeContainer::FillInt(std::string name, int x){

	std::map<std::string,int>::iterator it = int_branches.find(name);
	if (it!=int_branches.end()){
		(*it).second = x;
	} else {
		std::cerr << "TreeContainer -- No Double Tree " << name << std::endl;
	}
}

void TreeContainer::FillDouble(std::string name, double x){

	std::map<std::string,double>::iterator it = double_branches.find(name);
	if (it!=double_branches.end()){
		(*it).second = x;
	} else {	
		std::cerr << "TreeContainer -- No Double Tree " << name << std::endl;
	}
}

void TreeContainer::FillTree(){
	tr_->Fill();
	resetDefaults();
}

void TreeContainer::Save(){
	tr_->Write(0,TObject::kWriteDelete);
  	delete tr_; 
}

void TreeContainer::resetDefaults(){
	for (std::map<std::string,int>::iterator it = int_branches.begin();it!=int_branches.end() ;it++){
		(*it).second=-999;
	}
	for (std::map<std::string,float>::iterator it = float_branches.begin();it!=float_branches.end() ;it++){
		(*it).second=-999.;
	}
	for (std::map<std::string,double>::iterator it = double_branches.begin();it!=double_branches.end() ;it++){
		(*it).second=-999.;
	}
}
