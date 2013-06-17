#ifndef TREECONTAINER
#define TREECONTAINER

#include <TFile.h>
#include <TTree.h>
#include <map>
#include <string>

class TreeContainer {

 public:
  TreeContainer();
  TreeContainer(int,std::string, std::string);
  ~TreeContainer();

  void FillFloat(std::string,  float);
  void FillDouble(std::string, double);
  void FillUInt(std::string,  unsigned int);
  void FillInt(std::string,  int);
  void FillString(std::string, std::string);
  void FillBool(std::string, bool);

  void AddTreeBranch(std::string, int);
  template<class T> void AddExternalBranch(const char * name, T* addr) { tr_->Branch(name,addr); };
  template <class T> void AddExternalBranch(const char * name, T* addr, const char*  type) { tr_->Branch(name,addr,type); };
  template <class T> void AddExternalBranch(const char * name, T* addr, int bufsize, int splitlevel) { tr_->Branch(name,addr,bufsize,splitlevel); };
  void FillTree();
  
  int ncat(int n);
  void Save(TFile*);

  void setDirName(std::string);
  void setTreeVal(int);
  int getTreeVal();
  void setTreeNam(std::string);
  std::string ModifiedName(char*, int);
  float total_scale;

 private:
  int treeVal;
  std::string treeNam;
  std::string dirName;
   
  TTree *tr_;
  std::map<std::string, double>         double_branches;
  std::map<std::string, float>          float_branches;
  std::map<std::string, int>            int_branches;
  std::map<std::string, unsigned int>   uint_branches;
  std::map<std::string, std::string>    string_branches;
  std::map<std::string, bool>           bool_branches;

  void resetDefaults();

};

#endif
