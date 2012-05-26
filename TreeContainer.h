#ifndef TREECONTAINER
#define TREECONTAINER

#include <TTree.h>
#include <map>
#include <string>

class TreeContainer {

 public:
  TreeContainer();
  TreeContainer(int,std::string);
  ~TreeContainer();

  void FillFloat(std::string,  float);
  void FillDouble(std::string, double);
  void FillInt(std::string,  int);

  void AddTreeBranch(std::string, int);
  void FillTree();

  int ncat(int n);
  void Save();

  void setTreeVal(int);
  int getTreeVal();
  void setTreeNam(std::string);
  std::string ModifiedName(char*, int);

 private:
  int treeVal;
  std::string treeNam;
   
  TTree *tr_;
  std::map<std::string, double> double_branches;
  std::map<std::string, float>  float_branches;
  std::map<std::string, int>    int_branches;

  void resetDefaults();

};

#endif
