#ifndef HISTOCONTAINER
#define HISTOCONTAINER

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <map>
#include <string>

class HistoContainer {

 public:
  HistoContainer();
  HistoContainer(int);
  ~HistoContainer();
    
  void Add(char*, int, float, float);
  void Add(char*, int, float, float, int, float, float);
  void Add(char*, int, float, float, float, float);

  void Fill(char*, float);
  void Fill(char*, float, float);
  
  std::map<std::string, TH1F*>* getMap() { return h1; };
  
  void Save();
 
  int getHistVal();
  void setHistVal(int);
  std::string ModifiedName(char*);


 private:
  int histVal;
  std::map<std::string, TH1F*>* h1;
  std::map<std::string, TH2F> h2;
  std::map<std::string, TProfile> hp;
};

#endif
