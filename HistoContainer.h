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
    
  void Add(char *, int, int, float, float);
  void Add(char *, int, int, float, float, int, float, float);
  void Add(char *, int, int, float, float, float, float);

  void Fill(std::string, int, float);
  void Fill(std::string, int, float, float);
  
  void Fill2D(std::string, int, float, float);
  void Fill2D(std::string, int, float, float, float);
  
  void Save();
 
  void setHistVal(int);
  int getHistVal();
  std::string ModifiedName(char*, int);
  float total_scale;

 private:
  int histVal;
  std::map<std::string, std::vector<TH1F> > h1;
  std::map<std::string, std::vector<TH2F> > h2;
  std::map<std::string, std::vector<TProfile> > hp;
};

#endif
