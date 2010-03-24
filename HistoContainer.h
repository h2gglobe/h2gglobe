#ifndef HISTOCONTAINER
#define HISTOCONTAINER

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <map>

class HistoContainer {

 public:
  HistoContainer();
  ~HistoContainer();
    
  void Add(const char*, int, float, float);
  void Add(const char*, int, float, float, int, float, float);
  void Add(const char*, int, float, float, float, float);

  void Fill(const char*, float);
  void Fill(const char*, float, float);
  
  void Save();
  
 private:
  std::map<const char*, TH1F> h1;
  std::map<const char*, TH2F> h2;
  std::map<const char*, TProfile> hp;
};

#endif
