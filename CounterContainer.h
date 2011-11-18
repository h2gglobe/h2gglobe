#ifndef COUNTERCONTAINER
#define COUNTERCONTAINER

#include <map>
#include <string>
#include <vector>

class CounterContainer {

 public:
  CounterContainer();
  CounterContainer(int);
  ~CounterContainer();

  std::vector<int> operator[](unsigned int);
  void Add(std::string, int, std::string, std::string, std::string);
  void Fill(std::string, int);
  void Fill(std::string, int, float);
  void Save();
  unsigned int mapSize() { return c.size(); }
  unsigned int ncat(int);
  std::string name(int);

 private:
  int histVal;
  std::vector<int> denom1_, denom2_, denom3_;
  std::map<std::string, std::vector<int> > c;
};

#endif
