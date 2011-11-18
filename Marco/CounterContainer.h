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

  std::vector<float> operator[](unsigned int);
  void Add(std::string, int, std::string, std::string, std::string);
  void Fill(std::string, int);
  void Fill(std::string, int, float);
  void Save();
  unsigned int size() { return c.size(); }
  unsigned int ncat(unsigned int);
  std::string name(unsigned int);
  std::string denomName(unsigned int, unsigned int);
  float efficiency(unsigned int, unsigned int, unsigned int);
  float efficiency(unsigned int, unsigned int);
  float tot(unsigned int);

 private:
  int histVal;
  std::vector<std::vector<std::string> > denoms_;
  std::vector<std::string> names;
  std::vector<std::vector<float> > c;
};

#endif
