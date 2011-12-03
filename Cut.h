#ifndef CUT
#define CUT

#include <vector>
#include <string>

class Cut {
 public:
  Cut();
  ~Cut();
  
  void Print();

  std::string name;
  int fromright;
  int finalcut;
  int ncat;
  int useit;
  int index;
  std::vector<float> cut;
  std::vector<float> cutintervall;
  std::vector<float> cutintervalh;
  float* mycutvar;

};
#endif
