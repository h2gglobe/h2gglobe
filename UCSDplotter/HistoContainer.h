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
  HistoContainer(int,std::string);
  ~HistoContainer();
    

  /** books a set of 1D histograms for a plot */
  void Add(char *name, char *xaxisTitle, char *yaxisTitle, int categories, int bins, float xmin, float xmax);

  
  void Add(char *, char*, char*, int, int, float, float, int, float, float);
  void Add(char *, char*, char*, int, int, float, float, float, float);

  void Fill(std::string, int, float);
  void Fill(std::string, int, float, float);
  
  void Fill2D(std::string, int, float, float);
  void Fill2D(std::string, int, float, float, float);

  int ncat(int n);

  /** number of different plots ('names') booked */
  int size() { return (int)names.size(); }

  int nbins(int, bool);
  std::string axisName(int, bool);
  float max(int, bool);
  float min(int, bool);
  std::string getName(int n) { return names[n]; }
  
  /** write histograms to the output file */
  void Save();

  int getDimension(int);
  int getHistVal();
  void setHistVal(int);
  void setHistNam(std::string);
  void setScale(float);
  std::string ModifiedName(char*, int);

  /** overall scaling factor (with which
      the given weights are multiplied) */
  float total_scale;

 private:
  int histVal;
  std::string histNam;
  std::vector<std::string> names;
  
  /** list of 1D histograms. 
      First index is the 'plot name',
      the second index is the 'category / process' number */
  std::map<std::string, std::vector<TH1F> > h1;

  /** similar to h1 but for 2d histograms */
  std::map<std::string, std::vector<TH2F> > h2;

  /** similar to h1 but for profile histograms */
  std::map<std::string, std::vector<TProfile> > hp;
};

#endif
