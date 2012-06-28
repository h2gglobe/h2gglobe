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
  void Add(const std::string &name, const std::string &xaxisTitle,const std::string &yaxisTitle, 
           int categories, int bins, float xmin, float xmax);

  /** books a set of 2D histograms for a plot */
  void Add(const std::string &name, const std::string &xaxisTitle, const std::string &yaxisTitle, 
           int categories, 
           int xbins, float xmin, float xmax, 
           int ybins, float ymin, float ymax);


  void Add(char *, char*, char*, int, int, float, float, float, float);

  /** fills a histogram given by the name and category */
  void Fill(const std::string &name, int category, float value, float weight = 1.0);
  
  /** fills a 2D histogram or a profile histogram given by name and category */
  void Fill2D(const std::string &name, int category, float valuex, float valuey, float weight = 1.0);

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
  std::string ModifiedName(const char* name, int i);

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
