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

  /** fills a histogram given by the name and category

      @param itype is typically the process number (data, background MC, signal MC)
  */
  void Fill(const std::string &name, int itype, int category, float value, float weight = 1.0);
  
  /** fills a 2D histogram or a profile histogram given by name and category

      @param itype is typically the process number (data, background MC, signal MC)
  */
  void Fill2D(const std::string &name, int itype, int category, float valuex, float valuey, float weight = 1.0);

  /** @return -1 if the plot name is not found */
  int getNumCategories(const std::string &catname);

  /** number of different plots ('names') booked */
  int size() { return (int)names.size(); }

  // int nbins(int, bool);
  // std::string axisName(int, bool);
  // float max(int, bool);
  // float min(int, bool);
  std::string getName(int n) { return names[n]; }
  
  /** @return a list of the names of all booked plots */
  std::vector<std::string> getAllNames() { return names; }

  /** write histograms to the output file */
  void Save();

  // int getDimension(int);
  int getHistVal();
  void setHistVal(int);

  /** sets the global suffix for the histogram names */
  void setHistNam(std::string);

  void setScale(float);

  /** builds the full histogram name from the given name,
      the category etc. */
  std::string ModifiedName(const char* name, int itype, int cat);

  /** overall scaling factor (with which
      the given weights are multiplied) */
  float total_scale;

 private:
  /** returns an 1D histogram from the list of histograms
      or creates it if the corresponding process type
      was never seen before */
  TH1F *Get1Dhisto(const std::string &name, int itype, int category);

  TH2F *Get2Dhisto(const std::string &name, int itype, int category);

  TProfile *GetProfileHisto(const std::string &name, int itype, int category);

  //----------
  /** class with histogram definition parameters (except the name) */
  struct HistoDef
  {
    std::string xtitle;
    std::string ytitle;

    unsigned xbins;

    float xmin, xmax;

    unsigned ybins;
    float ymin, ymax;
    
  };

  /** we must keep these here in order to book the histograms
      per process on demand */
  std::map<std::string, HistoDef> histoDefs;
  //----------

  int histVal;
  std::string histNam;
  std::vector<std::string> names;
  
  /** list of 1D histograms. 
      First index is the 'plot name',
      the second index is the category number 
      the third index is the process number
  */
  std::map<std::string, std::vector< std::map<int, TH1F *> > > h1;

  /** similar to h1 but for 2d histograms */
  std::map<std::string, std::vector< std::map<int, TH2F *> > > h2;

  /** similar to h1 but for profile histograms */
  std::map<std::string, std::vector< std::map<int, TProfile *> > > hp;
};

#endif
