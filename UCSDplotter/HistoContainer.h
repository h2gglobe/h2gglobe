#ifndef HISTOCONTAINER
#define HISTOCONTAINER

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <map>
#include <string>

/** class keeping track of the histograms needed for plots
 *  (split by category and by processes).
 *
 *
 *  (technical note: we should use TH1 where possible because TH1F/TH1D as well as
 *   TH2F/TH2D and TProfile inherit from it)
 */
class HistoContainer {

 public:
  HistoContainer();
  HistoContainer(int,std::string);
  ~HistoContainer();
    
  enum HistoType
  {
    HISTO_1D = 0,
    HISTO_2D = 1,
    HISTO_PROFILE = 2,
  };

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

  /** get all itypes for a given name and category */
  std::vector<int> getItypes(const std::string &plotName, int cat);

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

  /** produces a sum of the histograms found for the given name and category
   *  and given itypes.
   *
   *  @param sumName should be something distinguishing this histogram
   *   from the existing histograms with the same name and category.
   *   Typical examples are "data", "bg", "sig".
   *   */
  TH1* makeHistoSum(const std::string &name, const std::string &sumName, const std::vector<int> &itypes, int cat);

  /** @return true iff this category was encountered for this plot */
  bool hasCategory(const std::string &plotName, int cat);

  /** overall scaling factor (with which
      the given weights are multiplied) */
  float total_scale;

 private:
  //----------

  /** class with histogram definition parameters (except the name) */
    struct HistoDef
    {
      HistoType type;
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

  /** helper for the other Add(..) methods */
  void Add(const std::string &name, const HistoDef &histoDef, int categories);

  std::map<int, TH1 *> *getMapForNameAndCategory(const std::string &name, int category);

  /** returns an 1D histogram from the list of histograms
      or creates it if the corresponding process type
      was never seen before */
  TH1F *Get1Dhisto(const std::string &name, int itype, int category);

  TH2F *Get2Dhisto(const std::string &name, int itype, int category);

  TProfile *GetProfileHisto(const std::string &name, int itype, int category);


  int histVal;
  std::string histNam;
  std::vector<std::string> names;
  
  /** list of all 1D, 2D and profile histograms.
      First index is the 'plot name',
      the second index is the category number 
      the third index is the process number
  */
  std::map<std::string, std::vector< std::map<int, TH1 *> > > histos;

//  /** similar to h1 but for 2d histograms */
//  std::map<std::string, std::vector< std::map<int, TH2F *> > > h2;
//
//  /** similar to h1 but for profile histograms */
//  std::map<std::string, std::vector< std::map<int, TProfile *> > > hp;

  /** builds the full histogram name from the given name,
      the category etc. */
  std::string fullName(const std::string &name, int itype, int cat);

  /** similar to fullName but uses 'suffix' instead of itype */
  std::string makeSumName(const std::string &name, const std::string &suffix, int cat);

  /** creates a histogram with the specified binning and axis labels etc.
   *  without adding it to any container.
   */
  TH1 *createHistogram(const std::string &name, const HistoDef &histoDef);

};

#endif
