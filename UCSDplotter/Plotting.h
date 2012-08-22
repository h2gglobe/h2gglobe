#ifndef _Plotting_h
#define _Plotting_h

#include <string>

#include "HistoContainer.h"
#include "RooContainer.h"

class Plotting
{

protected:
  HistoContainer *histoContainer;



public:
  Plotting(HistoContainer *histoContainer);

  /** add an itype to define what signal is.
   *  For the moment does NOT check whether
   *  a signal itype is added more than once. */
  void addSignalItype(int itype);

  /** draw all plots */
  void plotAll();

  /** plot all categories of a given plot name */
  void plot(const std::string &plotName);

  /** plot ONE category of a given plot name */
  void plot(const std::string &plotName, int cat);

protected:
  /** will check whether itype is in signalItypes.
   *  This is used to only plot one signal type,
   *  e.g. if we use different mass hypotheses in the same file.
   */
  bool isSignal(int itype);

  bool isBackground(int itype) { return itype > 0; }

  inline bool isData(int itype) { return itype == 0; }

  /** which itype are signal */
  std::vector<int> signalItypes;

};

#endif
