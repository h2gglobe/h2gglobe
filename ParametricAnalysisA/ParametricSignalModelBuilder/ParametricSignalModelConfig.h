#ifndef _ParametricSignalModelConfig_h
#define _ParametricSignalModelConfig_h

#include <string>
#include <map>
#include <vector>

/** class for reading / keeping the configuration given by the user */
class ParametricSignalModelConfig
{
public:
  /** name of the file with the signal MC events (and the background model) */
  std::string inputFname;

  /** true if this is to be used for a fermiophobic interpretation of the data */
  bool fermiophobic;

  /** overall number of categories. Should be set to the number of categories
   available in the input file, NOT the categories used for calculating the limits */
  unsigned numCategories;

  /** number of inclusive categories (used to determine which categories are 'special',
   *  i.e. VBF like)
   */
  unsigned numInclusiveCategories;

  /** this can be set to false in case the input file
   *  does NOT distinguish between right and wrong
   *  vertex. Normally left true.
   */
  bool useRightWrongVertex;

  /** can be used to merge signal process input datasets into a combined
   *  dataset for fitting (e.g. keep ggh as it is while merging the rest).
   *
   *  Maps from the suffix used for the RooDataSet (from the background
   *  workspace) to the name used for the merged dataset.
   */
  std::map<std::string, std::string> signalProcessMergingMapping;

  /** maximum mass (used e.g. for interpolating the fitted parameters) */
  double massmax;

  /** distance (in GeV) between masses used for interpolation */
  double massInterpolationWidth;

  /** weight used for scaling when fitting Gaussians */
  double weightscale;

  /** name of the workspace in the input file */
  std::string inputWorkspaceName;

  /** constants for smearing the mass resolution of
   *  MC (to adapt to the data)
   */
  std::vector<double> smearingv;

  /** center of mass energy (used for cross sections) */
  int sqrts;

  /** a suffix to be used in the names of the generated
      objects which can e.g. be used to distinguish
      the datasets of different CM energies
      in order to run them in combination */
  std::string nameSuffix;

  //----------------------------------------------------------------------

  /** constructor to initialize some default values */
  ParametricSignalModelConfig();

  /** reads the configuration from an XML file */
  static ParametricSignalModelConfig read(const std::string &configFname);

private:
  void initializeSignalProcessMergingMapping(bool fermiophobic);

};
#endif
