#include <sstream>
#include <iostream>
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooHistFunc.h"
#include "RooMoment.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooProdPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TEfficiency.h"
#include "RooConstVar.h"

#ifndef __CINT__
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost::assign;
#endif

// #include "XsectDataOrig.C"


using namespace RooFit;
#include "utils.h"
#include "ParametricSignalModelConfig.h"


//----------------------------------------------------------------------
#include "XsectDataDetailed.C"

//----------------------------------------------------------------------

//pt-reweighing (but only for gf samples, use extra bool flag for now)
TH1D *ptweights = 0;
float
ptweight(float genhpt, bool isgfsample)
{
  if (!isgfsample)
    return 1.0;
  if (genhpt < 0)
    return 1.0;
  return ptweights->GetBinContent(ptweights->FindFixBin(genhpt));
}

//----------------------------------------------------------------------

//efficiency scale factors (barrel/endcap for now, but can be extended to full kinematic binning)
float
effweight(bool iseb1, bool iseb2)
{
  const float ebscale = 0.992;
  const float eescale = 0.951;

  float effw = 1.0;

  if (iseb1)
    effw *= ebscale;
  else
    effw *= eescale;

  if (iseb2)
    effw *= ebscale;
  else
    effw *= eescale;

  return effw;
}

//----------------------------------------------------------------------

//append data from file to RooDataSet adding a column which has the weight by given cross section
//(Note that the weight is not enabled at this stage, only made available for subsequent use)
void
append(RooDataSet &data, TFile *infile, TCut sel, Double_t xsec, Bool_t isgfsample)
{

  TDirectory *hdirfwk = (TDirectory*) infile->FindObjectAny("AnaFwkMod");
  const TH1D *hDAllEvents = (TH1D*) hdirfwk->Get("hDAllEvents");

  TDirectory *hdir = (TDirectory*) infile->FindObjectAny("HGGMod");
  TTree *hdata = (TTree*) hdir->Get("hHggNtuple");

  RooRealVar xsecweight("xsecweight", "xsecweight", 1e6 * xsec / (double) hDAllEvents->GetEntries());
  RooRealVar isgf("isgf", "isgf", (double) isgfsample);

  RooArgSet varlistsmall = *data.get();
  varlistsmall.remove(xsecweight, kFALSE, kTRUE);
  varlistsmall.remove(isgf, kFALSE, kTRUE);

  RooDataSet newdata("newdata", "newdata", varlistsmall, RooFit::Import(*hdata), RooFit::Cut(sel));
  newdata.addColumn(xsecweight);
  newdata.addColumn(isgf);

  data.append(newdata);
}

//----------------------------------------------------------------------

/**
 * Creates a new RooDataSet(..), add the events contained in 'indata' with their weight
 * multiplied by 'weightscale'.
 *
 * @param name is the name of the returned RooDataSet */
RooDataSet *
createWeightedDataset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, Double_t weightscale)
{
  ASSERT(indata != NULL);

  RooDataSet *outdata = new RooDataSet(name, "", RooArgList(*mvar, *wvar), wvar->GetName());

  for (Int_t ient = 0; ient < indata->numEntries(); ++ient)
  {
    // get the argset of the ient'th event/entry
    const RooArgSet *ent = indata->get(ient);

    // set the value of mvar to the value of this entry in indata
    mvar->setVal(static_cast<RooAbsReal*> (ent->find(mvar->GetName()))->getVal());

    outdata->add(*mvar, weightscale * indata->weight());
  } // loop over the entries of the input dataset

  return outdata;
}

//----------------------------------------------------------------------

/** appends a weighted copy of 'indata' to 'outdata'. Similar to createWeightedDataset but APPENDS instead
 *  to an existing dataset
 *
 *  @param wvar is the weight variable
 *  @param mvar is the mass variable
 *
 *  @return outdata
 *  */
RooDataSet *
appendWeightedDataset(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, Double_t weightscale)
{
  ASSERT(indata != NULL);

  //RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient = 0; ient < indata->numEntries(); ++ient)
  {
    const RooArgSet *ent = indata->get(ient);
    mvar->setVal(static_cast<RooAbsReal*> (ent->find(mvar->GetName()))->getVal());
    outdata->add(*mvar, weightscale * indata->weight());
  } // loop over the entries of the input dataset

  return outdata;
}
//----------------------------------------------------------------------

/** used for RooFormulaVar: creates a string of the form @i + ... + @j */
string
makeArgSum(int first, int last)
{
  ASSERT(first <= last);
  string retval = "@" + boost::lexical_cast<string>(first);

  for (int i = first + 1; i <= last; ++i)
    retval += " + @" + boost::lexical_cast<string>(i);

  return retval;
}

//----------------------------------------------------------------------


// static const unsigned FITPARAM_INDEX_EFFACC = 14;
// #define FITPARAM_INDEX_EFFACC 14

/** implements the same functionality as found before in the single (big)
 *  function but allows splitting the code into several methods
 *  for easier maintainability.
 */
class FitHmassGlobe4CatClass
{
private:
  //----------------------------------------------------------------------
  // parameters
  //----------------------------------------------------------------------
  ParametricSignalModelConfig config;

  /** inverse mapping of the above */
  map<string, vector<string> > inverseSignalProcessMergingMapping;

  //----------------------------------------------------------------------

  /** the input workspace */
  RooWorkspace* win;

  std::vector<int> mhs;
  std::vector<TString> catnames;

  /** with recent input files, the signal comes in per Feynman-Diagram datasets */
  vector<string> inputSigProcessNames;

  /** signal process names after merging */
  vector<string> mergedSigProcessNames;

  // variable for the higss mass hypothesis
  RooRealVar *mnom;

  //----------

  /** names of the categories ('catX') for which the complex signal mass parametrisation
   *  should be used.
   */
  std::set<TString> complexset;

  /** names of the categories ('catX') for which the simple signal mass parametrisation
   *  should be used.
   */
  std::set<TString> simpleset;

  //----------------------------------------------------------------------

  map<string, float *> sm_xs_br;

  /** for adding together (and weighting) the signal datasets */
  map<string, RooDataSet *> mergedSignalDataSets;

  //----------------------------------------------------------------------


public:
  //----------------------------------------------------------------------
  FitHmassGlobe4CatClass(const ParametricSignalModelConfig &_config) :
    config(_config)
  {
    // open the input file (containing the output of fitter.py)
    TFile *fdata = new TFile(config.inputFname.c_str());

    if (!fdata->IsOpen())
    {
      cerr << "could not open input root file '" << config.inputFname << "'" << endl;
      exit(1);
    }

    win = (RooWorkspace*) fdata->Get(config.inputWorkspaceName.c_str());

    // if we don't do this, we'll get a segmentation fault at the end
    // of the program (when RooFit cleans up)
    gROOT->cd();

    if (win == NULL)
    {
      cerr << "workspace '" << config.inputWorkspaceName << "' not found in file " << config.inputFname << endl;
      exit(1);
    }

    // perform initializations
    initializeMassPoints();
    initializeProcessNames();

    //----------
    // calculate the inverse mapping of the signal process mappings
    //----------
    for (map<string, string>::const_iterator it = config.signalProcessMergingMapping.begin(); it != config.signalProcessMergingMapping.end(); ++it)
    {
      inverseSignalProcessMergingMapping[it->second].push_back(it->first);
    }

    //----------

    // create the list of process names AFTER merging
    for (map<string, string>::const_iterator it = config.signalProcessMergingMapping.begin(); it != config.signalProcessMergingMapping.end(); ++it)
    {
      if (find(mergedSigProcessNames.begin(), mergedSigProcessNames.end(), it->second) == mergedSigProcessNames.end())
        // element not yet seen
        mergedSigProcessNames.push_back(it->second);
    }
    //--------------------

    if (! config.fermiophobic)
    {
      sm_xs_br["ggh"] = SMCrossSections::xs_br_ggh;
      sm_xs_br["vbf"] = SMCrossSections::xs_br_vbf;
      sm_xs_br["wzh"] = SMCrossSections::xs_br_wzh;
      sm_xs_br["tth"] = SMCrossSections::xs_br_tth;
      sm_xs_br["sum"] = SMCrossSections::xs_br_sum;
    }
    else
    {
      sm_xs_br["vbf"] = FPCrossSections::xs_br_vbf;
      sm_xs_br["wzh"] = FPCrossSections::xs_br_wzh;
      sm_xs_br["sum"] = FPCrossSections::xs_br_sum;
    }
    //--------------------

    // TODO: these two groups could actually be merged
    for (unsigned i = 0; i < config.numInclusiveCategories; ++i)
      catnames.push_back(TString::Format("cat%d", i));

    {
      ASSERT(config.numCategories >= config.numInclusiveCategories);
      for (unsigned i = config.numInclusiveCategories; i < config.numCategories; ++i)
        catnames.push_back(TString::Format("cat%d", i));
    }

    //--------------------

    //--------------------

    // initializeMassResolutionSmearing();
    //--------------------

    mnom = new RooRealVar("MH", "m_{h}", 110.0, 100.0, 200.0, "GeV");
    mnom->setConstant();
    //mnom->setRange("plotrange",110,140);
    //--------------------

    initializeSimpleComplexSet();

  }
  //----------------------------------------------------------------------
private:

  void
  initializeProcessNames()
  {

    if (!config.fermiophobic)
    {
      inputSigProcessNames.push_back("ggh");
      inputSigProcessNames.push_back("tth");
    }

    inputSigProcessNames.push_back("vbf");
    inputSigProcessNames.push_back("wzh");

  }

  //----------------------------------------------------------------------

  void
  initializeMassPoints()
  {
    //define mass points
    if (0)
    {
      // older version
      mhs.push_back(110);
      mhs.push_back(115);
      mhs.push_back(120);
      //mhs.push_back(123);
      //mhs.push_back(125);
      mhs.push_back(130);
      //mhs.push_back(135);
      mhs.push_back(140);
      //mhs.push_back(145);
      mhs.push_back(150);
    }

    // newer version depends on the given parameters
    for (double mass = 110; mass <= 150; mass += config.massInterpolationWidth)
      mhs.push_back((int) (mass + 0.5));

    // testing / debugging
    // mhs.clear();
    // mhs.push_back(130);
  }
  //----------------------------------------------------------------------

  void
  initializeMassResolutionSmearing()
  {
    cerr << "this function should not be called any more, data should be taken from the XML file" << endl;
    exit(1);

    //Jun30 May10ReReco numbers
    //   smearingv.push_back(0.0063);
    //   smearingv.push_back(0.012);
    //   smearingv.push_back(0.023);
    //   smearingv.push_back(0.018);
    //   smearingv.push_back(0.0063);
    //   smearingv.push_back(0.012);
    //   smearingv.push_back(0.023);
    //   smearingv.push_back(0.018);

    //smearingv.push_back(0.0257);

    //Jul2 951/pb
    //   smearingv.push_back(0.0096);
    //   smearingv.push_back(0.012);
    //   smearingv.push_back(0.026);
    //   smearingv.push_back(0.022);
    //   smearingv.push_back(0.0096);
    //   smearingv.push_back(0.012);
    //   smearingv.push_back(0.026);
    //   smearingv.push_back(0.022);
    //
    //   smearingv.push_back(0.0065);
    //   smearingv.push_back(0.010);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.0065);
    //   smearingv.push_back(0.010);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.017);

    if (false)
    {
      // LP11
      // the values below seem NOT to depend on the pT split
      // and do NOT depend on R9 in the endcap
      //
      // note that these numbers seem NOT to be exactly the same as e.g. found
      // in https://twiki.cern.ch/twiki/bin/view/CMS/ReloadLP11Results
      config.smearingv.push_back(0.0071); // cat 0: BB     high pT   high R9
      config.smearingv.push_back(0.010); // cat 1: BB     high pT   low R9
      config.smearingv.push_back(0.017); // cat 2: mixed  high pT   high R9
      config.smearingv.push_back(0.017); // cat 3: mixed  high pT   low R9

      config.smearingv.push_back(0.0071); // cat 4: BB     low pT    high R9
      config.smearingv.push_back(0.010); // cat 5: BB     low pT    low R9
      config.smearingv.push_back(0.017); // cat 6: mixed  low pT    high R9
      config.smearingv.push_back(0.017); // cat 7: mixed  low pT    low R9
    }

    // new test: 2011-11-15 with full dataset: SM only has 4 categories
    config.smearingv.push_back(0.0071); // cat 0: BB     high pT   high R9
    config.smearingv.push_back(0.010); // cat 1: BB     high pT   low R9
    config.smearingv.push_back(0.017); // cat 2: mixed  high pT   high R9
    config.smearingv.push_back(0.017); // cat 3: mixed  high pT   low R9

    // discussed on 2011-10-14 with Marco and Chris:
    // for testing, for the moment, just put 1% smearing
    // for all additional categories
    while (config.smearingv.size() < catnames.size())
      config.smearingv.push_back(0.01);

    //lregression sept6 nomcats
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.019);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.019);

    //regression sept6 altcats
    //   smearingv.push_back(0.008);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.019);
    //   smearingv.push_back(0.008);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.017);
    //   smearingv.push_back(0.019);

    //oldcats
    //   smearingv.push_back(0.008);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.020);
    //   smearingv.push_back(0.018);
    //   smearingv.push_back(0.008);
    //   smearingv.push_back(0.009);
    //   smearingv.push_back(0.020);
    //   smearingv.push_back(0.018);

  }

  //----------------------------------------------------------------------

  void
  initializeSimpleComplexSet()
  {

    complexset.insert("cat0");
    complexset.insert("cat1");
    //complexset.insert("cat2");
    //complexset.insert("cat3");
    complexset.insert("cat4");
    //complexset.insert("cat5");
    //complexset.insert("cat6");
    //complexset.insert("cat4");
    //complexset.insert("cat5");


    //simpleset.insert("cat2");
    //simpleset.insert("cat3");
    //simpleset.insert("cat6");
    //simpleset.insert("cat7");

  }
  //----------------------------------------------------------------------

  /** this uses an internal map of RooDataSets */
  RooDataSet *
  createOrAppendToDataset(RooDataSet *dataToAppend, RooRealVar *massVar, RooRealVar *weightVar, const string &datasetName, Double_t weightscale)
  {
    // check whether the given dataset already exists or not
    map<string, RooDataSet *>::const_iterator it = mergedSignalDataSets.find(datasetName);
    if (it == mergedSignalDataSets.end())
    {
      // dataset not found, create a new one
      RooDataSet *retval = createWeightedDataset(dataToAppend, massVar, weightVar, datasetName, weightscale);
      mergedSignalDataSets[datasetName] = retval;
      return retval;
    }
    else
      // dataset found, append to it
      return appendWeightedDataset(it->second, dataToAppend, massVar, weightVar, weightscale);
  }
  //----------------------------------------------------------------------

  /** get a RooDataSet from mergedSignalDataSets (i.e. the map internal to this class).
   *  If it does not exist, exits the program. */
  RooDataSet *
  getMergedWeightedDataSet(const string &datasetName)
  {
    map<string, RooDataSet *>::const_iterator it = mergedSignalDataSets.find(datasetName);
    if (it == mergedSignalDataSets.end())
    {
      cerr << "dataset '" << datasetName << "' does not exist in internal map" << endl;
      exit(1);
    }

    return it->second;
  }

  //----------------------------------------------------------------------

  /** @return the dataset with the given name from the input (background) workspace.
   *  If such a workspace does not exist, exits the program.
   */
  RooDataSet *
  getDataSet(const string &name)
  {
    RooDataSet *retval = (RooDataSet*) win->data(name.c_str());
    if (retval == NULL)
    {
      cerr << "dataset '" << name << "' does not exist in workspace '" << config.inputWorkspaceName << "' in file " << config.inputFname << endl;
      exit(1);
    }

    return retval;
  }

  RooDataSet *
  getDataSet(const TString &name)
  {
    return getDataSet(string((const char*) name));
  }

  RooDataSet *
  getDataSet(const char *name)
  {
    return getDataSet(string(name));
  }

  //----------------------------------------------------------------------

  /** makes the 'uber'-PDF for the signal model, which includes
   *  all categories and nuisance parameter variables ('knobs for systematics')
   */
  void
  makeFinalPDF(
      // input parameters
      unsigned numFitParms, // actually unused

      // these depend on the merged process name
      map<unsigned, map<string, map<unsigned, RooHistFunc *> > > &fitparmfuncs,

      // these depend on the input process name
      map<unsigned, map<string, RooHistFunc *> > &effaccfuncs,

      RooRealVar *hmass,
      RooRealVar &nuissancedeltafracright,

      // output
      // first index is the category, the second the signal process name
      map<unsigned, RooRealVar *> &nuissanceDeltaSmears, map<unsigned, map<string, RooAddPdf *> > &combHVtxSlides, map<unsigned, map<string, RooAddPdf *> > &combhvtxminslides, map<unsigned, map<string, RooAddPdf *> > &combhvtxsimpleslides,
      map<unsigned, map<string, RooAbsPdf *> > &finalpdfslides, map<unsigned, map<string, RooAbsReal *> > &finalnormslides)
  {
    // index is the category number

    RooConstVar **mcToDataMassSmearingVar = new RooConstVar*[catnames.size()]; // constant !
    // nuissanceDeltaSmears = new RooRealVar*[catnames.size()];      // nuisance parameter
    map<unsigned, RooRealVar *> nuissanceDeltaMs; // nuisance parameter for the delta mu (1,2,3 fully correlated)


    RooFormulaVar **smearmods = new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > mean1slides;
    map<unsigned, map<string, RooFormulaVar *> > mean2slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > mean3slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > sigma1slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > sigma2slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > sigma3slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > wmean1slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > wmean2slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > wsigma1slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > wsigma2slides;//= new RooFormulaVar*[catnames.size()];
    map<unsigned, map<string, RooGaussian *> > g1slides;//= new RooGaussian*[catnames.size()];
    map<unsigned, map<string, RooGaussian *> > g2slides;//= new RooGaussian*[catnames.size()];
    map<unsigned, map<string, RooGaussian *> > g3slides;//= new RooGaussian*[catnames.size()];
    map<unsigned, map<string, RooAddPdf *> > combhslides;//= new RooAddPdf*[catnames.size()];
    map<unsigned, map<string, RooAddPdf *> > combhminslides;//= new RooAddPdf*[catnames.size()];
    map<unsigned, map<string, RooGaussian *> > wg1slides;//= new RooGaussian*[catnames.size()];
    map<unsigned, map<string, RooGaussian *> > wg2slides;//= new RooGaussian*[catnames.size()];
    map<unsigned, map<string, RooAddPdf *> > combhwrongslides;//= new RooAddPdf*[catnames.size()];
    map<unsigned, map<string, RooFormulaVar *> > fracrightmodslides;//= new RooFormulaVar*[catnames.size()];

    // combHVtxSlides       = new RooAddPdf*[catnames.size()];
    // combhvtxminslides    = new RooAddPdf*[catnames.size()];
    // combhvtxsimpleslides = new RooAddPdf*[catnames.size()];

    // finalpdfslides       = new RooAbsPdf*[catnames.size()];

    // finalnormslides      = new RooAbsReal*[catnames.size()];

    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      TString catname = catnames.at(icat);

      // a constant with the same value as in smearingv
      mcToDataMassSmearingVar[icat] = new RooConstVar(TString("smear") + catname, "", config.smearingv.at(icat));

      // these DO NOT depend on the signal process
      if (icat < 4)
      {
        nuissanceDeltaSmears[icat] = new RooRealVar(TString("CMS_hgg_nuissancedeltasmear") + catname, "", 0.0, - config.smearingv.at(icat), config.smearingv.at(icat));
        nuissanceDeltaSmears[icat]->setConstant();
        nuissanceDeltaMs[icat] = new RooRealVar(TString("CMS_hgg_nuissancedeltam") + catname, "", 0.0, -5.0, 5.0);
        nuissanceDeltaMs[icat]->setConstant();
      }
      else if (icat < config.numInclusiveCategories)
      {
        // take the quantities for categories 4..7 the same as for 0..3
        nuissanceDeltaSmears[icat] = nuissanceDeltaSmears[icat - 4];
        nuissanceDeltaMs[icat] = nuissanceDeltaMs[icat - 4];
      }
      else
      {
        // for categories >= 8 (VBF / VHAD in our case)
        nuissanceDeltaSmears[icat] = nuissanceDeltaSmears[0];
        nuissanceDeltaMs[icat] = nuissanceDeltaMs[0];
      }

      //                                                                   mhGen *(mcToDataMassSmearingVar + nuissancedeltasmears)
      smearmods[icat] = new RooFormulaVar(TString("smearmod") + catname, "", "@0 * (@1 + @2)", RooArgList(*mnom, *mcToDataMassSmearingVar[icat], *nuissanceDeltaSmears[icat]));

      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            // BOOST_FOREACH(string sigProcName, inputSigProcessNames)
            //for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              //----------------------------------------
              //                                               name                                        title  formula                                           dependents
              //----------------------------------------

              // we take the fit parameters (apart from the normalization which is NOT fitted)
              // from the merged processes while building everything else
              // per input process name
              string mergedSigProcName = this->config.signalProcessMergingMapping[inputSigProcName];

              // fitted delta means (w.r.t to generated Higgs mass), depending on nuisance parameter                                                          fitted delta mean                 nuisance parameter
              mean1slides[icat][inputSigProcName] = new RooFormulaVar(TString("mean1slide") + catname + "_" + inputSigProcName, "", "@0 + @1 + @0*@2", RooArgList(*mnom, *fitparmfuncs[icat][mergedSigProcName][0], *nuissanceDeltaMs[icat]));
              mean2slides[icat][inputSigProcName] = new RooFormulaVar(TString("mean2slide") + catname + "_" + inputSigProcName, "", "@0 + @1 + @0*@2", RooArgList(*mnom, *fitparmfuncs[icat][mergedSigProcName][1], *nuissanceDeltaMs[icat]));
              mean3slides[icat][inputSigProcName] = new RooFormulaVar(TString("mean3slide") + catname + "_" + inputSigProcName, "", "@0 + @1 + @0*@2", RooArgList(*mnom, *fitparmfuncs[icat][mergedSigProcName][2], *nuissanceDeltaMs[icat]));

              // sqrt[ (fitted sigma)^2 - (mHgen * smears)^2 + smearmods^2 ]
              sigma1slides[icat][inputSigProcName] = new RooFormulaVar(TString("sigma1slide") + catname + "_" + inputSigProcName, "", "TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",
                  RooArgList(*fitparmfuncs[icat][mergedSigProcName][3], *smearmods[icat], *mcToDataMassSmearingVar[icat], *mnom));
              sigma2slides[icat][inputSigProcName] = new RooFormulaVar(TString("sigma2slide") + catname + "_" + inputSigProcName, "", "TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",
                  RooArgList(*fitparmfuncs[icat][mergedSigProcName][4], *smearmods[icat], *mcToDataMassSmearingVar[icat], *mnom));
              sigma3slides[icat][inputSigProcName] = new RooFormulaVar(TString("sigma3slide") + catname + "_" + inputSigProcName, "", "TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",
                  RooArgList(*fitparmfuncs[icat][mergedSigProcName][5], *smearmods[icat], *mcToDataMassSmearingVar[icat], *mnom));

              wmean1slides[icat][inputSigProcName] = new RooFormulaVar(TString("wmean1slide") + catname + "_" + inputSigProcName, "", "@0 + @1 + @0*@2", RooArgList(*mnom, *fitparmfuncs[icat][mergedSigProcName][8], *nuissanceDeltaMs[icat]));
              wmean2slides[icat][inputSigProcName] = new RooFormulaVar(TString("wmean2slide") + catname + "_" + inputSigProcName, "", "@0 + @1 + @0*@2", RooArgList(*mnom, *fitparmfuncs[icat][mergedSigProcName][9], *nuissanceDeltaMs[icat]));
              wsigma1slides[icat][inputSigProcName] = new RooFormulaVar(TString("wsigma1slide") + catname + "_" + inputSigProcName, "", "TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))",
                  RooArgList(*fitparmfuncs[icat][mergedSigProcName][10], *smearmods[icat], *mcToDataMassSmearingVar[icat], *mnom));
              wsigma2slides[icat][inputSigProcName] = new RooFormulaVar(TString("wsigma2slide") + catname, "", "TMath::Max(0.01,sqrt(@0*@0-@3*@3*@2*@2 +@1*@1))", RooArgList(*fitparmfuncs[icat][mergedSigProcName][11], *smearmods[icat], *mcToDataMassSmearingVar[icat], *mnom));

              //                                     name                                              indep var
              //                                                                                    title      mean               sigma
              g1slides[icat][inputSigProcName] = new RooGaussian(TString("g1slide") + catname + "_" + inputSigProcName, "", *hmass, *mean1slides[icat][inputSigProcName], *sigma1slides[icat][inputSigProcName]);
              g2slides[icat][inputSigProcName] = new RooGaussian(TString("g2slide") + catname + "_" + inputSigProcName, "", *hmass, *mean2slides[icat][inputSigProcName], *sigma2slides[icat][inputSigProcName]);
              g3slides[icat][inputSigProcName] = new RooGaussian(TString("g3slide") + catname + "_" + inputSigProcName, "", *hmass, *mean3slides[icat][inputSigProcName], *sigma3slides[icat][inputSigProcName]);

              combhslides[icat][inputSigProcName] = new RooAddPdf(TString("combhslide") + catname + "_" + inputSigProcName, "", RooArgList(*g1slides[icat][inputSigProcName], *g2slides[icat][inputSigProcName], *g3slides[icat][inputSigProcName]),
                  RooArgList(*fitparmfuncs[icat][mergedSigProcName][6], *fitparmfuncs[icat][mergedSigProcName][7]), kTRUE);
              
              // 2 gaussians is the new standard in complex
              combhminslides[icat][inputSigProcName] = new RooAddPdf(TString("combhminslide") + catname + "_" + inputSigProcName, "", RooArgList(*g1slides[icat][inputSigProcName], *g2slides[icat][inputSigProcName]), RooArgList(*fitparmfuncs[icat][mergedSigProcName][6]), kTRUE);

              wg1slides[icat][inputSigProcName] = new RooGaussian(TString("wg1slide") + catname + "_" + inputSigProcName, "", *hmass, *wmean1slides[icat][inputSigProcName], *wsigma1slides[icat][inputSigProcName]);
              wg2slides[icat][inputSigProcName] = new RooGaussian(TString("wg2slide") + catname + "_" + inputSigProcName, "", *hmass, *wmean2slides[icat][inputSigProcName], *wsigma2slides[icat][inputSigProcName]);
              combhwrongslides[icat][inputSigProcName] = new RooAddPdf(TString("combhwrongslide") + catname + "_" + inputSigProcName, "", RooArgList(*wg1slides[icat][inputSigProcName], *wg2slides[icat][inputSigProcName]), RooArgList(*fitparmfuncs[icat][mergedSigProcName][12]));

              fracrightmodslides[icat][inputSigProcName] = new RooFormulaVar(TString("fracrightmodslide") + catname + "_" + inputSigProcName, "", "@0*@1", RooArgList(nuissancedeltafracright, *fitparmfuncs[icat][mergedSigProcName][13]));
              combHVtxSlides[icat][inputSigProcName] = new RooAddPdf(TString("combhvtxslide") + catname + "_" + inputSigProcName, "", RooArgList(*combhminslides[icat][inputSigProcName], *combhwrongslides[icat][inputSigProcName]),
                  RooArgList(*fracrightmodslides[icat][inputSigProcName]), kTRUE);
              
              combhvtxminslides[icat][inputSigProcName] = new RooAddPdf(TString("combhvtxminslide") + catname + "_" + inputSigProcName, "", RooArgList(*combhminslides[icat][inputSigProcName], *wg1slides[icat][inputSigProcName]),
                  RooArgList(*fracrightmodslides[icat][inputSigProcName]), kTRUE);
              combhvtxsimpleslides[icat][inputSigProcName] = new RooAddPdf(TString("combhvtxsimpleslide") + catname + "_" + inputSigProcName, "", RooArgList(*g1slides[icat][inputSigProcName], *wg1slides[icat][inputSigProcName]),
                  RooArgList(*fracrightmodslides[icat][inputSigProcName]), kTRUE);

              //                                   14th fit parameter
              //                                   of category icat
              // normalization ?

              // what is happening here ? overwriting of signal processes ?
              // ASSERT(fitparmfuncs[icat][mergedSigProcName][FITPARAM_INDEX_EFFACC] != NULL);
              // finalnormslides[icat][inputSigProcName] = fitparmfuncs[icat][mergedSigProcName][FITPARAM_INDEX_EFFACC];
              ASSERT(effaccfuncs[icat][inputSigProcName] != NULL);
              finalnormslides[icat][inputSigProcName] = effaccfuncs[icat][inputSigProcName];
              Bool_t usecomplexmodel = complexset.count(catname);
              Bool_t usesimplemodel = simpleset.count(catname);

              if (usecomplexmodel)
                finalpdfslides[icat][inputSigProcName] = combHVtxSlides[icat][inputSigProcName];
              else if (usesimplemodel)
                finalpdfslides[icat][inputSigProcName] = combhvtxsimpleslides[icat][inputSigProcName];
              else
                finalpdfslides[icat][inputSigProcName] = combhvtxminslides[icat][inputSigProcName];

              finalpdfslides[icat][inputSigProcName]->SetName(TString("hggpdf_") + catname + "_" + inputSigProcName);
            } // loop over signal processes
    } // loop over categories
  }

  //----------------------------------------------------------------------

  /** creates a histogram (as function of the mass) with the given name.
   *  This is mostly used to create an interpolating function for
   *  fitted shape parameters afterwards. */
  TH1F *
  makeMassFunctionHist(const std::string &name)
  {
    return new TH1F(name.c_str(), name.c_str(),

    (config.massmax - 100.0 - config.massInterpolationWidth) / config.massInterpolationWidth, // number of bins
        100 + config.massInterpolationWidth / 2.0, // lower edge
        config.massmax - config.massInterpolationWidth / 2.0 // upper edge
    );
  }
  //----------------------------------------------------------------------

  //  inline unsigned getIndex(int cat, int paramIndex, int prodMechIndex)
  //  {
  //
  //  }

  /** number of times each parameter is fit per mass point */
  inline unsigned
  numCases()
  {
    return config.numCategories * mergedSigProcessNames.size();
  }

  //----------------------------------------------------------------------

  /** @return the sum of the cross sections of all relevant production mechanisms
   *  for the given mass.
   */
  double
  getTotalCrossSection(int mass)
  {
    std::stringstream numstringstr;
    numstringstr << mass;
    TString numstring(numstringstr.str());

    double totalCrossSectionValue = 0;

    BOOST_FOREACH(string sigProcName, inputSigProcessNames)
          {
            // cross section times BR variable
            RooRealVar *mtotalxsec = NULL;

            if (!config.fermiophobic)
              // standard model
              mtotalxsec = win->var(TString("XSBR_") + sigProcName + "_" + numstring);
            else
              // fermiophobic cross section times BR
              mtotalxsec = win->var(TString("ff_XSBR_") + sigProcName + "_" + numstring);

            if (mtotalxsec == NULL)
              cerr << "WARNING: cross section does not exist for process " + sigProcName + " and mass " + numstring << endl;
            ASSERT(mtotalxsec != NULL);

            totalCrossSectionValue += mtotalxsec->getVal();

          } // loop over signal Feynman diagrams

    return totalCrossSectionValue;
  }

  //----------------------------------------------------------------------

  RooHistFunc*
  makeCrossSectionFunction(const string &sigProcName)
  {
    // this is actually cross section times branching ratio

    // histogram of cross section as function of the Higgs mass.
    TH1F *crossSectionHistogram;

    // TODO: choose a more appropriate naming for these for the case of Fermiophobic Higgs.
    RooDataHist *dsmxsecs;

    crossSectionHistogram = new TH1F(Form("hsmxsecs_%s", sigProcName.c_str()), "", 81, 109.75, 150.25);

    for (unsigned ipoint = 0; ipoint < SMCrossSections::numPoints; ++ipoint)
    {
      //hsmxsecs->Fill(smmasses[ipoint],smxsecs[ipoint]*smbrs[ipoint]);

      // cout << "FILLING XSECT, fermiophobic=" << fermiophobic << endl;

      if (!config.fermiophobic)
      {
        // STANDARD MODEL
        crossSectionHistogram->Fill(SMCrossSections::masses[ipoint], sm_xs_br[sigProcName][ipoint]);
      }
      else
      {
        // FERMIOPHOBIC
        crossSectionHistogram->Fill(FPCrossSections::masses[ipoint], sm_xs_br[sigProcName][ipoint]);
         // crossSectionHistogram[sigProcName]->Fill(SMCrossSections::smmasses[ipoint], SMCrossSections::ffxsbr[ipoint]);
      }

      // FOUR GENERATIONS STANDARD MODEL
      //hsmxsecs->Fill(smmasses[ipoint],sm4xsbr[ipoint]);

    } // loop over all mass points

    dsmxsecs = new RooDataHist(Form("dsmxsecs_%s", sigProcName.c_str()), "", RooArgList(*mnom), crossSectionHistogram);

    RooHistFunc *funcXsecNorm = new RooHistFunc(Form("fsmxsecs_%s", sigProcName.c_str()), "", RooArgList(*mnom), *dsmxsecs, 1);
    //--------------------
    // plot the cross section times branching ratio
    //--------------------
    {
      TCanvas *csmxsecs = new TCanvas;
      RooPlot *plotsmxsecs = mnom->frame();

      dsmxsecs->plotOn(plotsmxsecs);
      funcXsecNorm->plotOn(plotsmxsecs);

      plotsmxsecs->Draw();
      csmxsecs->SaveAs(Form("xsecs_%s.png", sigProcName.c_str()));
    }

    return funcXsecNorm;
  }

  //----------------------------------------------------------------------

  /** assigns a ratio with ClopperPearson uncertainty to a RooRealVar */
  static void
  assignRatio(RooRealVar &var, double numerator, double denominator)
  {
    double ratio = numerator / denominator;
    double errorLow = TEfficiency::ClopperPearson(Int_t(denominator), Int_t(numerator), 0.683, kFALSE) - ratio;
    double errorHigh = TEfficiency::ClopperPearson(Int_t(denominator), Int_t(numerator), 0.683, kTRUE) - ratio;

    var.setVal(ratio);
    var.setAsymError(errorLow, errorHigh);
  }
  //----------------------------------------------------------------------

  /** creates the objects named hggpdf_cat%d_%s_norm */
  RooFormulaVar*
  makeNsigCat(unsigned cat, const string &inputSigProcName, RooRealVar &intlumi, RooRealVar *totalxsec,
              RooFormulaVar *effaccbarrel,
              RooFormulaVar *effaccmixed,
              RooFormulaVar *r9fracbarrel,
              RooFormulaVar *r9fracmixed,

              // pt fraction variables
              RooAbsReal *ptfracbarrelhighr9,
              RooAbsReal *ptfracbarrelmixedr9,
              RooAbsReal *ptfracmixedhighr9,
              RooAbsReal *ptfracmixedmixedr9,

              RooAbsReal *ptfrac_B_barrelhighr9,
              RooAbsReal *ptfrac_B_barrelmixedr9,
              RooAbsReal *ptfrac_B_mixedhighr9,
              RooAbsReal *ptfrac_B_mixedmixedr9,

              map<unsigned, map<string, map<unsigned, RooHistFunc*> > > &fitparmfuncsX,
              map<unsigned, map<string, RooHistFunc*> > &effaccfuncs
  )
  {
    // need to programmatically generate the list of all processes to normalize to

    // 2012-01-14: no reason why the following code should not work for the fermiophobic case ?
    // assert(! config.fermiophobic);

    const vector<string> *signalProcessNames = NULL;
    // when we take the normalization per signal process but
    // the shapes from merged categories
    signalProcessNames = &inputSigProcessNames;

    // when taking the fits and normalizations from the same signal production mechanism groups
    // signalProcessNames = &mergedSigProcessNames;

    RooArgList args(intlumi, *totalxsec);

    string varExp;

    string ptFractionExpr;

    const unsigned numNumeratorArgs = 7;
    // e.g. (@6 + @7 + @8 + @9)
    string denomExp = "(" + makeArgSum(numNumeratorArgs, numNumeratorArgs + signalProcessNames->size() - 1) + ")";

    switch (cat)
    {
    case 0: // BB / high pT / high R9
      //       lumi      effacc         ptfrac(r9cat)
      //             xsect     r9frac
      varExp = "@0 * @1 * @2 * (@3)    *(@4)  * (@6) / " + denomExp;
      args.add(*effaccbarrel);       // @2
      args.add(*r9fracbarrel);       // @3
      args.add(*ptfracbarrelhighr9); // @4
      args.add(*ptfrac_B_barrelhighr9); // @5
      break;

    case 1: // BB     high pT   low R9
      varExp = "@0 * @1 * @2 * (1.0-@3)*(@4)  * (@6) / " + denomExp;
      args.add(*effaccbarrel);
      args.add(*r9fracbarrel);
      args.add(*ptfracbarrelmixedr9);
      args.add(*ptfrac_B_barrelmixedr9);
      break;

    case 2: // mixed  high pT   high R9
      varExp = "@0 * @1 * @2 * (@3)    *(@4)  * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedhighr9);
      args.add(*ptfrac_B_mixedhighr9);
      break;

    case 3: // mixed  high pT   low R9
      varExp = "@0 * @1 * @2 * (1.0-@3)*(@4)  * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedmixedr9);
      args.add(*ptfrac_B_mixedmixedr9);
      break;
 
      //----------
      // low pT categories
      //----------

    case 4: // BB     low pT    high R9
      varExp = "@0*@1*@2*(@3)* (1.0 - @4 - @5) * (@6) / " + denomExp;
      args.add(*effaccbarrel);
      args.add(*r9fracbarrel);
      args.add(*ptfracbarrelhighr9);
      args.add(*ptfrac_B_barrelhighr9);
      break;

    case 5: // BB     low pT    low R9
      varExp = "@0*@1*@2*(1.0-@3)*(1.0-@4 - @5) * (@6) / " + denomExp;
      args.add(*effaccbarrel);
      args.add(*r9fracbarrel);
      args.add(*ptfracbarrelmixedr9);
      args.add(*ptfrac_B_barrelmixedr9);
      break;

    case 6: // mixed  low pT    high R9
      varExp = "@0*@1*@2*(@3)*(1.0-@4 - @5) * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedhighr9);
      args.add(*ptfrac_B_mixedhighr9);
      break;

    case 7: // mixed  low pT    low R9 
      varExp = "@0*@1*@2*(1.0-@3)*(1.0-@4 - @5) * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedmixedr9);
      args.add(*ptfrac_B_mixedmixedr9);
      break;

    //----------
    // highest pT categories
    //----------

    case 8: // BB     highest pT    high R9
      varExp = "@0*@1*@2*(@3)* (@5) * (@6) / " + denomExp;
      args.add(*effaccbarrel);
      args.add(*r9fracbarrel);
      args.add(*ptfracbarrelhighr9);
      args.add(*ptfrac_B_barrelhighr9);
      break;

    case 9: // BB     highest pT    low R9
      varExp = "@0*@1*@2*(1.0-@3)*(@5) * (@6) / " + denomExp;
      args.add(*effaccbarrel);
      args.add(*r9fracbarrel);
      args.add(*ptfracbarrelmixedr9);
      args.add(*ptfrac_B_barrelmixedr9);
      break;

    case 10: // mixed  highest pT    high R9
      varExp = "@0*@1*@2*(@3)*(@5) * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedhighr9);
      args.add(*ptfrac_B_mixedhighr9);
      break;

    case 11: // mixed  highest pT    low R9
      varExp = "@0*@1*@2*(1.0-@3)*(@5) * (@6) / " + denomExp;
      args.add(*effaccmixed);
      args.add(*r9fracmixed);
      args.add(*ptfracmixedmixedr9);
      args.add(*ptfrac_B_mixedmixedr9);
      break;

    default:
      ASSERT(false)
      ;
      }

    // fraction of this process' efficiency (@6)
    ASSERT(effaccfuncs[cat][inputSigProcName] != NULL);
    args.add(*effaccfuncs[cat][inputSigProcName]);

    // consistency check
    ASSERT((unsigned)args.getSize() == numNumeratorArgs);

    // These are used in the denominator
    // note: must loop o
    BOOST_FOREACH(string proc, (*signalProcessNames))
          {
            // cout << "YYYYYYY ADDING " << proc << endl;
            RooHistFunc* arg = effaccfuncs[cat][proc];
            ASSERT(arg != NULL);
            args.add(*arg);
          }

    // cout << "YYYYYYY ARGS=" << args.getSize() << endl;

    //      args.add(*fitparmfuncsX[cat]["ggh"][FITPARAM_INDEX_EFFACC]);
    //      args.add(*fitparmfuncsX[cat]["vbf"][FITPARAM_INDEX_EFFACC]);
    //      args.add(*fitparmfuncsX[cat]["wzh"][FITPARAM_INDEX_EFFACC]);
    //      args.add(*fitparmfuncsX[cat]["tth"][FITPARAM_INDEX_EFFACC]);

    return new RooFormulaVar(Form("hggpdf_cat%d_%s_norm", cat, inputSigProcName.c_str()), // name
        "", // title
        varExp.c_str(), // formula
        args);
  }

  //----------------------------------------------------------------------

  void
  makeWeightedDataSets(RooRealVar *weight, RooRealVar *hmass, RooAbsPdf *allpdf)
  {
    RooDataSet *mcsigdata, *mcsigwrongdata, // before weighting
        *mcsigwdata, *mcsigwrongwdata; // weighted

    // loop over categories for fits
    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      string catname = (const char *) catnames.at(icat);

      // loop over input process names but potentially merge to output process names
      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              // apply the mapping between input signal process names and output
              // signal process names
              string outputSigProcName = config.signalProcessMergingMapping[inputSigProcName];

              // loop over mass points for fits
              for (UInt_t i = 0; i < mhs.size(); ++i)
              {
                std::stringstream numstringstr;
                numstringstr << mhs.at(i);
                string numstring(numstringstr.str());

                if (config.useRightWrongVertex)
                {
                  //----------
                  //right vertex
                  string name = "sig_" + inputSigProcName + "_mass_m" + numstring + "_rv_" + catname;
                  mcsigdata = getDataSet(name);

                  mcsigwdata = createOrAppendToDataset(mcsigdata, hmass, weight, "mcsigwdata" + numstring + catname + "_" + outputSigProcName, config.weightscale);

                  //----------
                  //wrong vertex
                  name = "sig_" + inputSigProcName + "_mass_m" + numstring + "_wv_" + catname;
                  mcsigwrongdata = getDataSet(name);

                  mcsigwrongwdata = createOrAppendToDataset(mcsigwrongdata, hmass, weight, "mcsigwrongwdata" + numstring + catname + "_" + outputSigProcName, config.weightscale);

                }
                else
                {
                  // for the sake of simplicity, just assume ALL signal events
                  // have the correctly assigned vertex

                  string name = "sig_" + inputSigProcName + "_mass_m" + numstring + "_" + catname;
                  mcsigdata = getDataSet(name);

                  mcsigwdata = createOrAppendToDataset(mcsigdata, hmass, weight, "mcsigwdata" + numstring + catname + "_" + outputSigProcName, config.weightscale);
                } // if no distinction between right and wrong vertex

                //----------
                //combined (right plus wrong vertex)

                // for IC inputs: loop over different signal Feynman diagrams
                // this is just used for plotting apparently..
                RooDataSet *mcsigallwdata = NULL;
                // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
                {
                  // string sigProcName = sigProcessNames[sigProcessIndex];

                  string name = "sig_" + inputSigProcName + "_mass_m" + numstring + "_" + catname;
                  RooDataSet *mcsigalldata = getDataSet(name);

                  mcsigallwdata = createOrAppendToDataset(mcsigalldata, hmass, weight, "mcsigallwdata" + numstring + catname, config.weightscale);

                } // loop over signal Feynman diagrams

                //----------

                //track some test datasets for later tests
#ifdef TEST_SINGLE_MASS_POINT
                if (mhs.at(i)==115)
                {
                  testdsets.push_back(mcsigwdata);
                  testdsetswrong.push_back(mcsigwrongwdata);
                  testdsetsall.push_back(mcsigallwdata);

                  //continue;
                }
#endif

#if 0
                // commented out for the moment
                //--------------------
                // plotting
                //--------------------
                if (allpdf != NULL)
                {
                  TCanvas *chfitall = new TCanvas;
                  TString plotnameall = TString("allvtx") + numstring + catnames.at(icat) + "_" + sigProcName + TString(".png");
                  RooPlot *hplotall = hmass->frame(Bins(100), Range("plotrange"));
                  mcsigallwdata->plotOn(hplotall);
                  allpdf->plotOn(hplotall, RooFit::LineColor(kBlue), Range("higgsrange"), NormRange("higgsrange"));
                  hplotall->SetTitle("");
                  hplotall->Draw();
                  chfitall->SaveAs(plotnameall);
                }
                //--------------------
#endif

              } // loop over mass points
            } // loop over signal input processes
    } // loop over categories
  }

  //----------------------------------------------------------------------
  void
  makeCombCatDatasets(RooRealVar *weight, RooRealVar *hmass, RooWorkspace *wextra)
  {
    for (UInt_t i = 0; i < mhs.size(); ++i)
    {
      std::stringstream numstringstr;
      numstringstr << mhs.at(i);
      TString numstring(numstringstr.str());
      //--------------------
      // take the datasets with the signal events from the input workspace
      //--------------------

      // create a new RooDataSet for signal MC, correct vertex assignment, all categories combined ?
      // (not clear where this is used here, seems only to be imported into the workspace ?)
      RooDataSet *rvdata = new RooDataSet(TString::Format("sig_mass_m%i_rv_combcat", mhs.at(i)), "", RooArgList(*hmass, *weight), weight->GetName());

      // create a new RooDataSet for signal MC, wrong vertex assignment, all categories combined ?
      // (not clear where this is used here, seems only to be imported into the workspace ?)
      RooDataSet *wvdata = new RooDataSet(TString::Format("sig_mass_wv_m%i_combcat", mhs.at(i)), "", RooArgList(*hmass, *weight), weight->GetName());

      // create a new RooDataSet for signal MC, both correct and wrong vertex assignments ? all categories combined ?
      // (not clear where this is used here, seems only to be imported into the workspace ?)
      RooDataSet *alldata = new RooDataSet(TString::Format("sig_mass_m%i_combcat", mhs.at(i)), "", RooArgList(*hmass, *weight), weight->GetName());

      //--------------------
      // merging (again) the input datasets
      BOOST_FOREACH(string sigProcName, inputSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              for (UInt_t icat = 0; icat < catnames.size(); ++icat)
              {
                // Note that the IC inputs seem not to distinguish between
                // right and wrong vertex assignment ?! (i.e. these
                // are fit together, not separately)

                if (config.useRightWrongVertex)
                {
                  RooDataSet *mcsigdata = getDataSet(TString("sig_") + sigProcName + "_mass_m" + numstring + TString("_rv_") + catnames.at(icat));
                  RooDataSet *mcsigwrongdata = getDataSet(TString("sig_") + sigProcName + "_mass_m" + numstring + TString("_wv_") + catnames.at(icat));

                  appendWeightedDataset(rvdata, mcsigdata, hmass, weight, 1.0);
                  appendWeightedDataset(wvdata, mcsigwrongdata, hmass, weight, 1.0);
                }
                else
                {
                  // for the sake of simplicity, just assign ALL events to the
                  // 'right' vertex datasets
                  RooDataSet *mcsigdata = getDataSet(TString("sig_") + sigProcName + "_mass_m" + numstring + TString("_") + catnames.at(icat));
                  appendWeightedDataset(rvdata, mcsigdata, hmass, weight, 1.0);
                }

                RooDataSet *mcsigalldata = getDataSet(TString("sig_") + sigProcName + "_mass_m" + numstring + TString("_") + catnames.at(icat));

                //rvdata->append(*mcsigdata);
                //wvdata->append(*mcsigwrongdata);
                //alldata->append(*mcsigalldata);
                appendWeightedDataset(alldata, mcsigalldata, hmass, weight, 1.0);
              } // loop over categories
            } // loop over signal process names

      // add these combined datasets to the workspace
      win->import(*rvdata);
      win->import(*wvdata);
      win->import(*alldata);

      // also to the extra workspace
      wextra->import(*rvdata);
      wextra->import(*wvdata);
      wextra->import(*alldata);
      //--------------------


    } // loop over categories

  }
  //----------------------------------------------------------------------

  /** @param categories: barrel is 0 and 1, endcap is 2 and 3 (in the case of 4 categories) */
  RooFormulaVar *
  makeEffAccVar(const string &varname, map<unsigned, map<string, RooAbsReal *> > &finalnormslides, RooRealVar &nuissancedeltaeffaccbarrel, const list<int> &categories)
  {
    const vector<string> *sigProcessNames = NULL;

    // use this when taking the normalization per
    // input signal process but the signal shape
    // from merged groups
    sigProcessNames = &inputSigProcessNames;

    // use this when taking the normalization from
    // the merged groups
    // sigProcessNames = &mergedSigProcessNames;

    RooArgList args(nuissancedeltaeffaccbarrel);

    // // cat 0 and 1 are barrel
    // for (unsigned cat = 0; cat <= 1; ++cat)
    BOOST_FOREACH(int cat, categories)
          {
            BOOST_FOREACH(string sigProcName, *sigProcessNames)
                  {
                    ASSERT(finalnormslides[cat][sigProcName] != NULL);
                    args.add(*finalnormslides[cat][sigProcName]);
                  }
          }

    // e.g. "@0*(@1+@2+@3+@4 + @5+@6+@7+@8)"
    string expr = "@0 * (" + makeArgSum(1, sigProcessNames->size() * categories.size()) + ")";
    cout << "XXX making effAccVar: " << expr << " sigProcessNames.size()=" << sigProcessNames->size() << endl;
    return new RooFormulaVar(varname.c_str(), "", expr.c_str(), args);
  }

  //----------------------------------------------------------------------

  /** @param categories: for barrel and 4 categories: highR9cat = 0, lowR9cat = 1
   *                         mixed  and 4 categories: highR9cat = 2, lowR9cat = 3
   *
   * builds a 'fraction' variable with nuisance parameter (e.g. the fraction of high R9 with respect to the total)
   * where the fraction is   numerator / ( numerator + other )
   *
   * can be used to build the R9 fraction variable and the pt fraction variable
   *
   * builds sums over all signal production mechanisms
   *
   * numeratorCats can e.g. be the categories corresponding to high R9 (in the given detector region)
   * and otherCats can e.g. be the categories corresponding to low R9 (in the same region)
   */
  RooFormulaVar *
  makeFracVar(const string &varname, map<unsigned, map<string, RooAbsReal *> > &finalnormslides, RooRealVar &deltaRnuisance,
                const list<int> &numeratorCats,  const list<int> &otherCats)
  {
    const vector<string> *sigProcessNames = NULL;

    // use this when taking the normalization per input signal process
    sigProcessNames = &inputSigProcessNames;
    // use this when taking the normalization per merged signal process
    // sigProcessNames = &mergedSigProcessNames;

    unsigned n = sigProcessNames->size();

    //    new RooFormulaVar(, "", "@0*(@1+@2+@3+@4)/(@1+@2+@3+@4 + @5+@6+@7+@8)", RooArgList(nuissancedeltar9fracbarrel,
    //          // cat 0
    //              *finalnormslides[0][sigProcessNames[0]], *finalnormslides[0][sigProcessNames[1]], *finalnormslides[0][sigProcessNames[2]], *finalnormslides[0][sigProcessNames[3]],
    //
    //              // cat 1
    //              *finalnormslides[1][sigProcessNames[0]], *finalnormslides[1][sigProcessNames[1]], *finalnormslides[1][sigProcessNames[2]], *finalnormslides[1][sigProcessNames[3]]));

    RooArgList args(deltaRnuisance);

    // // cat 0 and 1 are barrel
    // for (unsigned cat = 0; cat <= 1; ++cat)
    BOOST_FOREACH(int cat, numeratorCats)
    {
      BOOST_FOREACH(string sigProcName, *sigProcessNames)
      {
          ASSERT(finalnormslides[cat][sigProcName] != NULL);
          args.add(*finalnormslides[cat][sigProcName]);
      }
    }

    BOOST_FOREACH(int cat, otherCats)
    {
      BOOST_FOREACH(string sigProcName, *sigProcessNames)
      {
        ASSERT(finalnormslides[cat][sigProcName] != NULL);
        args.add(*finalnormslides[cat][sigProcName]);
      }
    }

    // n elements in the numerator
    // 2n elements in the denominator

    string exprNumerator = makeArgSum(1, n * numeratorCats.size()); // signal processes in the high R9 category
    string exprOther = makeArgSum(n * numeratorCats.size() + 1,
                                  n * (numeratorCats.size() + otherCats.size()) ); // signal processes in the low R9 category

    // e.g. @0*(@1+@2+@3+@4)/(@1+@2+@3+@4 + @5+@6+@7+@8)
    // i.e.
    string expr = "@0 * (" + exprNumerator + ") / (" + exprNumerator + "  +  " + exprOther + ")";
    return new RooFormulaVar(varname.c_str(), "", expr.c_str(), args);
  }

  //----------------------------------------------------------------------

  void
  fillEffAcc(int mh, unsigned icat, const string &inputSigProcName, map<unsigned, map<string, TH1F *> > &additionalUCSDcatEffaccHist, RooRealVar *IntLumi, RooRealVar &effacc, map<unsigned, map<string, TH1F *> > &effaccHistsX)
  {
    // TODO: note that the total cross section is NOT used as a variable
    // (but only its value -- this was the case already before we at UCSD
    // manipulated the macro), so this is currently NOT suitable
    // for determining RELATIVE limits (just absolute ones where
    // the theoretical cross section uncertainty does not enter)

    // IC also has separated this into the cross section per Feynman diagram.
    // for the moment, just add them together (no possibility to vary them
    // independently for the moment)

    double totalCrossSectionValue = getTotalCrossSection(mh);

    //compute acceptance * efficiency and right vertex fraction
    double effAccNumerator = 0;

    // now get the original datasets from the background input working space
    // (the weighted ones are merged already so do not contain the per
    // signal process information any more)

    if (config.useRightWrongVertex)
    {
      effAccNumerator =
      // right vertex
          getDataSet(Form("sig_%s_mass_m%d_rv_cat%d", inputSigProcName.c_str(), mh, icat))->sumEntries() +

          // wrong vertex
              getDataSet(Form("sig_%s_mass_m%d_wv_cat%d", inputSigProcName.c_str(), mh, icat))->sumEntries();
    }
    else
    {
      // just use the 'correct' ones
      // ASSERT(mcsigwrongwdata == NULL);
      effAccNumerator = getDataSet(Form("sig_%s_mass_m%d_cat%d", inputSigProcName.c_str(), mh, icat))->sumEntries();

    }

    double effAccDenominator = IntLumi->getVal() * totalCrossSectionValue;
    //double eaccden = (gfxsecs.at(i)+vbfxsecs.at(i)+vhxsecs.at(i))*1e6;
    //double eaccden = (vbfxsecs.at(i)+vhxsecs.at(i))*1e6;

    // set the value of efficiency times acceptance
    assignRatio(effacc, effAccNumerator, effAccDenominator);

    cout << "UUU mass=" << mh << " cat=" << icat << " proc=" << inputSigProcName << " effAccNumerator=" << effAccNumerator << " effAccDenominator=" << effAccDenominator << " ratio=" << effAccNumerator / effAccDenominator << endl;
    cout << "UUU0 totalCrossSectionValue=" << totalCrossSectionValue << endl;

    // make sure we do NOT fill values into a bin more than once
    if ((mh % ((int) (config.massInterpolationWidth + 0.5))) == 0)
    {
      // fill the efficiency (and its asymmetric error) into
      // the efficiency histogram for category 'numInclusiveCategories'
      cout << "FILLING " << Form("effaccCat%dHist", icat) << " at mh=" << mh << " value=" << effacc.getVal() << endl;

      TH1F *histo = NULL;

      // special treatment for UCSD-specific categories
      if (icat >= config.numInclusiveCategories)
        // this is a VBF category
        // MUST USE INPUT SIGNAL PROCESSES HERE
        histo = additionalUCSDcatEffaccHist[icat][inputSigProcName];
      else
        // non-VBF category
        histo = effaccHistsX[icat][inputSigProcName];

      // fill the histogram
      histo->Fill(mh, effacc.getVal());

      // set the error in this bin (average of positive and negative error bars)
      double avgerror = (effacc.getErrorLo() + effacc.getErrorHi()) / 2.0;
      histo->SetBinError(histo->FindFixBin(mh), avgerror);

    } // if this is a mass point used in the histogram


  }

  //----------------------------------------------------------------------

  /** plot the comparison between the fitted functional form and the
   *  mass points after fitting the mass spectrum
   *
   *  @param plotNamePart is typically "rightvtx" or "wrongvtx"
   */
  void
  drawFittedMass(const string &plotNamePart, int mh, unsigned icat, const string &mergedSigProcName, int bins, RooRealVar *hmass, RooDataSet *mcsigwdata, RooAbsPdf *rightpdf, RooGaussian *g1, RooGaussian *g2, RooGaussian *g3)
  {
    //plot fit results for this mass point
    TCanvas *chfit = new TCanvas;
    TString plotname = "mass_";
    plotname += plotNamePart;
    plotname += "_";

    plotname += "m" + boost::lexical_cast<string>(mh) + "_";
    plotname += catnames.at(icat) + "_" + mergedSigProcName + TString(".png");

    RooPlot *hplot = hmass->frame(Bins(bins), Range("plotrange"));

    mcsigwdata->plotOn(hplot);
    rightpdf->plotOn(hplot, Components(*g1), LineColor(kOrange), Range("higgsrange"), NormRange("higgsrange"));
    rightpdf->plotOn(hplot, Components(*g2), LineColor(kMagenta), Range("higgsrange"), NormRange("higgsrange"));

    if (g3 != NULL)
      rightpdf->plotOn(hplot, Components(*g3), LineColor(kRed), Range("higgsrange"), NormRange("higgsrange"));

    rightpdf->plotOn(hplot, RooFit::LineColor(kBlue), Range("higgsrange"), NormRange("higgsrange"));
    hplot->SetTitle(TString("mass fit comparison ") + plotNamePart + " " + catnames.at(icat) + " " + mergedSigProcName);
    hplot->Draw();

    // add a legend (see e.g. http://root.cern.ch/phpBB3/viewtopic.php?p=31694 and
    // the attached code http://root.cern.ch/phpBB3/download/file.php?id=6019 )
    //    for (int i=0; i<hplot->numItems(); i++)
    //    {
    //      TString obj_name = hplot->nameOf(i);
    //      if (obj_name=="")
    //        continue;
    //      cout << Form("ITEM: %d. '%s'\n",i,obj_name.Data());
    //    }

    // the items are of the form (example):
    // combh_Norm[CMS_hgg_mass]_Comp[g1]_Range[higgsrange]_NormRange[higgsrange]
    // or they seem to be items 1..3

    TLegend *leg = new TLegend(0.15, 0.5, 0.3, 0.75);

    leg->AddEntry(hplot->getObject(1), "Gaussian 1", "l");
    leg->AddEntry(hplot->getObject(2), "Gaussian 2", "l");
    if (g3 != NULL)
      leg->AddEntry(hplot->getObject(3), "Gaussian 3", "l");

    leg->Draw();

    chfit->SaveAs(plotname);

  }
  //----------------------------------------------------------------------

  void
  fitSinglePoint(int mh, unsigned icat, const string &catname,

  const string &mergedSigProcName, RooRealVar *mnom, RooRealVar *hmass,

  map<unsigned, map<string, map<unsigned, TH1F *> > > &fitparmhistsX,

  TH2F *histfitstatus, TH2F *histfitstatuswrong,

  std::vector<RooRealVar*> &fitparms, std::vector<double> &fitparmsinit,

  RooAddPdf &combh, RooAddPdf &combhwrong, RooAddPdf &combhmin,

  RooAddPdf &combhvtx, RooAddPdf &combhvtxsimple, RooAddPdf &combhvtxmin,

  RooRealVar &f1, RooRealVar &dm1, RooRealVar &sigma1, RooGaussian &g1, RooGaussian &wg1,

  RooRealVar &f2, RooRealVar &dm2, RooRealVar &sigma2, RooGaussian &g2, RooGaussian &wg2,

  RooRealVar &dm3, RooRealVar &sigma3, RooGaussian &g3,

  RooRealVar &fracright)
  {
    std::stringstream numstringstr;
    numstringstr << mh;
    TString numstring(numstringstr.str());
    TString vbfstring = numstring;

    //define aggregate weight, so far using xsec, pileup, pt-reweighting and efficiency scale factors
    //RooFormulaVar totweight("totweight","totweight","xsecweight*puweight(ngenvtx-1)*ptweight(genhpt,isgf)*effweight(iseb1,iseb2)",RooArgList(xsecweight,ngenvtx,genhpt,isgf,iseb1,iseb2));

    //          RooDataSet *mcsigdata = NULL;
    //          RooDataSet *mcsigwdata = NULL;
    //          RooDataSet *mcsigwrongdata = NULL;
    //          RooDataSet *mcsigwrongwdata = NULL;


    //--------------------
    // set higgs mass hypothesis and fit range
    mnom->setVal(mh);
    hmass->setRange("higgsrange", 100.0, config.massmax);
    hmass->setRange("plotrange", TMath::Max(100.0, mh - 30.0), TMath::Min(config.massmax, mh + 20.0));

    //reset parameters to initial values
    ASSERT(fitparms.size() == fitparmsinit.size());
    for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
    {
      fitparms.at(iparm)->setVal(fitparmsinit.at(iparm));
    }
    // cout << "AFTER RESET: icat=" << icat << " fitparms[FITPARAM_INDEX_EFFACC]= " << fitparms.at(FITPARAM_INDEX_EFFACC)->getVal() << endl;

    // the PDF to be fitted for the events with right vertex assignment
    RooAbsPdf *rightpdf = 0;

    // the PDF to be fitted for the events with wrong vertex assignment
    RooAbsPdf *wrongpdf = 0;
    RooAbsPdf *allpdf = 0;

    Bool_t usecomplexmodel = complexset.count(catnames.at(icat));
    Bool_t usesimplemodel = simpleset.count(catnames.at(icat));

    if (usecomplexmodel)
    {
      rightpdf = &combhmin;
      wrongpdf = &combhwrong;
      allpdf = &combhvtx;
      f1.setRange(0.6, 0.95);
      sigma2.setVal(4.0);
      dm1.setRange(-1.0, 0.05);
      dm1.setVal(-0.01);
    }
    else if (usesimplemodel)
    {
      rightpdf = &g1;
      wrongpdf = &wg1;
      allpdf = &combhvtxsimple;
      sigma1.setVal(5.0);
    }
    else
    {
      rightpdf = &combhmin;
      wrongpdf = &wg1;
      allpdf = &combhvtxmin;
      sigma1.setVal(2.0);
      sigma2.setVal(3.0);
      f1.setRange(0.0, 1.0);
    }
    
    dm2.setRange(-5.0, 1.0);

    //if (mh == 110 && !usecomplexmodel)
    //{
    //  dm2.setRange(-8.0, 5.0);
    //}
    //if (mh == 110)
    //{
    //  dm1.setVal(0.0);
    //}
    //if (mh == 150)
    //{
    //  dm1.setVal(-0.5);
    //  dm2.setVal(-1.5);
    //}

    RooDataSet *mcsigwdata = getMergedWeightedDataSet("mcsigwdata" + boost::lexical_cast<string>(mh) + catname + "_" + mergedSigProcName);
    RooDataSet *mcsigwrongwdata = NULL;
    if (config.useRightWrongVertex)
      mcsigwrongwdata = getMergedWeightedDataSet("mcsigwrongwdata" + boost::lexical_cast<string>(mh) + catname + "_" + mergedSigProcName);

    //----------------------------------------
    // IT'S HERE WHERE FITTING IS HAPPENING !!!
    // perform fits
    //----------------------------------------
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    rightpdf->Print();
    rightpdf->Print("V");
    rightpdf->printCompactTree();

    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    // fit the (mass) distribution of events
    // with CORRECT vertex assignment
    RooFitResult *fitres = rightpdf->fitTo(*mcsigwdata, Strategy(0), Minimizer("Minuit2", ""), Minos(kFALSE), SumW2Error(kFALSE), Save(kTRUE), NumCPU(8));
    if (!fitres->status())
      // try again
      fitres = rightpdf->fitTo(*mcsigwdata, Strategy(1), Minimizer("Minuit2", ""), Minos(kFALSE), SumW2Error(kFALSE), Save(kTRUE), NumCPU(8));

    // fit the (mass) distribution of events
    // with WRONG vertex assignment
    RooFitResult *fitreswrong = NULL;
    if (config.useRightWrongVertex)
    {
      ASSERT(mcsigwrongwdata != NULL);
      fitreswrong = wrongpdf->fitTo(*mcsigwrongwdata, Strategy(1), Minimizer("Minuit2", ""), Minos(kFALSE), SumW2Error(kFALSE), Save(kTRUE), NumCPU(8));
      if (!fitreswrong->status())
        // try again
        fitreswrong = wrongpdf->fitTo(*mcsigwrongwdata, Strategy(0), Minimizer("Minuit2", ""), Minos(kFALSE), SumW2Error(kFALSE), Save(kTRUE), NumCPU(8));
    }
    //----------------------------------------

    dm2.removeRange();

    // testing
    //if (usecomplexmodel)
    //{
    //  // more Gaussians in the sum for some of the categories
    //  if (TMath::Abs(sigma2.getVal()) > TMath::Abs(sigma3.getVal()) || f2.getVal() < 0.5)
    //  {
    //    double sigma2tmp = sigma2.getVal();
    //    double dm2tmp = dm2.getVal();

    //    sigma2.setVal(sigma3.getVal());
    //    dm2.setVal(dm3.getVal());
    //    sigma3.setVal(sigma2tmp);
    //    dm3.setVal(dm2tmp);
    //    f2.setVal(1.0 - f2.getVal());
    //  }
    //} // if usecomplexmodel


    if (!usesimplemodel)
    {
      //printf("testing sigmas, sigma1 = %5f, sigma2 = %5f\n",sigma1.getVal(),sigma2.getVal());
      if (TMath::Abs(sigma1.getVal()) > TMath::Abs(sigma2.getVal()) || f1.getVal() < 0.5 || dm1.getVal() < dm2.getVal())
      {
        //printf("swapping gaussians\n");
        double sigma1tmp = sigma1.getVal();
        double dm1tmp = dm1.getVal();

        sigma1.setVal(sigma2.getVal());
        dm1.setVal(dm2.getVal());
        sigma2.setVal(sigma1tmp);
        dm2.setVal(dm1tmp);
        f1.setVal(1.0 - f1.getVal());
      }

    } // if ! usesimplemodel
    //----------------------------------------


    // the fraction of signal events with their vertex correctly assigned ?
    //double fright = mcsigwdata->sumEntries()/(mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
    double fractionRightNumerator = mcsigwdata->sumEntries();
    double fractionRightDenominator = 0;

    if (config.useRightWrongVertex)
      fractionRightDenominator = mcsigwdata->sumEntries() + mcsigwrongwdata->sumEntries();
    else
    {
      ASSERT(mcsigwrongwdata == NULL);
      fractionRightDenominator = mcsigwdata->sumEntries();
    }

    cout << "mass=" << mh << " cat=" << icat << " frightden=" << fractionRightDenominator << " frightnum=" << fractionRightNumerator << endl;

    // fraction of events with correctly assigned vertices
    assignRatio(fracright, fractionRightNumerator, fractionRightDenominator);

    //----------

    // correct negative resolution terms which screw up interpolation
    // (convert to absolute values)
    fitparms.at(3)->setVal(TMath::Abs(fitparms.at(3)->getVal()));
    fitparms.at(4)->setVal(TMath::Abs(fitparms.at(4)->getVal()));
    fitparms.at(5)->setVal(TMath::Abs(fitparms.at(5)->getVal()));
    fitparms.at(10)->setVal(TMath::Abs(fitparms.at(10)->getVal()));
    fitparms.at(11)->setVal(TMath::Abs(fitparms.at(11)->getVal()));

    //----------------------------------------
    //fill histograms with fit parameters in 5 (10 ?) GeV steps
    // (fill fitparmhists from fitparms)
    //----------------------------------------
    if ((mh % ((int) (config.massInterpolationWidth + 0.5))) == 0)
    {
      for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
      {
        // fitparmhists[icat*numFitParms + iparm]
        // WATCH OUT NOT TO FILL THE SAME HISTOGRAM MORE THAN ONCE !
        fitparmhistsX[icat][mergedSigProcName][iparm]->Fill(mh, fitparms.at(iparm)->getVal());

        fitparmhistsX[icat][mergedSigProcName][iparm]->SetBinError(fitparmhistsX[icat][mergedSigProcName][iparm]->FindFixBin(mh), fitparms.at(iparm)->getError());
        if (fitparms.at(iparm)->hasAsymError())
        {
          // calculate the average error from the asymmetric error bars
          double avgerror = (fitparms.at(iparm)->getErrorLo() + fitparms.at(iparm)->getErrorHi()) / 2.0;
          fitparmhistsX[icat][mergedSigProcName][iparm]->SetBinError(fitparmhistsX[icat][mergedSigProcName][iparm]->FindFixBin(mh), avgerror);
        }
        printf("filling histogram named: %s, variable named %s, proc=%s, val = %5f, err = %5f\n", fitparmhistsX[icat][mergedSigProcName][iparm]->GetName(), fitparms.at(iparm)->GetName(), mergedSigProcName.c_str(), fitparms.at(iparm)->getVal(), fitparms.at(iparm)->getError());
      } // loop over parameters

      histfitstatus->Fill(mh, icat, fitres->status());
      if (config.useRightWrongVertex)
      {
        ASSERT(fitreswrong != NULL);
        histfitstatuswrong->Fill(mh, icat, fitreswrong->status());
      }

    } // if mass point

    // draw the mass distribution for the events with correctly assigned vertex
    drawFittedMass("rightvtx", mh, icat, mergedSigProcName, 100, hmass, mcsigwdata, rightpdf, &g1, &g2, &g3);

    if (config.useRightWrongVertex)
    {
      // also draw the mass distribution comparison for the events with a wrongly assigned veretx
      drawFittedMass("wrongvtx", mh, icat, mergedSigProcName, 40, hmass, mcsigwrongwdata, wrongpdf, &wg1, &wg2, NULL);
    }

    // commented out for the moment as mcsigwrongwdata may be NULL
    // printf ("right = %5f, wrong = %5f, all = %5f, right+wrong = %5f\n",mcsigwdata->sumEntries(),mcsigwrongwdata->sumEntries(), mcsigallwdata->sumEntries(),mcsigwdata->sumEntries()+mcsigwrongwdata->sumEntries());
    printf("data weights = %5e\n", mcsigwdata->sumEntries());

    // commented out for the moment as fitreswrong may be NULL
    // printf("mass = %i, status = %i, statuswrong = %i\n",mhs.at(i),fitres->status(),fitreswrong->status());
  }

  //----------------------------------------------------------------------

  list<int> getBarrelCategories()
  {
    switch (config.numInclusiveCategories)
    {
      case 4:
        return list_of(0)(1);

      case 8:
        return list_of(0)(1) (4)(5);

      case 12:
        return list_of(0)(1) (4)(5) (8)(9);

      default:
        cerr << "numInclusiveCategories = " << config.numInclusiveCategories << " is not yet supported here" << endl;
        exit(1);
    }
  }

  //----------------------------------------------------------------------

  /** @return a list of the non-Barrel-Barrel categories */
  list<int> getEndcapCategories()
  {
    switch (config.numInclusiveCategories)
    {
      case 4:
        return list_of(2)(3);

      case 8:
        return list_of(2)(3) (6)(7);

      case 12:
        return list_of(2)(3) (6)(7) (10)(11);

      default:
        cerr << "numInclusiveCategories = " << config.numInclusiveCategories << " is not yet supported here" << endl;
        exit(1);
    }
  }

  //----------------------------------------------------------------------

public:

  /** perform the calculations */
  void
  run()
  {
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    //--------------------

    gStyle->SetOptStat(1110);

    RooRealVar *hmass = win->var("CMS_hgg_mass");

    //RooRealVar *hmass = win->var("mass");
    hmass->setRange(100, config.massmax);
    hmass->setBins(2.0 * (config.massmax - 100.0));
    hmass->SetTitle("m_{#gamma#gamma}");
    hmass->setUnit("GeV");

    std::vector<RooRealVar*> fitparms;
    std::vector<double> fitparmsinit;

    //define signal pdf variables

    const double dm1init = 0.5;
    RooRealVar dm1("dm1", "", dm1init, -12.0, 5.0);
    dm1.removeRange();

    const double dm2init = -1.0;
    RooRealVar dm2("dm2", "", dm2init, -12.0, 5.0);
    dm2.removeRange();

    const double dm3init = -2.0;
    RooRealVar dm3("dm3", "", dm3init, -9.0, 5.0);
    //dm3.removeRange();

    RooFormulaVar mean1("mean1", "", "@0+@1", RooArgList(*mnom, dm1));
    RooFormulaVar mean2("mean2", "", "@0+@1", RooArgList(*mnom, dm2));
    RooFormulaVar mean3("mean3", "", "@0+@1", RooArgList(*mnom, dm3));

    const double sigma1init = 1.0;
    RooRealVar sigma1("sigma1", "", sigma1init, 0.5, 5.0);
    sigma1.removeRange();

    const double sigma2init = 1.2;
    RooRealVar sigma2("sigma2", "", sigma2init, 0.8, 7.0);
    sigma2.removeRange();

    const double sigma3init = 2.0;
    RooRealVar sigma3("sigma3", "", sigma3init, 1.0, 10.0);
    sigma3.removeRange();

    const double f1init = 0.60;
    RooRealVar f1("f1", "", f1init, 0.52, 1.0);

    const double f2init = 0.90;
    RooRealVar f2("f2", "", f2init, 0.0, 1.0);

    RooGaussian g1("g1", "", *hmass, mean1, sigma1);
    RooGaussian g2("g2", "", *hmass, mean2, sigma2);
    RooGaussian g3("g3", "", *hmass, mean3, sigma3);

    //              name   title    pdfs             coefficients for pdfs
    //                                                                use recursive fractions
    RooAddPdf combh("combh", "", RooArgList(g1, g2, g3), RooArgList(f1, f2), kTRUE);
    RooAddPdf combhmin("combhmin", "", RooArgList(g1, g2), RooArgList(f1), kTRUE);
    // RooGaussian &combh = g1;

    const double wdm1init = 0.5;
    RooRealVar wdm1("wdm1", "", wdm1init, -12.0, 5.0);
    wdm1.removeRange();

    const double wdm2init = -0.7;
    RooRealVar wdm2("wdm2", "", wdm2init, -12.0, 5.0);
    wdm2.removeRange();

    //--------------------

    RooFormulaVar wmean1("wmean1", "", "@0+@1", RooArgList(*mnom, wdm1));
    RooFormulaVar wmean2("wmean2", "", "@0+@1", RooArgList(*mnom, wdm2));

    const double wsigma1init = 4.0;
    RooRealVar wsigma1("wsigma1", "", wsigma1init, 0.0, 10.0); //2.0
    wsigma1.removeRange();

    const double wsigma2init = 2.0;
    RooRealVar wsigma2("wsigma2", "", wsigma2init, 0.0, 10.0); //3.0
    wsigma2.removeRange();

    // some fraction for adding wg1 and wg2
    const double wf1init = 0.6;
    RooRealVar wf1("wf1", "", wf1init, 0.52, 1.0);

    //                      indep var      sigma
    //                              mean
    RooGaussian wg1("wg1", "", *hmass, wmean1, wsigma1);
    RooGaussian wg2("wg2", "", *hmass, wmean2, wsigma2);

    // a sum of two Gaussians
    RooAddPdf combhwrong("combhwrong", "", RooArgList(wg1, wg2), RooArgList(wf1), kTRUE);
    //RooGaussian &combhwrong = wg1;

    //--------------------

    // fraction of signal events assigned to the
    // correct vertex (seems to be fitted later on ?)
    const double fracrightinit = 0.9;
    RooRealVar fracright("fracright", "fracright", fracrightinit, 0.0, 1.0);

    const double effaccinit = 0.9;
    RooRealVar effacc("effacc", "effacc", effaccinit, 0.0, 1.0);

    // histogram for keeping the efficiency of category 8
    // see also initialization of fitparmhists
    //
    // first index is the category number, second index
    // is the name of the production mechanism
    // std::vector<TH1F *> additionalUCSDcatEffaccHist;
    map<unsigned, map<string, TH1F *> > additionalUCSDcatEffaccHist;

    // AH: we could actually use funceffacc_catX instead...
    for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
    {
      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            //for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              additionalUCSDcatEffaccHist[j][inputSigProcName] = makeMassFunctionHist(Form("effaccCat%dHist_%s", j, inputSigProcName.c_str()));
            }
    } // loop over all categories

    //signal pdf
    // sum of right and wrong vertex distributions
    RooAddPdf combhvtx("combhvtx", "combhvtx", RooArgList(combh, combhwrong), RooArgList(fracright));
    RooAddPdf combhvtxmin("combhvtxmin", "combhvtxmin", RooArgList(combhmin, wg1), RooArgList(fracright));
    RooAddPdf combhvtxsimple("combhvtxsimple", "combhvtxsimple", RooArgList(g1, wg1), RooArgList(fracright));

    //----------------------------------------
    // parameters to fit
    //----------------------------------------
    fitparms.push_back(&dm1); // 0: delta mean of 1st Gaussian (delta w.r.t. generated Higgs mass)
    fitparms.push_back(&dm2); // 1: delta mean of 2nd Gaussian (delta w.r.t. generated Higgs mass)
    fitparms.push_back(&dm3); // 2: delta mean of 3rd Gaussian (delta w.r.t. generated Higgs mass)
    fitparms.push_back(&sigma1); // 3: sigma of 1st Gaussian
    fitparms.push_back(&sigma2); // 4: sigma of 2nd Gaussian
    fitparms.push_back(&sigma3); // 5: sigma of 3rd Gaussian
    fitparms.push_back(&f1); // 6
    fitparms.push_back(&f2); // 7
    fitparms.push_back(&wdm1); // 8
    fitparms.push_back(&wdm2); // 9
    fitparms.push_back(&wsigma1); // 10
    fitparms.push_back(&wsigma2); // 11
    fitparms.push_back(&wf1); // 12
    fitparms.push_back(&fracright); // 13

    // not used any more in this array, treated separately because not fitted but calculated
    //    fitparms.push_back(&effacc); // FITPARAM_INDEX_EFFACC: efficiency times acceptance ? = effacc
    //                                 //                        IS CALCULATED FROM SUM OF MC WEIGHTS DIVIDED
    //                                 //                        BY CROSS SECTION TIMES LUMI, I.E. NOT FITTED

    // initial guesses for the fit parameters ?!
    fitparmsinit.push_back(dm1init);
    fitparmsinit.push_back(dm2init);
    fitparmsinit.push_back(dm3init);
    fitparmsinit.push_back(sigma1init);
    fitparmsinit.push_back(sigma2init);
    fitparmsinit.push_back(sigma3init);
    fitparmsinit.push_back(f1init);
    fitparmsinit.push_back(f2init);
    fitparmsinit.push_back(wdm1init);
    fitparmsinit.push_back(wdm2init);
    fitparmsinit.push_back(wsigma1init);
    fitparmsinit.push_back(wsigma2init);
    fitparmsinit.push_back(wf1init);
    fitparmsinit.push_back(fracrightinit);
    // fitparmsinit.push_back(effaccinit);

    RooRealVar *weight = new RooRealVar("weight", "", 1.0);
    weight->removeRange();

    RooRealVar *IntLumi = win->var("IntLumi");

    //loop over categories for fits
    RooWorkspace *wextra = new RooWorkspace("wextra", "");

    makeCombCatDatasets(weight, hmass, wextra);

    //catnames.push_back("singlecat");


    const unsigned numFitParms = fitparms.size();

    //----------------------------------------
    // define histograms to keep track of fit parameters for each mass point
    //----------------------------------------

    // indexing is:
    //   first index is the category number
    //   second index is the signal process name
    //   third index is the number of the parameter to be fit

    // TH1F **fitparmhists = new TH1F*[numCases() * numFitParms];
    map<unsigned, map<string, map<unsigned, TH1F *> > > fitparmhistsX;
    map<unsigned, map<string, TH1F *> > effaccHistsX;

    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      BOOST_FOREACH(string mergedSigProcName, mergedSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
              {
                TString histname = TString("hist") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat) + mergedSigProcName;
                fitparmhistsX[icat][mergedSigProcName][iparm] = makeMassFunctionHist((const char*) histname);
              } // loop over fit parameters


            } // loop over signal processes

      // also create histograms for effacc
      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            {
              TString histname = TString("hist") + "effacc" + catnames.at(icat) + inputSigProcName;
              effaccHistsX[icat][inputSigProcName] = makeMassFunctionHist((const char*) histname);
            }

    } // loop over categories

    // overview over which fits (as function of mass hypothesis and parameter ?) failed ?!
    TH2F *histfitstatus = new TH2F("histfitstatus", "fit failures for right vertex",

    (config.massmax - 100.0 - config.massInterpolationWidth) / config.massInterpolationWidth, 100 + config.massInterpolationWidth / 2, config.massmax - config.massInterpolationWidth / 2,

    numCases(), -0.5, numCases() - 0.5);

    TH2F *histfitstatuswrong = new TH2F("histfitstatuswrong", "fit failures for wrong vertex",

    (config.massmax - 100.0 - config.massInterpolationWidth) / config.massInterpolationWidth, 100 + config.massInterpolationWidth / 2, config.massmax - config.massInterpolationWidth / 2,

    numCases(), -0.5, numCases() - 0.5);

#ifdef TEST_SINGLE_MASS_POINT
    std::vector<RooDataSet*> testdsets;
    std::vector<RooDataSet*> testdsetswrong;
    std::vector<RooDataSet*> testdsetsall;
#endif

    makeWeightedDataSets(weight, hmass, NULL);

    // loop over categories for fits
    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      string catname = (const char *) catnames.at(icat);

      // TODO: can we move the mass loop outside the loop over signal processes ?

      // BOOST_FOREACH(string sigProcName, mergedSigProcessNames)
      // BOOST_FOREACH(string sigProcName, sigProcessNames)
      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            {
              for (UInt_t i = 0; i < mhs.size(); ++i)
              {
                int mh = mhs.at(i);
                fillEffAcc(mh, icat, inputSigProcName, additionalUCSDcatEffaccHist, IntLumi, effacc, effaccHistsX);
              } // loop over Higgs masses

            } // loop over input signal Feynman diagrams

      BOOST_FOREACH(string mergedSigProcName, mergedSigProcessNames)
            {
              // loop over mass points for fits
              for (UInt_t i = 0; i < mhs.size(); ++i)
              {
                int mh = mhs.at(i);

                // perform the fit to the signal mass distribution
                fitSinglePoint(mh, icat, catname,

                mergedSigProcName, mnom, hmass,

                fitparmhistsX, histfitstatus, histfitstatuswrong,

                fitparms, fitparmsinit, combh, combhwrong, combhmin,

                combhvtx, combhvtxsimple, combhvtxmin,

                f1, dm1, sigma1, g1, wg1,

                f2, dm2, sigma2, g2, wg2,

                dm3, sigma3, g3,

                fracright);
              } // loop over merged signal Feynman diagrams


            } // loop over Higgs masses
    } // loop over categories

    //----------------------------------------------------------------------
    // construct interpolation construction for each fit parameter in each category
    // (fill fitparmfuncs from fitparmhists via fitparmdatas)
    //----------------------------------------------------------------------

    /** first index is category, second index is signal process name, third index is parameter number */
    map<unsigned, map<string, map<unsigned, RooDataHist *> > > fitparmdatasX;
    map<unsigned, map<string, map<unsigned, RooHistFunc*> > > fitparmfuncsX;

    map<unsigned, map<string, RooDataHist *> > effaccDatasX;
    map<unsigned, map<string, RooHistFunc*> > effaccFuncsX;

    // RooDataHist **fitparmdatas = new RooDataHist*[numCases()*numFitParms];
    // RooHistFunc **fitparmfuncs = new RooHistFunc*[numCases()*numFitParms];
    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      BOOST_FOREACH(string mergedSigProcName, mergedSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
              {
                TString dataname = TString("data") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat) + "_" + mergedSigProcName;
                TString funcname = TString("func") + TString(fitparms.at(iparm)->GetName()) + catnames.at(icat) + "_" + mergedSigProcName;

                fitparmdatasX[icat][mergedSigProcName][iparm] = new RooDataHist(dataname, dataname, RooArgList(*mnom), fitparmhistsX[icat][mergedSigProcName][iparm]);
                fitparmfuncsX[icat][mergedSigProcName][iparm] = new RooHistFunc(funcname, funcname, RooArgList(*mnom), *fitparmdatasX[icat][mergedSigProcName][iparm], 1);
              } // loop over fit parameters


            } // loop over signal Feynman diagrams

      // also create the histogram functions for the efficiency*acceptance
      BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
            {
              TString dataname = TString("data") + TString("effacc") + catnames.at(icat) + "_" + inputSigProcName;
              TString funcname = TString("func") + TString("effacc") + catnames.at(icat) + "_" + inputSigProcName;

              effaccDatasX[icat][inputSigProcName] = new RooDataHist(dataname, dataname, RooArgList(*mnom), effaccHistsX[icat][inputSigProcName]);
              effaccFuncsX[icat][inputSigProcName] = new RooHistFunc(funcname, funcname, RooArgList(*mnom), *effaccDatasX[icat][inputSigProcName], 1);
            }

    } // loop over categories

    //----------------------------------------------------------------------
    // plot evolution of each fit parameter in each category as function of higgs mass
    //----------------------------------------------------------------------
    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      BOOST_FOREACH(string sigProcName, mergedSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
              {
                TString plotname = TString("func_") + TString(fitparms.at(iparm)->GetName()) + "_" + catnames.at(icat) + "_" + sigProcName + TString(".png");

                TCanvas *cfunctest = new TCanvas;
                // fitparmhists[icat*nparms + iparm]->Draw();
                RooPlot *hploteffacc = mnom->frame(Bins(90), Range(105, 145));

                cout << "----------------------------------------" << endl;
                cout << "printing icat=" << icat << " iparm=" << iparm << endl;
                fitparmdatasX[icat][sigProcName][iparm]->Print();
                cout << "sum entries:" << fitparmdatasX[icat][sigProcName][iparm]->sumEntries() << endl;
                cout << "----------------------------------------" << endl;

                // TODO: workaround to avoid crash: disabled drawing of errors
                fitparmdatasX[icat][sigProcName][iparm]->plotOn(hploteffacc, RooFit::DataError(RooAbsData::None));
                fitparmfuncsX[icat][sigProcName][iparm]->plotOn(hploteffacc, RooFit::LineColor(kBlue), RooFit::DataError(RooAbsData::None));
                hploteffacc->SetTitle("");
                hploteffacc->GetYaxis()->SetTitle(TString(fitparms.at(iparm)->GetName()) + " " + catnames.at(icat) + " " + sigProcName);
                hploteffacc->Draw();
                cfunctest->SaveAs(plotname);
                delete cfunctest;
              } // loop over parameters
            } // loop over signal Feynman diagrams
    } // loop over categories

    //----------------------------------------------------------------------

    TCanvas *cfitstatus = new TCanvas;
    histfitstatus->Draw("TEXT" /* "COL" */);
    cfitstatus->SaveAs("fitstatus_right_vertex.png");

    TCanvas *cfitstatuswrong = new TCanvas;
    histfitstatuswrong->Draw("TEXT" /*"COL"*/);
    cfitstatuswrong->SaveAs("fitstatus_wrong_vertex.png");

    RooRealVar nuissancedeltafracright("CMS_hgg_nuissancedeltafracright", "", 1.0, 0.1, 10.0);
    nuissancedeltafracright.setConstant();

    //   RooRealVar *nuissancedeltaescalebarrel = new RooRealVar("nuissancedeltaescalebarrel","",0.0,-0.1,0.1);
    //   nuissanceescalebarrel->setConstant();
    //   RooRealVar *nuissancedeltaescaleendcap = new RooRealVar("nuissancedeltaescaleendcap","",0.0,-0.1,0.1);
    //   nuissanceescaleendcap->setConstant();
    //
    //   RooFormulaVar *deltambarrel = new RooRealVar("deltambarrel","","@0*@1",RooArgList(mnom,*nuissancedeltaescalebarrel));
    //   RooFormulaVar *deltammixedhighr9 = new RooFormulaVar("deltammixedhighr9","","@0*(0.7*(sqrt((@0-1.0)*(@1-1.0))-1.0) + 0.3*@2)",RooArgList(mnom,*nuissancedeltaescalebarrel,*nuissancedeltaescaleendcap));
    //   RooFormulaVar *deltammixedmixedr9 = new RooFormulaVar("deltammixedmixedr9","","@0*(0.8*(sqrt((@0-1.0)*(@1-1.0))-1.0) + 0.2*@2)",RooArgList(mnom,*nuissancedeltaescalebarrel,*nuissancedeltaescaleendcap));
    //   RooRealVar **nuissancedeltamcors = new RooRealVar*[8];
    //   nuissancedeltamcors[0] = deltambarrel;
    //   nuissancedeltamcors[1] = deltambarrel;
    //   nuissancedeltamcors[2] = deltammixedhighr9;
    //   nuissancedeltamcors[3] = deltammixedmixedr9;
    //   nuissancedeltamcors[4] = deltambarrel;
    //   nuissancedeltamcors[5] = deltambarrel;
    //   nuissancedeltamcors[6] = deltammixedhighr9;
    //   nuissancedeltamcors[7] = deltammixedmixedr9;

    //----------------------------------------------------------------------
    // consistency checks
    //----------------------------------------------------------------------
    if (config.smearingv.size() != catnames.size())
    {
      cerr << "smearingv.size() must be the same as catnames.size()" << endl;
      ASSERT(config.smearingv.size() == catnames.size());
    }

    //----------------------------------------------------------------------
    // define final pdfs in each category
    //----------------------------------------------------------------------
    map<unsigned, RooRealVar *> nuissancedeltasmears;
    map<unsigned, map<string, RooAddPdf *> > combhvtxslides;
    map<unsigned, map<string, RooAddPdf *> > combhvtxminslides;
    map<unsigned, map<string, RooAddPdf *> > combhvtxsimpleslides;
    map<unsigned, map<string, RooAbsPdf *> > finalpdfslides;
    map<unsigned, map<string, RooAbsReal *> > finalnormslides;

    // this works on the merged signal process names
    this->makeFinalPDF(numFitParms, fitparmfuncsX, effaccFuncsX,

    hmass, nuissancedeltafracright,
    // output variables
        nuissancedeltasmears, combhvtxslides, combhvtxminslides, combhvtxsimpleslides, finalpdfslides, finalnormslides // efficiency times acceptance ?!
    );

    //----------------------------------------------------------------------
    // comparison tests of interpolated shape vs actual fit
    //----------------------------------------------------------------------
    for (UInt_t icat = 0; icat < catnames.size(); ++icat)
    {
      BOOST_FOREACH(string sigProcName, mergedSigProcessNames)
            // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
            {
              //nuissancedeltasmears[icat]->setVal(-smearingv[icat]);

              for (UInt_t iparm = 0; iparm < fitparms.size(); ++iparm)
                fitparms.at(iparm)->setVal(fitparmsinit.at(iparm));

              RooAbsPdf *rightpdf = 0;
              RooAbsPdf *wrongpdf = 0;
              RooAbsPdf *allpdf = 0;
              RooAbsPdf *interpdf = 0;

              Bool_t usecomplexmodel = complexset.count(catnames.at(icat));
              Bool_t usesimplemodel = simpleset.count(catnames.at(icat));

              if (usecomplexmodel)
              {
                rightpdf = &combhmin;
                wrongpdf = &combhwrong;
                allpdf = &combhvtx;
                interpdf = combhvtxslides[icat][sigProcName];
                f1.setRange(0.52, 1.0);
              }
              else if (usesimplemodel)
              {
                rightpdf = &g1;
                wrongpdf = &wg1;
                allpdf = &combhvtxsimple;
                interpdf = combhvtxsimpleslides[icat][sigProcName];
              }
              else
              {
                rightpdf = &combhmin;
                wrongpdf = &wg1;
                allpdf = &combhvtxmin;
                interpdf = combhvtxminslides[icat][sigProcName];
                sigma1.setVal(2.0);
                sigma2.setVal(3.0);
                f1.setRange(0.0, 1.0);
              }

              mnom->setVal(115);
              hmass->setRange("higgsrange", 100.0, config.massmax);
              hmass->setRange("plotrange", 100.0, 115 + 20.0);

              //----------------------------------------
              // another fit is happening here !!!
              //----------------------------------------
#ifdef TEST_SINGLE_MASS_POINT
              RooFitResult *fitres = rightpdf->fitTo(*testdsets.at(icat),
                  Strategy(1), Minimizer("Minuit2", ""), Minos(kFALSE),
                  Range("higgsrange"), SumW2Error(kTRUE), Save(kTRUE), NumCPU(8));
              RooFitResult *fitreswrong = NULL;
              if (useRightWrongVertex)
              {
                ASSERT(testdsetswrong.at(icat) != NULL);
                wrongpdf->fitTo(*testdsetswrong.at(icat), Strategy(1),
                    Minimizer("Minuit2", ""), Minos(kFALSE), Range("higgsrange"),
                    SumW2Error(kTRUE), Save(kTRUE), NumCPU(8));
              }
              //----------------------------------------

              double fright = 0;
              if (useRightWrongVertex)
              {
                ASSERT(testdsetswrong.at(icat) != NULL);
                fright = testdsets.at(icat)->sumEntries()
                / (testdsets.at(icat)->sumEntries()
                    + testdsetswrong.at(icat)->sumEntries());
              }
              else
              fright = 1;

              fracright.setVal(fright);

              TCanvas *ccompint = new TCanvas;
              TString plotname = TString("inttest") + catnames.at(icat) + "_" + sigProcName + TString(".png");
              RooPlot *hplotcompint = hmass->frame(Bins(100), Range("plotrange"));
              testdsetsall.at(icat)->plotOn(hplotcompint);
              allpdf->plotOn(hplotcompint, RooFit::LineColor(kBlue),
                  NormRange("higgsrange"), Range("higgsrange"));
              interpdf->plotOn(hplotcompint, RooFit::LineColor(kRed),
                  RooFit::LineStyle(kDashed), NormRange("higgsrange"),
                  Range("higgsrange"));
              hplotcompint->SetTitle("");
              hplotcompint->Draw();

              TLegend *legmc = new TLegend(0.62, 0.75, 0.92, 0.9);
              legmc->AddEntry(hplotcompint->getObject(0), "MC", "LPE");
              legmc->AddEntry(hplotcompint->getObject(1), "Fit", "L");
              legmc->AddEntry(hplotcompint->getObject(2), "Interpolated Fit", "L");
              legmc->SetBorderSize(0);
              legmc->SetFillStyle(0);
              legmc->Draw();

              ccompint->SaveAs(plotname);
              TCanvas *csmear = new TCanvas;
              TString smearplotname = TString("smeartest") + catnames.at(icat) + "_"
              + sigProcName + TString(".png");
              RooPlot *hplotsmear = hmass->frame(Bins(100), Range("plotrange"));
              testdsetsall.at(icat)->plotOn(hplotsmear);
              interpdf->plotOn(hplotsmear, RooFit::LineColor(kBlue),
                  NormRange("higgsrange"), Range("higgsrange"));
              nuissancedeltasmears[icat]->setVal(smearingv[icat]);
              //nuissancedeltasmears[icat]->setVal(0.0);
              interpdf->plotOn(hplotsmear, RooFit::LineColor(kRed),
                  NormRange("higgsrange"), Range("higgsrange"));
              hplotsmear->Draw();

              TLegend *legsmear = new TLegend(0.62, 0.75, 0.92, 0.9);
              legsmear->AddEntry(hplotsmear->getObject(0), "MC", "LPE");
              legsmear->AddEntry(hplotsmear->getObject(1), "Fit", "L");
              legsmear->AddEntry(hplotsmear->getObject(2), "Smeared Fit", "L");
              legsmear->SetBorderSize(0);
              legsmear->SetFillStyle(0);
              legsmear->Draw();
              csmear->SaveAs(smearplotname);
#endif

              nuissancedeltasmears[icat]->setVal(0.0);
            } // loop over all signal processes
    } // loop over categories

    // integrated luminosity (that's actually just the scaling parameter with respect to the actual luminosity
    // set in the input data cards for combine ?!)
    RooRealVar intlumi("intlumi", "intlumi", 1.0, 0.0, 100e3);
    intlumi.setConstant();

    // total higgs cross-section -- why is this one and not changed ?!
    RooRealVar *totalxsec = new RooRealVar("totalxsec", "", 1.0, 0.0, 100.0);
    totalxsec->setConstant();

    //--------------------

    // a 'histogram function' containing the (interpolated) cross sections
    map<string, RooHistFunc *> funcXsecNorm;

    // TODO: for the moment, we loop over the input signal processes, in the future we
    // could sum the input contributions of the output groups
    // (currently, these are not needed)
    BOOST_FOREACH(string sigProcName, inputSigProcessNames)
          //for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
          {
            funcXsecNorm[sigProcName] = makeCrossSectionFunction(sigProcName);
          } // loop over signal diagrams
    funcXsecNorm["sum"] = makeCrossSectionFunction("sum");
    //--------------------

    RooRealVar nuissancedeltaeffaccbarrel("CMS_hgg_nuissancedeltaeffaccbarrel", "", 1.0, 0.1, 10.0);
    nuissancedeltaeffaccbarrel.setConstant();
    RooRealVar nuissancedeltaeffaccmixed("CMS_hgg_nuissancedeltaeffaccmixed", "", 1.0, 0.1, 10.0);
    nuissancedeltaeffaccmixed.setConstant();

    RooRealVar nuissancedeltar9fracbarrel("CMS_hgg_nuissancedeltar9fracbarrel", "", 1.0, 0.1, 10.0);
    nuissancedeltar9fracbarrel.setConstant();
    RooRealVar nuissancedeltar9fracmixed("CMS_hgg_nuissancedeltar9fracmixed", "", 1.0, 0.1, 10.0);
    nuissancedeltar9fracmixed.setConstant();

    RooRealVar nuissancehighptfrac("CMS_hgg_nuissancehighptfrac", "", 1.0, 0.1, 10.0);
    nuissancehighptfrac.setConstant();

    RooFormulaVar *effaccbarrel = 0;
    RooFormulaVar *effaccmixed = 0;
    RooFormulaVar *r9fracbarrel = 0;
    RooFormulaVar *r9fracmixed = 0;
    RooAbsReal *ptfracbarrelhighr9 = 0;
    RooAbsReal *ptfracbarrelmixedr9 = 0;
    RooAbsReal *ptfracmixedhighr9 = 0;
    RooAbsReal *ptfracmixedmixedr9 = 0;

    RooAbsReal *ptfrac_B_barrelhighr9 = 0;
    RooAbsReal *ptfrac_B_barrelmixedr9 = 0;
    RooAbsReal *ptfrac_B_mixedhighr9 = 0;
    RooAbsReal *ptfrac_B_mixedmixedr9 = 0;

    // these (for the moment) do NOT depend on the signal production mechanism
    // sum over all production mechanism normalizations

    // Note that this adds up all finalnormslides functions
    // 4 categories: cats 0 and 1 are barrel
    // 8 categories: cats 0,1,4,5
    effaccbarrel = makeEffAccVar("effaccbarrel", finalnormslides, nuissancedeltaeffaccbarrel, getBarrelCategories());

    // 4 categories: cats 2 and 3 are mixed
    // 8 categories: cats 2,3,6,7 are mixed
    effaccmixed = makeEffAccVar("effaccmixed", finalnormslides, nuissancedeltaeffaccmixed, getEndcapCategories());

    if (config.numInclusiveCategories == 4)
    {
      // cat 0 and 1 (barrel)
      // cat 0 [BB, high R9] / (cat 0 [BB, high R9] + cat 1 [BB, low R9])
      r9fracbarrel = makeFracVar("r9fracbarrel", finalnormslides, nuissancedeltar9fracbarrel, list_of(0), list_of(1));

      r9fracmixed = makeFracVar("r9fracmixed", finalnormslides, nuissancedeltar9fracmixed, list_of(2), list_of(3));
      //      r9fracmixed = new RooFormulaVar("r9fracmixed", "", "@0*(@1+@2+@3+@4)/(@1+@2+@3+@4 + @5+@6+@7+@8)", RooArgList(nuissancedeltar9fracmixed,
      //          // cat 2
      //          *finalnormslides[2][sigProcessNames[0]], *finalnormslides[2][sigProcessNames[1]], *finalnormslides[2][sigProcessNames[2]], *finalnormslides[2][sigProcessNames[3]],
      //
      //          // cat 3
      //          *finalnormslides[3][sigProcessNames[0]], *finalnormslides[3][sigProcessNames[1]], *finalnormslides[3][sigProcessNames[2]], *finalnormslides[3][sigProcessNames[3]]));

      // there is no split by pT here
      ptfracbarrelhighr9 = new RooConstVar("ptfracbarrelhighr9", "", 1.0);
      ptfracbarrelmixedr9 = new RooConstVar("ptfracbarrelmixedr9", "", 1.0);
      ptfracmixedhighr9 = new RooConstVar("ptfracmixedhighr9", "", 1.0);
      ptfracmixedmixedr9 = new RooConstVar("ptfracmixedmixedr9", "", 1.0);

      ptfrac_B_barrelhighr9 = new RooConstVar("ptfrac_B_barrelhighr9", "", 0.0);
      ptfrac_B_barrelmixedr9 = new RooConstVar("ptfrac_B_barrelmixedr9", "", 0.0);
      ptfrac_B_mixedhighr9 = new RooConstVar("ptfrac_B_mixedhighr9", "", 0.0);
      ptfrac_B_mixedmixedr9 = new RooConstVar("ptfrac_B_mixedmixedr9", "", 0.0);

    }
    else if (config.numInclusiveCategories == 8)
    {
      // build the quantities we want to use for systematic uncertainties from the per-category efficiencies
      // determined above (a kind of 'base transformation')

      /*
        r9fracbarrel = new RooFormulaVar("r9fracbarrel","","@0 * (@1 + @2)/(@1 +@2 +@3 + @4)",RooArgList(nuissancedeltar9fracbarrel,
              *finalnormslides[0],
              *finalnormslides[4],
              *finalnormslides[1],
              *finalnormslides[5]) );
      */
      r9fracbarrel = makeFracVar("r9fracbarrel", finalnormslides, nuissancedeltar9fracbarrel, list_of(0)(4), list_of(1)(5));

      /*      r9fracmixed  = new RooFormulaVar("r9fracmixed", "","@0 * (@1 + @2)/(@1 +@2 +@3 + @4)",RooArgList(nuissancedeltar9fracmixed,
                  *finalnormslides[2],
                  *finalnormslides[6],
                  *finalnormslides[3],
                  *finalnormslides[7]));
      */
      r9fracmixed = makeFracVar("r9fracmixed", finalnormslides, nuissancedeltar9fracmixed, list_of(2)(6), list_of(3)(7));

      // pt fraction variables (separate per (Barrel/Endcap categories) x (R9 categories) )
      /*
      ptfrac_barrel_highr9  = new RooFormulaVar("ptfracbarrelhighr9","","@0*@1/(@1+@2)", RooArgList(nuissancehighptfrac,
                                         *finalnormslides[0],*finalnormslides[4]));
      */
      //                                                                                        BB high R9 high pT
      //                                                                                                       BB high R9 low pT
      ptfracbarrelhighr9 = makeFracVar("ptfracbarrelhighr9", finalnormslides, nuissancehighptfrac, list_of(0), list_of(4));

      /*      ptfrac_barrel_mixedr9 = new RooFormulaVar("ptfracbarrelmixedr9","","@0*@1/(@1+@2)",RooArgList(nuissancehighptfrac,
                                      *finalnormslides[1],*finalnormslides[5])); */
      //                                                                                        BB low R9 high pT
      //                                                                                                       BB low R9 low pT
      ptfracbarrelmixedr9 = makeFracVar("ptfracbarrelmixedr9", finalnormslides, nuissancehighptfrac, list_of(1), list_of(5));

      /*      ptfrac_mixed_highr9   = new RooFormulaVar("ptfracmixedhighr9","","@0*@1/(@1+@2)",  RooArgList(nuissancehighptfrac,
                   *finalnormslides[2],*finalnormslides[6])); */
      //                                                                                        BE high R9 high pT
      //                                                                                                       BE high R9 low pT
      ptfracmixedhighr9 = makeFracVar("ptfracmixedhighr9", finalnormslides, nuissancehighptfrac, list_of(2), list_of(6));

      /*      ptfrac_mixed_mixedr9  = new RooFormulaVar("ptfracmixedmixedr9","","@0*@1/(@1+@2)", RooArgList(nuissancehighptfrac,
               *finalnormslides[3],*finalnormslides[7])); */
      //                                                                                        BE low R9 high pT
      //                                                                                                       BE low R9 low pT
      ptfracmixedmixedr9 = makeFracVar("ptfracmixedmixedr9", finalnormslides, nuissancehighptfrac, list_of(3), list_of(7));

      ptfrac_B_barrelhighr9 = new RooConstVar("ptfrac_B_barrelhighr9", "", 0.0);
      ptfrac_B_barrelmixedr9 = new RooConstVar("ptfrac_B_barrelmixedr9", "", 0.0);
      ptfrac_B_mixedhighr9 = new RooConstVar("ptfrac_B_mixedhighr9", "", 0.0);
      ptfrac_B_mixedmixedr9 = new RooConstVar("ptfrac_B_mixedmixedr9", "", 0.0);

    }
    else if (config.numInclusiveCategories == 12)
    {

      // for 12, one would has to introduce two pt fractions (with the same nuisance parameter ?)
      // pt frac A for cats 0..3
      // pt frac B for cats 4..7
      // pt frac for cats 8..11 is 1 - (pt frac A + pt frac B)

      r9fracbarrel = makeFracVar("r9fracbarrel", finalnormslides, nuissancedeltar9fracbarrel, list_of(0)(4)(8), list_of(1)(5)(9));

      /*      r9fracmixed  = new RooFormulaVar("r9fracmixed", "","@0 * (@1 + @2)/(@1 +@2 +@3 + @4)",RooArgList(nuissancedeltar9fracmixed,
                  *finalnormslides[2],
                  *finalnormslides[6],
                  *finalnormslides[3],
                  *finalnormslides[7]));
      */
      r9fracmixed = makeFracVar("r9fracmixed", finalnormslides, nuissancedeltar9fracmixed, list_of(2)(6)(10), list_of(3)(7)(11));

#warning CHECK THIS IMPLEMENTATION OF THE PT FRACTIONS HERE

      //----------
      // pt fraction variables (separate per (Barrel/Endcap categories) x (R9 categories) )
      /*
      ptfrac_barrel_highr9  = new RooFormulaVar("ptfracbarrelhighr9","","@0*@1/(@1+@2)", RooArgList(nuissancehighptfrac,
                                         *finalnormslides[0],*finalnormslides[4]));
      */
      //                                                                                        BB high R9 high pT
      //                                                                                                       BB high R9 low pT
      //                                                                                                       and highest pT
      ptfracbarrelhighr9 = makeFracVar("ptfracbarrelhighr9", finalnormslides, nuissancehighptfrac, list_of(0), list_of(4)(8));

      /*      ptfrac_barrel_mixedr9 = new RooFormulaVar("ptfracbarrelmixedr9","","@0*@1/(@1+@2)",RooArgList(nuissancehighptfrac,
                                      *finalnormslides[1],*finalnormslides[5])); */
      //                                                                                        BB low R9 high pT
      //                                                                                                       BB low R9 low pT
      //                                                                                                       and highest pT
      ptfracbarrelmixedr9 = makeFracVar("ptfracbarrelmixedr9", finalnormslides, nuissancehighptfrac, list_of(1), list_of(5)(9));

      /*      ptfrac_mixed_highr9   = new RooFormulaVar("ptfracmixedhighr9","","@0*@1/(@1+@2)",  RooArgList(nuissancehighptfrac,
                   *finalnormslides[2],*finalnormslides[6])); */
      //                                                                                        BE high R9 high pT
      //                                                                                                       BE high R9 low pT
      //                                                                                                       and highest pT
      ptfracmixedhighr9 = makeFracVar("ptfracmixedhighr9", finalnormslides, nuissancehighptfrac, list_of(2), list_of(6)(10));

      /*      ptfrac_mixed_mixedr9  = new RooFormulaVar("ptfracmixedmixedr9","","@0*@1/(@1+@2)", RooArgList(nuissancehighptfrac,
               *finalnormslides[3],*finalnormslides[7])); */
      //                                                                                        BE low R9 high pT
      //                                                                                                       BE low R9 low pT
      //                                                                                                       and highest pT
      ptfracmixedmixedr9 = makeFracVar("ptfracmixedmixedr9", finalnormslides, nuissancehighptfrac, list_of(3), list_of(7)(11));

      //--------------------
      // do the same pT fractions for categories 8..11
      //--------------------
      ptfrac_B_barrelhighr9 = makeFracVar("ptfrac_B_barrelhighr9", finalnormslides, nuissancehighptfrac, list_of(8), list_of(0)(4));

      ptfrac_B_barrelmixedr9 = makeFracVar("ptfrac_B_barrelmixedr9", finalnormslides, nuissancehighptfrac, list_of(9), list_of(1)(5));

      ptfrac_B_mixedhighr9 = makeFracVar("ptfrac_B_mixedhighr9", finalnormslides, nuissancehighptfrac, list_of(10), list_of(2)(6));

      ptfrac_B_mixedmixedr9 = makeFracVar("ptfrac_B_mixedmixedr9", finalnormslides, nuissancehighptfrac, list_of(11), list_of(3)(7));

    }
    else
    {
      cerr << "numInclusiveCategories = " << config.numInclusiveCategories << " is not yet supported here" << endl;
      exit(1);
    }

    map<unsigned, map<string, RooAbsReal*> > nsigcats;

    // define RooCategory for simultaneous fit
    RooCategory finalcat("finalcat", "finalcat");
    for (unsigned j = 0; j < catnames.size(); ++j)
      finalcat.defineType(Form("cat%d", j));

    // one pdf to rule them all
    RooSimultaneous fullsigpdf("fullsigpdf", "fullsigpdf", finalcat);
    RooSimultaneous fullsigpdfrel("fullsigpdfrel", "fullsigpdfrel", finalcat);

    vector<RooAbsArg *> pdfsToAdd;

    //----------
    map<unsigned, map<string, RooFormulaVar *> > nsigrelcat;

    // produce quantities which are per unmerged signal process
    BOOST_FOREACH(string sigProcName, inputSigProcessNames)
          {
            map<unsigned, RooFormulaVar *> nsigcat;

            for (unsigned cat = 0; cat < config.numInclusiveCategories; ++cat)
            {
              RooFormulaVar *tmp = makeNsigCat(cat, sigProcName, intlumi, totalxsec,

              effaccbarrel, effaccmixed, r9fracbarrel, r9fracmixed,

                // pt fraction variables
                ptfracbarrelhighr9,    ptfracbarrelmixedr9,    ptfracmixedhighr9,    ptfracmixedmixedr9,
                ptfrac_B_barrelhighr9, ptfrac_B_barrelmixedr9, ptfrac_B_mixedhighr9, ptfrac_B_mixedmixedr9,

                fitparmfuncsX, effaccFuncsX);
              nsigcat[cat] = tmp;
              nsigcats[cat][sigProcName] = tmp;
              pdfsToAdd.push_back(nsigcat[cat]);

              //                                                 xsecnorm
              //                                                     nsigcatX
              // TODO: clarify why using the 'sum' cross section (instead of sigProcName) here gives the right values
              nsigrelcat[cat][sigProcName] = new RooFormulaVar(Form("hggpdfrel_cat%d_%s_norm", cat, sigProcName.c_str()), "", "@0*@1", RooArgList(*funcXsecNorm["sum"], *nsigcat[cat]));

            } // loop over categories
          } // loop over unmerged signal process names

    //----------------------------------------
    // create the PDFs
    //----------------------------------------
    BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
          // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
          {

            map<unsigned, RooExtendPdf *> sigpdfcat;
            map<unsigned, RooAbsPdf *> hggpdfrel_cat;
            map<unsigned, RooExtendPdf *> sigpdfrelcat;

            string mergedSigProcName = config.signalProcessMergingMapping[inputSigProcName];

            for (unsigned cat = 0; cat < config.numInclusiveCategories; ++cat)
            {
              // final extended pdf's in each category
              //                                                                                            normalized pdf     normalization
              sigpdfcat[cat] = new RooExtendPdf(Form("sigpdfcat%d_%s", cat, inputSigProcName.c_str()), "", *finalpdfslides[cat][inputSigProcName], *nsigcats[cat][inputSigProcName]);

              // relative cross section version of pdfs and normalization
              hggpdfrel_cat[cat] = (RooAbsPdf*) finalpdfslides[cat][inputSigProcName]->Clone(Form("hggpdfrel_cat%d_%s", cat, inputSigProcName.c_str()));
              pdfsToAdd.push_back(hggpdfrel_cat[cat]);

              sigpdfrelcat[cat] = new RooExtendPdf(Form("sigpdfrelcat%d_%s", cat, inputSigProcName.c_str()), "", *hggpdfrel_cat[cat], *nsigrelcat[cat][inputSigProcName]);
              pdfsToAdd.push_back(sigpdfrelcat[cat]);

              // does this work: adding to each category more than once ? No, it does NOT (gives RooFit error messages).
              fullsigpdf.addPdf(*sigpdfcat[cat], Form("cat%d", cat));
              pdfsToAdd.push_back(sigpdfcat[cat]);

              fullsigpdfrel.addPdf(*sigpdfrelcat[cat], Form("cat%d", cat));

            } // loop over categories
          } // loop over signal processes

    //----------------------------------------
    // additional UCSD (VBF/VHAD) categories
    //----------------------------------------

    map<unsigned, map<string, RooFormulaVar *> > additionalUCSDnsigcat;

    // BOOST_FOREACH(string sigProcName, mergedSigProcessNames)
    BOOST_FOREACH(string inputSigProcName, inputSigProcessNames)
          // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
          {
            string mergedSigProcName = config.signalProcessMergingMapping[inputSigProcName];

            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
            {

              // create a RooHistFunc from a RooDataHist from the TH1F filled with the total
              // efficiencies in category 'numInclusiveCategories'
              RooDataHist *dataHist = new RooDataHist(Form("effaccCat%dDataHist_%s", j, inputSigProcName.c_str()), Form("effaccCat%dDataHist_%s", j, inputSigProcName.c_str()), RooArgList(*mnom), additionalUCSDcatEffaccHist[j][inputSigProcName]);
              RooHistFunc *histFunc = new RooHistFunc(Form("effaccCat%dHistFunc_%s", j, inputSigProcName.c_str()), Form("effaccCat%dHistFunc_%s", j, inputSigProcName.c_str()), RooArgList(*mnom), *dataHist, 1);

              // TODO: CHECK WHETHER THIS GIVES THE CORRECT NORMALIZATION
              // request 2011-12-19: use the merged PDF but use per unmerged pdf normalization factors
              additionalUCSDnsigcat[j][inputSigProcName] = new RooFormulaVar(Form("hggpdf_cat%d_%s_norm", j, inputSigProcName.c_str()), "", "@0*@1*@2", RooArgList(intlumi, *totalxsec, *histFunc));
              pdfsToAdd.push_back(additionalUCSDnsigcat[j][inputSigProcName]);
            } // loop over categories

            map<unsigned, map<string, RooExtendPdf*> > additionalUCSDsigpdfcat;
            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
              additionalUCSDsigpdfcat[j][inputSigProcName] = new RooExtendPdf(Form("sigpdfcat%d_%s", j, inputSigProcName.c_str()), "", *finalpdfslides[0][inputSigProcName], *additionalUCSDnsigcat[j][inputSigProcName]);

            map<unsigned, map<string, RooAbsPdf*> > additionalUCSDhggpdfrel;
            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
            {
              additionalUCSDhggpdfrel[j][inputSigProcName] = ((RooAbsPdf*) finalpdfslides[j][inputSigProcName]->Clone(Form("hggpdfrel_cat%d_%s", j, inputSigProcName.c_str())));
              pdfsToAdd.push_back(additionalUCSDhggpdfrel[j][inputSigProcName]);
            }

            map<unsigned, map<string, RooFormulaVar*> > additionalUCSDnsigrel;
            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
              // TODO: clarify why using the 'sum' cross section (instead of sigProcName) here gives the right values
              additionalUCSDnsigrel[j][inputSigProcName] = new RooFormulaVar(Form("hggpdfrel_cat%d_%s_norm", j, inputSigProcName.c_str()), "", "@0*@1", RooArgList(*funcXsecNorm["sum"], *additionalUCSDnsigcat[j][inputSigProcName]));

            map<unsigned, map<string, RooExtendPdf*> > additionalUCSDsigpdfrel;
            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
            {
              additionalUCSDsigpdfrel[j][inputSigProcName] = new RooExtendPdf(Form("sigpdfrelcat%d_%s", j, inputSigProcName.c_str()), "", *additionalUCSDhggpdfrel[j][inputSigProcName], *additionalUCSDnsigrel[j][inputSigProcName]);
              pdfsToAdd.push_back(additionalUCSDsigpdfrel[j][inputSigProcName]);
            }

            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
            {
              fullsigpdf.addPdf(*additionalUCSDsigpdfcat[j][inputSigProcName], Form("cat%d", j));
              pdfsToAdd.push_back(additionalUCSDsigpdfcat[j][inputSigProcName]);
            }

            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
              fullsigpdfrel.addPdf(*additionalUCSDsigpdfrel[j][inputSigProcName], Form("cat%d_%s", j, inputSigProcName.c_str()));

            for (unsigned j = config.numInclusiveCategories; j < catnames.size(); ++j)
              nsigcats[j][inputSigProcName] = additionalUCSDnsigcat[j][inputSigProcName];

          } // loop over signal processes

    //----------------------------------------
    //   double nsigtotal = 0.0;
    //   for (int i=0; i<nsigcats.size(); ++i) {
    //     nsigtotal += nsigcats.at(i)->getVal();
    //   }

    // the overall normalization sum (?!)
    RooFormulaVar *hggpdf_combcat_norm = 0;

    RooArgList addnorm;
    {
      std::string buf;
      bool isFirst = true;

      unsigned k = 0;
      assert(catnames.size() == nsigcats.size());

      for (unsigned cat = 0; cat < catnames.size(); ++cat)
      {
        BOOST_FOREACH(string sigProcName, mergedSigProcessNames)
              // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
              {
                // compose the overall sum expression
                if (isFirst)
                  isFirst = false;
                else
                  buf += "+";

                buf += Form("@%d", k++);

                addnorm.add(*nsigcats[cat][sigProcName]);

              } // loop over signal processes
      } // loop over categories

      hggpdf_combcat_norm = new RooFormulaVar("hggpdf_combcat_norm", "", buf.c_str(), addnorm);

    } // end of block
    //----------

#if 0
    RooArgList addpdfs;
    RooArgList addcoeffs;
    cout << "nsigcats.size()=" << nsigcats.size() << endl;
    {
      int k = 0;
      for (int i = 0; i < nsigcats.size(); ++i)
      {
        for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
        {
          string sigProcName = sigProcessNames.at(sigProcessIndex);

          addpdfs.add(*finalpdfslides[i][sigProcName]);

#error not clear how this should generalize to the case with per production mechanism quantities
          if (k < (nsigcats.size() - 1))
          {
            // divides by the overall norm calculated just above
            RooFormulaVar *addcoef = new RooFormulaVar(
                TString("addcoef_" + catnames.at(i)), "", "@0/@1",
                RooArgList(*nsigcats.at(i), *hggpdf_combcat_norm));
            addcoeffs.add(*addcoef);
          }

          ++k;
        } // loop over signal processes

      } // loop over categories

      cout << "addpfs.size()=" << addpdfs.getSize() << " addcoeffs.size()="
      << addcoeffs.getSize() << endl;
    } // end of block
    RooAddPdf hggpdf_combcat("hggpdf_combcat","",addpdfs,addcoeffs);
    RooExtendPdf sigpdfcombcat("sigpdfcombcat","",hggpdf_combcat,*hggpdf_combcat_norm);
#endif
    //----------

    //   RooAbsReal *hggpdf_singlecat_norm = 0;
    //   if (catnames.at(catnames.size()-1)=="singlecat") {
    //     hggpdf_singlecat_norm = finalnormslides[catnames.size()-1];
    //   }
    //   hggpdf_singlecat_norm->SetName("hggpdf_singlecat_norm");
    //   RooExtendPdf sigpdfsinglecat("sigpdfsinglecat","",*finalpdfslides[catnames.size()-1],*hggpdf_singlecat_norm);

    // set the default Higgs mass hypothesis
    mnom->setVal(130);

    //save everything to file with RooWorkspace
    RooWorkspace *w = new RooWorkspace("wsig", "");

#warning commented out for the moment
    // w->import(fullsigpdf,RecycleConflictNodes());
    // w->import(fullsigpdfrel,RecycleConflictNodes());

    for (unsigned k = 0; k < pdfsToAdd.size(); ++k)
      w->import(*pdfsToAdd[k], RecycleConflictNodes());

    //w->import(sigpdfsinglecat,RecycleConflictNodes());

#if 0
    w->import(sigpdfcombcat,RecycleConflictNodes());
#endif
    //w->import(*combhvtxslides[0]);
    for (unsigned icat = 0; icat < catnames.size(); ++icat)
    {
      BOOST_FOREACH(string sigProcName, inputSigProcessNames)
            {
              // funceffacc_catX
              w->import(*effaccFuncsX[icat][sigProcName], RecycleConflictNodes());
            }
    }

    w->Print();

    w->writeToFile("ubersignalmodel.root");
    wextra->writeToFile("extra.root");

    //  if (catnames.size()>=8)
    BOOST_FOREACH(string sigProcName, inputSigProcessNames)
          // for (unsigned sigProcessIndex = 0; sigProcessIndex < sigProcessNames.size(); ++sigProcessIndex)
          {
            for (unsigned j = 0; j < catnames.size(); ++j)
              printf("cat%d proc %s: finalnormslides = %5f, n = %5f\n", j, sigProcName.c_str(), finalnormslides[j][sigProcName]->getVal(), nsigcats[j][sigProcName]->getVal());
          }

    return;

  }

  //----------------------------------------------------------------------
};

//----------------------------------------------------------------------
int
main(int argc, char **argv)
{
  if (argc != 2)
  {
    cerr << "usage: makeParametricSignalModel config.xml" << endl;
    exit(1);
  }

  ParametricSignalModelConfig config = ParametricSignalModelConfig::read(argv[1]);

  gROOT->SetStyle("Plain");
  FitHmassGlobe4CatClass fitter(config);
  fitter.run();

  cout << "----------------------------------------" << endl;
  cout << "signal fitting done" << endl;
  cout << "----------------------------------------" << endl;

  gROOT->cd();
}
//----------------------------------------------------------------------

