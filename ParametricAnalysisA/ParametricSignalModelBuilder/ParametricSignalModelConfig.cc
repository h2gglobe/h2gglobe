#include "ParametricSignalModelConfig.h"
#include <TDOMParser.h>
#include "XMLUtils.h"

#include "utils.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <cstdlib>



using namespace std;

//----------------------------------------------------------------------
ParametricSignalModelConfig::ParametricSignalModelConfig() :
  fermiophobic(false), numCategories(0), numInclusiveCategories(0), useRightWrongVertex(true), massmax(-1), massInterpolationWidth(-1), weightscale(2e3),
  sqrts(0)
{
  // initialize a default signal process mapping
  initializeSignalProcessMergingMapping(fermiophobic);
}

//----------------------------------------------------------------------
ParametricSignalModelConfig
ParametricSignalModelConfig::read(const std::string &configFname)
{
  // parse the XML file
  TDOMParser parser;

  // don't compare to a DTD for the moment
  parser.SetValidate(false);

  int res = parser.ParseFile(configFname.c_str());
  if (res != 0)
  {
    cerr << "parsing configuration file " << configFname << " failed with error " << res << ", exiting" << endl;
    exit(1);
  }

  TXMLDocument *doc = parser.GetXMLDocument();
  TXMLNode *root = doc->GetRootNode();

  ParametricSignalModelConfig config;

  //--------------------
  try
  {
    config.inputFname = XMLUtils::getAttribute(XMLUtils::getSingleChild(root, "inputFile"), "name");

    config.inputWorkspaceName = XMLUtils::getAttribute(XMLUtils::getSingleChild(root, "inputFile"), "workspace", "cms_hgg_workspace");

    config.fermiophobic = XMLUtils::getBooleanContent(XMLUtils::getSingleChild(root, "fermiophobic"));

    config.initializeSignalProcessMergingMapping(config.fermiophobic);

    config.numCategories = XMLUtils::getIntegerContent(XMLUtils::getSingleChild(root, "numCategories"));

    config.numInclusiveCategories = XMLUtils::getIntegerContent(XMLUtils::getSingleChild(root, "numInclusiveCategories"));

    config.useRightWrongVertex = XMLUtils::getBooleanContent(root, "useRightWrongVertex", true);

    config.massmax = XMLUtils::getDoubleContent(root, "maxMassHypothesis", 160);

    config.massInterpolationWidth = XMLUtils::getDoubleContent(root, "massInterpolationWidth", 10);

    config.sqrts = XMLUtils::getIntegerContent(XMLUtils::getSingleChild(root, "sqrts"));

    config.nameSuffix = XMLUtils::getTextContent(XMLUtils::getSingleChild(root, "nameSuffix"));

    //----------
    // read the signal process mapping (which input signal processes should be merged
    // into what output process names for fitting the signal shape)
    TXMLNode *signalProcessMerging = XMLUtils::getOptionalSingleChild(root, "signalProcessMerging");
    if (signalProcessMerging != NULL)
    {
      config.signalProcessMergingMapping.clear();

      vector<TXMLNode *> children = XMLUtils::getChildren(signalProcessMerging, "map");

      BOOST_FOREACH(TXMLNode *child, children)
            {
              string src = XMLUtils::getAttribute(child, "from");
              string dest = XMLUtils::getAttribute(child, "to");

              // make sure that 'from' is not already existing in the map
              if (config.signalProcessMergingMapping.count(src) > 0)
              {
                cerr << "error in signalProcessMerging specification: source '" << src << "' appears more than once" << endl;
                exit(1);
              }
              config.signalProcessMergingMapping[src] = dest;

            } // loop over all children (mapping elements)

      // TODO: we should check here whether there was any mapping at
      // all specified

    }

    //----------
    // read signal MC smearing values
    //----------
    TXMLNode *signalMCsmearing = XMLUtils::getSingleChild(root, "signalMCsmearing");
    {
      config.smearingv.clear();

      // first read the values into a map and check later that
      // each category has exactly one value
      map<unsigned, double> tmp;

      vector<TXMLNode *> children = XMLUtils::getChildren(signalMCsmearing, "smearing");

      BOOST_FOREACH(TXMLNode *child, children)
      {
        unsigned cat = boost::lexical_cast<unsigned>(XMLUtils::getAttribute(child, "cat"));
        double value = boost::lexical_cast<double>(XMLUtils::getAttribute(child, "value"));

        // check that the category is within the expected range
        if (cat < 0 || cat >= config.numCategories)
        {
          cerr << "error in signalMCsmearing specification: category '" << cat << "' is out of range" << endl;
          exit(1);
        }

        // make sure that 'from' is not already existing in the map
        if (tmp.count(cat) > 0)
        {
          cerr << "error in signalMCsmearing specification: category '" << cat << "' appears more than once" << endl;
          exit(1);
        }
        tmp[cat] = value;
        //
      } // loop over all children (mapping elements)

      // now fill them in and insist that there is at least one value for each category
      for (unsigned cat = 0; cat < config.numCategories; ++cat)
      {
        if (tmp.count(cat) < 1)
        {
          cerr << "error in signalMCsmearing specification: category '" << cat << "' not specified" << endl;
          exit(1);
        }

        ASSERT(config.smearingv.size() == cat);
        config.smearingv.push_back(tmp[cat]);
      } // loop over all categories

    }




    //----------

  }
  catch (...)
  {
    cerr << "exception caught while reading configuration file " << configFname << endl;
    exit(1);
  }
  //--------------------

  return config;
}

//----------------------------------------------------------------------

void
ParametricSignalModelConfig::initializeSignalProcessMergingMapping(bool fermiophobic)
{
  signalProcessMergingMapping.clear();

  if (false)
  {
    // default mapping
    signalProcessMergingMapping["ggh"] = "ggh";
    signalProcessMergingMapping["tth"] = "tth";
    signalProcessMergingMapping["vbf"] = "vbf";
    signalProcessMergingMapping["wzh"] = "wzh";
  }

  if (true)
  {
    // grouping into ggh and the rest
    if (! fermiophobic)
    {
      signalProcessMergingMapping["ggh"] = "ggh";
      signalProcessMergingMapping["tth"] = "other";
      signalProcessMergingMapping["vbf"] = "other";
      signalProcessMergingMapping["wzh"] = "other";
    }
    else
    {
      // fermiophobic (add wzh to vbf for fitting)
      signalProcessMergingMapping["vbf"] = "vbf";
      signalProcessMergingMapping["wzh"] = "vbf";
    }

  }

}
