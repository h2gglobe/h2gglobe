#include "parser.h"

//----------------------------------------------------------------------

std::string getString(const map<string, string> &values, const std::string &paramName) {

  map<string, string>::const_iterator it = values.find(paramName);
  if (it == values.end())
    throw std::runtime_error("no parameter named '" + paramName + "' found");
  
  return it->second;
}

//----------------------------------------------------------------------

std::string  getString(const map<string, string> &values, const std::string &paramName, const std::string &defaultValue) {

  map<string, string>::const_iterator it = values.find(paramName);
  if (it == values.end())
    return defaultValue;
  else
    return it->second;
}

//----------------------------------------------------------------------

unsigned getUint(const map<string, string> &values, const std::string &paramName) {
  try {
    return boost::lexical_cast<unsigned>(getString(values, paramName));
  } catch (const boost::bad_lexical_cast &ex) {
    throw std::runtime_error("value of parameter '" + paramName + "' is not a valid unsigned integer");
  }
}

//----------------------------------------------------------------------

unsigned getFloat(const map<string, string> &values, const std::string &paramName) {
  try {
    return boost::lexical_cast<float>(getString(values, paramName));
  } catch (const boost::bad_lexical_cast &ex) {
    throw std::runtime_error("value of parameter '" + paramName + "' is not a valid float");
  }
}

//----------------------------------------------------------------------

#include <boost/algorithm/string.hpp>

std::vector<int> getVint(const map<string, string> &values, const std::string &paramName) {
  try {
    vector<string> stringValues;
    string buf = getString(values, paramName);
    boost::split(stringValues, buf, boost::is_any_of(","));
    
    vector<int> retval;
    
    BOOST_FOREACH(string str, stringValues) {
      retval.push_back(boost::lexical_cast<int>(str));
    }

    return retval;
  } catch (const boost::bad_lexical_cast &ex) {
    throw std::runtime_error("values of parameter '" + paramName + "' contain invalid integers");
  }
}

//----------------------------------------------------------------------

std::vector<float> getVfloat(const map<string, string> &values, const std::string &paramName) {
  try {
    vector<string> stringValues;
    string buf = getString(values, paramName);
    boost::split(stringValues, buf, boost::is_any_of(","));
    
    vector<float> retval;
    
    BOOST_FOREACH(string str, stringValues) {
      retval.push_back(boost::lexical_cast<float>(str));
    }

    return retval;
  } catch (const boost::bad_lexical_cast &ex) {
    throw std::runtime_error("values of parameter '" + paramName + "' contain invalid integers");
  }
}

//----------------------------------------------------------------------

map<string, string> parseConfigFile(const std::string &configFname) {

  map<string, string> values;
  
  ifstream infile(configFname.c_str());
  string line;

  int lineNum = 0;

  while ( infile.good() ) {
    try {
      getline (infile,line);
      ++lineNum;
      //cout << line << endl;
      
      // remove everything after the first #
      size_t pos = line.find('#');
      if (pos != string::npos)
	line.erase(pos);
      
      // remove leading a trailing white space
      boost::algorithm::trim(line);
      
      if (line == "")
	continue;
      
      // split line into name=value pairs, separated by whitespace
      vector<string> parts;
      
      boost::algorithm::split_regex(parts, line,
				    boost::regex( "\\s+"));
      
      BOOST_FOREACH(std::string part, parts) {
	// loop over all parts of this line
	pos = part.find('=');
	assert(pos != string::npos);
	
	string key = part.substr(0,pos);
	string value = part.substr(pos+1);
	values[key] = value;
      }
      
      //if (bookHist)
	//bookHistogram(h, values);
      //histogramDefinitions.push_back(ConfigLineData(lineNum, values));
    }  catch (const std::exception &ex) {
      cerr << "exception caught while reading line " << lineNum << " of file " << configFname << ": " << ex.what() << endl;
      exit(1);
    }
  } // loop over lines

  return values;
}
