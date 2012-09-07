#ifndef __PARSER__
#define __PARSER__

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>

#include <dlfcn.h>

#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
std::string getString(const map<string, string> &values, const std::string &paramName);
std::string getString(const map<string, string> &values, const std::string &paramName, const std::string &defaultValue);
unsigned getUint(const map<string, string> &values, const std::string &paramName);
unsigned getFloat(const map<string, string> &values, const std::string &paramName);
std::vector<int> getVint(const map<string, string> &values, const std::string &paramName);
map<string, string> parseConfigFile(const std::string &configFname);

#endif
