#ifndef __TRIGGERSELECTION__ 
#define __TRIGGERSELECTION__  

#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------------------
class TriggerSelection 
{
public:
	TriggerSelection(int first, int last) : firstrun(first), lastrun(last) {};
	bool operator == (int run) { return run>=firstrun && ( lastrun<0 || run<=lastrun); }; 
	void addpath(const std::string & x) { paths.push_back(x); };
	
	bool pass(const std::vector<std::string> & menu, const std::vector<unsigned short> & bits);  
	
	int firstrun,  lastrun;
	std::vector<std::string> paths;
};

#endif
