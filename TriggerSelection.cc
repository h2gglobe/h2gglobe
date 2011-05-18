#include "TriggerSelection.h"

#include <iostream>
#include <algorithm>

using namespace std;

// ----------------------------------------------------------------------------------------------------
class IsSubstring : public binary_function<std::string, std::string, bool>
{
public:
	bool operator () (const std::string & lh, const std::string & rh) const { 
		return lh.find(rh) != std::string::npos; 
	}
};

// ----------------------------------------------------------------------------------------------------
bool TriggerSelection::pass(const std::vector<std::string> & menu, const std::vector<unsigned short> & bits) 
{
	// loop over requestedpath names
	for(std::vector<std::string>::iterator it=paths.begin(); it!=paths.end(); ++it ) {
		// is the path in the menu?
		std::vector<std::string>::const_iterator jt=find_if(menu.begin(),menu.end(), bind2nd(IsSubstring(),*it) );
		if( jt != menu.end() ) {
			// if yes check if the bit fired
			int ibit = jt - menu.begin();
			if( find(bits.begin(),bits.end(),ibit) != bits.end() ) { 
				return true; 
			}
		}
	}
	return false;
}
