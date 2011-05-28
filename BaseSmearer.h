#ifndef __BASESMEARER__
#define __BASESMEARER__

#include "LoopAll.h"

// ------------------------------------------------------------------------------------
/**
 * \class BaseSmearer
 *
 * Defines the interface to be implemented by any smearing function on photon informations
 *
 *
 */

class PhotonReducedInfo;

class BaseSmearer 
{
public:
	// ! C-TOR
	BaseSmearer();
	
	// ! D-TOR
	virtual ~BaseSmearer();
	
	// ! the class name
	virtual const std::string & name() const = 0; 
	
	// !  
	operator const std::string & () const { return this->name(); };
	
	// ! return smeared photon informations
	virtual bool smearPhoton( PhotonReducedInfo & pho, float & weight, float syst_shift=0. ) const = 0;
	
};

// ! Used to search analyzers by name 
bool operator == (BaseSmearer * a, const std::string & b);

#endif
